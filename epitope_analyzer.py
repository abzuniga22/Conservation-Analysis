#!/usr/bin/env python3
"""
Epitope Mapping & Visualization Tool

Combines epitope-to-sequence mapping with automated heatmap generation.
Maps epitopes to sequences, calculates similarity scores, and creates visualizations.
"""

import csv
import re
import sys
from pathlib import Path

import pandas as pd
import numpy as np
from Bio import SeqIO

try:
    import seaborn as sns
    import matplotlib.pyplot as plt
    VISUALIZATION_AVAILABLE = True
except ImportError:
    VISUALIZATION_AVAILABLE = False
    print("⚠️  Warning: seaborn/matplotlib not available. Visualization will be skipped.")


# ============================================
# Core Functions - Epitope Mapping
# ============================================

def hamming_distance(s1, s2):
    """Calculate the Hamming distance between two strings."""
    if len(s1) != len(s2):
        raise ValueError("The strings must have the same length.")
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def extract_subgenus(header: str):
    """
    Extract genus and species from FASTA header.

    Example:
    [species=Altai virus]
    [genus=Orthohantavirus]

    Returns:
        "Orthohantavirus_Altai virus"
    """
    genus_match = re.search(r"\[genus=([^\]]+)\]", header)
    species_match = re.search(r"\[species=([^\]]+)\]", header)

    genus = genus_match.group(1).strip() if genus_match else ""
    species = species_match.group(1).strip() if species_match else ""

    if genus and species:
        return f"{genus}_{species}"
    elif genus:
        return genus
    elif species:
        return species
    else:
        return ""


def process_sequences(
    epitope_file,
    fasta_file,
    output_file,
    flank_size=None,
    epitope_seq_col="Sequence",
    start_col="Start",
    end_col="End",
):
    """
    Map epitopes to sequences and save results to Excel.

    flank_size=None → full search
    flank_size=int → windowed search ± flank_size aa
    """
    epitope_file = Path(epitope_file)
    fasta_file = Path(fasta_file)
    output_file = Path(output_file)

    # -------------------------
    # Load epitope CSV (robust header resolution)
    # -------------------------
    epitopes = []
    with epitope_file.open("r", newline="") as f:
        reader = csv.DictReader(f)
        raw_fields = reader.fieldnames or []
        norm_map = {name.strip().lower(): name for name in raw_fields}

        def resolve(col):
            key = col.strip().lower()
            if key in norm_map:
                return norm_map[key]
            raise KeyError(f"Column '{col}' not found. Available: {raw_fields}")

        seq_col_name = resolve(epitope_seq_col)
        start_col_name = resolve(start_col)
        end_col_name = resolve(end_col)

        for row in reader:
            epitopes.append({
                "seq": row[seq_col_name],
                "start": int(row[start_col_name]),
                "end": int(row[end_col_name]),
            })

    # -------------------------
    # Load FASTA
    # -------------------------
    sequences, headers, subgenera = [], [], []
    for record in SeqIO.parse(str(fasta_file), "fasta"):
        headers.append(record.id)
        sequences.append(str(record.seq))
        subgenera.append(extract_subgenus(record.description))

    # -------------------------
    # Output columns
    # -------------------------
    columns = ["epitope"]
    for h, sg in zip(headers, subgenera):
        columns += [
            f"{h}_{sg}_position",
            f"{h}_{sg}_score",
            f"{h}_{sg}_best_match",
        ]

    results_data = []
    multiple_best = []

    # -------------------------
    # Compute matching
    # -------------------------
    for epi in epitopes:
        epitope = epi["seq"]
        start_0 = epi["start"] - 1
        end_0 = epi["end"] - 1
        L = len(epitope)

        row = [epitope]

        for seq, h, sg in zip(sequences, headers, subgenera):
            seqlen = len(seq)

            if seqlen < L:
                row += [None, float("-inf"), None]
                continue

            if flank_size is None:
                # Full search
                w_start, w_end = 0, seqlen
            else:
                # Windowed
                w_start = max(0, start_0 - flank_size)
                w_end = min(seqlen, end_0 + flank_size + 1)

            if w_end - w_start < L:
                row += [None, float("-inf"), None]
                continue

            max_score = float("-inf")
            best_pos = None
            best_sub = None
            ties_pos = []
            ties_sub = []

            # sliding window
            for i in range(w_start, w_end - L + 1):
                subseq = seq[i:i+L]
                score = (L - hamming_distance(epitope, subseq)) / L

                if score > max_score:
                    max_score = score
                    best_pos = i
                    best_sub = subseq
                    ties_pos = [i]
                    ties_sub = [subseq]
                elif score == max_score:
                    ties_pos.append(i)
                    ties_sub.append(subseq)

            row += [best_pos, max_score, best_sub]

            if len(ties_pos) > 1:
                multiple_best.append({
                    "epitope": epitope,
                    "header": h,
                    "subgenus": sg,
                    "max_score": max_score,
                    "all_positions": ";".join(map(str, ties_pos)),
                    "all_matches": ";".join(ties_sub),
                    "chosen_position": best_pos,
                    "chosen_match": best_sub,
                    "flank_size": flank_size if flank_size is not None else "full",
                })

        results_data.append(row)

    # -------------------------
    # Save Excel
    # -------------------------
    df = pd.DataFrame(results_data, columns=columns)

    # ✅ Convert main *_position columns to 1-based
    for col in df.columns:
        if col.endswith("_position"):
            df[col] = df[col].apply(lambda x: x + 1 if pd.notna(x) else x)

    # ✅ Convert Multiple_Best_Matches positions to 1-based
    for rec in multiple_best:
        if rec["chosen_position"] is not None:
            rec["chosen_position"] += 1

        if rec["all_positions"]:
            rec["all_positions"] = ";".join(
                str(int(p) + 1) for p in rec["all_positions"].split(";")
            )

    df_scores = df[["epitope"] + [c for c in df.columns if c.endswith("_score")]]

    if multiple_best:
        df_ties = pd.DataFrame(multiple_best)
    else:
        df_ties = pd.DataFrame(columns=[
            "epitope","header","subgenus","max_score",
            "all_positions","all_matches",
            "chosen_position","chosen_match","flank_size"
        ])

    with pd.ExcelWriter(output_file) as w:
        df.to_excel(w, sheet_name="Original", index=False)
        df_scores.to_excel(w, sheet_name="All_Scores", index=False)
        df_ties.to_excel(w, sheet_name="Multiple_Best_Matches", index=False)

    print(f"✅ Saved: {output_file}")
    return output_file


# ============================================
# Visualization Functions
# ============================================

def protein_name_from_filename(filename: str) -> str:
    """Determine protein name from filename."""
    name = filename.lower()

    if name.startswith("rdrp_"):
        return "RdRp Protein"
    elif name.startswith("nonstructural_"):
        return "NS Protein"
    elif name.startswith("nucleocapsid_"):
        return "N Protein"
    elif name.startswith("gngc_"):
        return "GnGc Protein"
    else:
        return "Protein"


def extract_species_genus(col):
    """
    Extract species and genus from column name.
    Format: QOI08666.1_Orthohantavirus_Altai virus_score
    """
    name = str(col).removesuffix("_score")
    parts = name.split("_")

    if len(parts) >= 3:
        genus = parts[-2]
        species = parts[-1]
        return species, genus

    return "Unknown", "Unknown"


def generate_heatmap_visualization(excel_file, protein_name=None, colormap="RdBu"):
    """
    Generate heatmap visualization from epitope similarity scores.
    
    Args:
        excel_file: Path to Excel file with All_Scores sheet
        protein_name: Custom protein name for title
        colormap: Colormap to use (default: RdBu)
    
    Returns:
        Path to SVG output file or None if failed
    """
    if not VISUALIZATION_AVAILABLE:
        print("⚠️  Skipping visualization: seaborn/matplotlib not available")
        return None

    excel_file = Path(excel_file)

    try:
        results_df = pd.read_excel(excel_file, sheet_name="All_Scores")
    except Exception as e:
        print(f"❌ Error reading {excel_file}: {e}")
        return None

    score_columns = [c for c in results_df.columns if str(c).endswith("_score")]
    if not score_columns:
        print(f"⚠️  No score columns found in {excel_file.name}")
        return None

    heatmap_data = results_df[score_columns].copy()

    mask_empty = results_df["epitope"].isna() | (
        results_df["epitope"].astype(str).str.strip() == ""
    )
    heatmap_data.loc[mask_empty, :] = np.nan

    species_genus = [extract_species_genus(c) for c in score_columns]
    species = [x[0] for x in species_genus]
    genera = [x[1] for x in species_genus]

    score_columns_sorted = score_columns
    species_sorted = species
    genera_sorted = genera

    heatmap_data = heatmap_data[score_columns_sorted]

    heatmap_matrix = heatmap_data.to_numpy(dtype=float)
    mask = np.isnan(heatmap_matrix)
    nrows, ncols = heatmap_matrix.shape

    # Species separator boundaries
    species_change_points = [
        i for i in range(1, ncols)
        if species_sorted[i] != species_sorted[i - 1]
    ]
    species_boundaries = [0] + species_change_points + [ncols]

    # Genus label centers
    genus_change_points = [
        i for i in range(1, ncols)
        if genera_sorted[i] != genera_sorted[i - 1]
    ]
    genus_boundaries = [0] + genus_change_points + [ncols]

    genus_names = []
    genus_centers = []

    for start, end in zip(genus_boundaries[:-1], genus_boundaries[1:]):
        genus_names.append(genera_sorted[start])
        genus_centers.append((start + end) / 2.0)

    epitope_numbers = [str(i) for i in range(1, nrows + 1)]

    ystep = 1 if nrows <= 50 else (2 if nrows <= 100 else 5)

    plt.figure(figsize=(18, 18))

    ax = sns.heatmap(
        heatmap_matrix,
        mask=mask,
        cmap=colormap,
        vmin=0,
        vmax=1,
        cbar_kws={"label": "Similarity Score"},
        linewidths=0
    )

    sns.heatmap(
        np.zeros_like(heatmap_matrix),
        mask=~mask,
        cmap=["#4d4d4d"],
        cbar=False,
        ax=ax
    )

    # Species separator lines
    for b in species_boundaries[1:-1]:
        ax.axvline(b, color="black", linewidth=0.35)

    # Hide individual sequence x labels
    ax.set_xticks([])
    ax.set_xticklabels([])

    # Genus labels on x axis
    secax = ax.secondary_xaxis("bottom")
    secax.set_xticks(genus_centers)
    secax.set_xticklabels(
        genus_names,
        rotation=90,
        ha="center",
        va="top",
        fontsize=6,
        fontweight="bold"
    )
    secax.spines["bottom"].set_visible(False)
    secax.tick_params(axis="x", length=0, pad=20)
    secax.set_xlabel("Genus", labelpad=35)

    # Y-axis peptide labels
    yticks = np.arange(0, nrows, ystep) + 0.5
    ax.set_yticks(yticks)
    ax.set_yticklabels(
        [epitope_numbers[i] for i in range(0, nrows, ystep)],
        fontsize=9
    )

    ax.set_ylabel("Peptide Index", fontsize=12)

    title_name = excel_file.stem.replace("flank", "window")
    if protein_name is None:
        protein_name = protein_name_from_filename(excel_file.name)

    ax.set_title(
        f"{protein_name} Peptide Similarity Grouped by Genus\n{title_name}",
        fontsize=14,
        pad=20
    )

    plt.subplots_adjust(bottom=0.42)

    output_svg = excel_file.with_suffix(".svg")
    plt.savefig(output_svg, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"📊 Saved visualization: {output_svg}")
    return output_svg


# ============================================
# User Input Functions
# ============================================

def get_user_inputs():
    """Prompt user for all required inputs."""
    print("\n" + "="*60)
    print("EPITOPE MAPPING & VISUALIZATION TOOL")
    print("="*60 + "\n")

    # Ask if multiple proteins
    multi_protein = input("Process multiple proteins? (y/n, default: n): ").strip().lower()
    multi_protein = multi_protein == "y"

    proteins_config = {}

    if multi_protein:
        # Multiple proteins mode
        num_proteins = input("How many proteins? (e.g., 3): ").strip()
        try:
            num_proteins = int(num_proteins)
        except ValueError:
            print("Invalid input. Using 1 protein.")
            num_proteins = 1

        for i in range(num_proteins):
            print(f"\n--- Protein {i+1} of {num_proteins} ---")
            
            protein_name = input(f"🧪 Enter protein name #{i+1}: ").strip()
            if not protein_name:
                protein_name = f"protein_{i+1}"

            epitope_file = input(f"📄 Epitope CSV file for {protein_name}: ").strip()
            while not Path(epitope_file).exists():
                print(f"❌ File not found: {epitope_file}")
                epitope_file = input(f"📄 Epitope CSV file for {protein_name}: ").strip()

            fasta_file = input(f"🧬 FASTA file for {protein_name}: ").strip()
            while not Path(fasta_file).exists():
                print(f"❌ File not found: {fasta_file}")
                fasta_file = input(f"🧬 FASTA file for {protein_name}: ").strip()

            flank_input = input(f"📏 Flank sizes for {protein_name} (e.g., '25,50,100' or leave empty): ").strip()
            if flank_input:
                try:
                    flanks = [int(x.strip()) for x in flank_input.split(",")]
                except ValueError:
                    print("⚠️  Invalid flank input. Using full search only.")
                    flanks = [None]
            else:
                flanks = [None]

            proteins_config[protein_name] = {
                "epitope_file": epitope_file,
                "fasta_file": fasta_file,
                "flanks": flanks,
            }

    else:
        # Single protein mode
        protein_name = input("🧪 Enter protein name: ").strip()
        if not protein_name:
            protein_name = "protein"

        epitope_file = input("📄 Enter epitope CSV file path: ").strip()
        while not Path(epitope_file).exists():
            print(f"❌ File not found: {epitope_file}")
            epitope_file = input("📄 Enter epitope CSV file path: ").strip()

        fasta_file = input("🧬 Enter FASTA sequence file path: ").strip()
        while not Path(fasta_file).exists():
            print(f"❌ File not found: {fasta_file}")
            fasta_file = input("🧬 Enter FASTA sequence file path: ").strip()

        flank_input = input("📏 Enter flank sizes (comma-separated, e.g. '25,50,100' or leave empty): ").strip()
        if flank_input:
            try:
                flanks = [int(x.strip()) for x in flank_input.split(",")]
            except ValueError:
                print("⚠️  Invalid flank input. Using full search only.")
                flanks = [None]
        else:
            flanks = [None]

        proteins_config[protein_name] = {
            "epitope_file": epitope_file,
            "fasta_file": fasta_file,
            "flanks": flanks,
        }

    # Output directory (shared)
    output_dir = input("\n📁 Enter output directory (default: results): ").strip()
    output_dir = output_dir or "results"
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    # Column names (shared)
    epitope_col = input("📋 Epitope sequence column name (default: Sequence): ").strip()
    epitope_col = epitope_col or "Sequence"

    start_col = input("📋 Start position column name (default: Start): ").strip()
    start_col = start_col or "Start"

    end_col = input("📋 End position column name (default: End): ").strip()
    end_col = end_col or "End"

    # Visualization
    print("\n📊 Visualization Options:")
    generate_viz = input("Generate heatmap visualization? (y/n, default: n): ").strip().lower()
    generate_viz = generate_viz == "y"

    colormap = "RdBu"
    if generate_viz:
        colormap_input = input("   Colormap (default: RdBu): ").strip()
        colormap = colormap_input or "RdBu"

    return {
        "proteins": proteins_config,
        "output_dir": output_dir,
        "epitope_col": epitope_col,
        "start_col": start_col,
        "end_col": end_col,
        "generate_viz": generate_viz,
        "colormap": colormap,
    }


# ============================================
# Main Execution
# ============================================

def main():
    """Main execution function."""
    try:
        # Get user inputs
        config = get_user_inputs()

        print("\n" + "="*60)
        print("RUNNING EPITOPE MAPPING")
        print("="*60 + "\n")

        proteins = config["proteins"]
        output_dir = config["output_dir"]
        epitope_col = config["epitope_col"]
        start_col = config["start_col"]
        end_col = config["end_col"]
        generate_viz = config["generate_viz"]
        colormap = config["colormap"]

        # Process each protein
        total_proteins = len(proteins)
        for idx, (protein_name, protein_data) in enumerate(proteins.items(), 1):
            print(f"\n📊 Processing protein {idx}/{total_proteins}: {protein_name}")
            print("-" * 60)

            epitope_file = protein_data["epitope_file"]
            fasta_file = protein_data["fasta_file"]
            flanks = protein_data["flanks"]

            # Run full search
            print(f"🔄 Processing full search...")
            out_full = output_dir / f"{protein_name}_fullsearch.xlsx"
            excel_full = process_sequences(
                epitope_file=epitope_file,
                fasta_file=fasta_file,
                output_file=out_full,
                flank_size=None,
                epitope_seq_col=epitope_col,
                start_col=start_col,
                end_col=end_col,
            )

            # Generate visualization for full search
            if generate_viz and VISUALIZATION_AVAILABLE:
                generate_heatmap_visualization(excel_full, protein_name=protein_name, colormap=colormap)

            # Run flank searches
            for flank in flanks:
                if flank is not None:
                    print(f"🔄 Processing flank={flank}...")
                    out_flank = output_dir / f"{protein_name}_flank{flank}.xlsx"
                    excel_flank = process_sequences(
                        epitope_file=epitope_file,
                        fasta_file=fasta_file,
                        output_file=out_flank,
                        flank_size=flank,
                        epitope_seq_col=epitope_col,
                        start_col=start_col,
                        end_col=end_col,
                    )

                    # Generate visualization for flank search
                    if generate_viz and VISUALIZATION_AVAILABLE:
                        generate_heatmap_visualization(excel_flank, protein_name=protein_name, colormap=colormap)

        print("\n" + "="*60)
        print("✅ ANALYSIS COMPLETE!")
        print("="*60)
        print(f"📁 Results saved to: {output_dir}")
        print(f"🧪 Processed {total_proteins} protein(s)")

    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
