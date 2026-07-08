# Epitope Analyzer - UPDATED README

## ✨ NEW: Multiple Protein Support!

Process **2, 3, 5, or more proteins** in a single run!

---

## What It Does

1. **Maps epitopes to sequences** - Finds where each epitope matches in your sequences
2. **Calculates similarity scores** - Uses Hamming distance (0-1 scale)
3. **Creates Excel output** - All results in organized sheets
4. **Generates heatmaps** - Beautiful visualizations grouped by genus/species
5. **🆕 Handles multiple proteins** - Process many proteins at once!

---

## Installation

### Requirements
```bash
pip install pandas biopython seaborn matplotlib numpy
```

### Make Executable (Linux/Mac)
```bash
chmod +x epitope_analyzer.py
```

---

## How to Run

```bash
python epitope_analyzer.py
```

Or:
```bash
./epitope_analyzer.py
```

---

## TWO MODES

### Mode 1: Single Protein

```
Process multiple proteins? (y/n, default: n): n
🧪 Enter protein name: RdRp
📄 Enter epitope CSV file path: rdrp_epitopes.csv
🧬 Enter FASTA sequence file path: rdrp_sequences.fasta
📏 Enter flank sizes (e.g., '25,50,100'): 50,100
...
```

### Mode 2: Multiple Proteins (NEW!)

```
Process multiple proteins? (y/n, default: n): y
How many proteins? (e.g., 3): 3

--- Protein 1 of 3 ---
🧪 Enter protein name #1: RdRp
📄 Epitope CSV file for RdRp: rdrp_epitopes.csv
🧬 FASTA file for RdRp: rdrp_sequences.fasta
📏 Flank sizes for RdRp: 50,100

--- Protein 2 of 3 ---
🧪 Enter protein name #2: GnGc
📄 Epitope CSV file for GnGc: gngc_epitopes.csv
🧬 FASTA file for GnGc: gngc_sequences.fasta
📏 Flank sizes for GnGc: 50,100

--- Protein 3 of 3 ---
🧪 Enter protein name #3: Nucleocapsid
📄 Epitope CSV file for Nucleocapsid: n_epitopes.csv
🧬 FASTA file for Nucleocapsid: n_sequences.fasta
📏 Flank sizes for Nucleocapsid: 50

📁 Enter output directory (default: results): analysis
📊 Generate heatmap visualization? (y/n): y
   Colormap (default: RdBu): RdBu
```

**All proteins processed in ONE RUN!** ✅

---

## Example Output

Running 3 proteins produces:

```
analysis/
├── RdRp_fullsearch.xlsx          <- Full search
├── RdRp_fullsearch.svg           <- Heatmap
├── RdRp_flank50.xlsx
├── RdRp_flank50.svg
├── RdRp_flank100.xlsx
├── RdRp_flank100.svg
│
├── GnGc_fullsearch.xlsx
├── GnGc_fullsearch.svg
├── GnGc_flank50.xlsx
├── GnGc_flank50.svg
├── GnGc_flank100.xlsx
├── GnGc_flank100.svg
│
├── Nucleocapsid_fullsearch.xlsx
├── Nucleocapsid_fullsearch.svg
├── Nucleocapsid_flank50.xlsx
└── Nucleocapsid_flank50.svg
```

**18 files created** in one execution! 🎉

---

## Input File Formats

### Epitope CSV File

Required columns: `Sequence`, `Start`, `End` (case-insensitive)

Example:
```csv
Sequence,Start,End
MFFLLLLAAA,1,10
LLLAAAMMMM,15,24
```

### FASTA Sequence File

Standard FASTA format:

```fasta
>seq_001 [genus=Orthohantavirus] [species=Altai virus]
MFFLLLLAAA...

>seq_002 [genus=Orthohantavirus] [species=Seoul virus]
LLLAAAMMMM...
```

---

## Output Files

### Excel Files (per protein)
- `protein_fullsearch.xlsx` - Full search
- `protein_flank50.xlsx` - ±50 aa window
- `protein_flank100.xlsx` - ±100 aa window

Each Excel file has 3 sheets:
1. **Original** - Position, score, match
2. **All_Scores** - Score columns only
3. **Multiple_Best_Matches** - Ties and ambiguities

### SVG Visualizations (if enabled)
- `protein_fullsearch.svg` - Heatmap
- `protein_flank50.svg` - Heatmap
- etc.

---

## Key Features

✅ **Single or multiple proteins** - Your choice  
✅ **Different flanks per protein** - Customize each one  
✅ **Shared output folder** - All results in one place  
✅ **All visualizations** - Automatic heatmaps for all proteins  
✅ **Progress tracking** - Shows "X of Y proteins"  
✅ **Same column names** - Use one naming scheme for all  
✅ **No code editing** - Fully interactive  
✅ **Error handling** - Clear messages if something fails  

---

## Tips for Multiple Proteins

1. **Same column names** - All CSVs should have same headers (Sequence, Start, End)
2. **Different flanks** - Can specify different flank sizes for each protein
3. **One output folder** - All proteins' results go to same directory
4. **Progress display** - Script shows "Processing protein 2/5" etc.
5. **Consistent format** - FASTA headers should follow same pattern

---

## Example: Processing 4 Proteins

```bash
python epitope_analyzer.py

Process multiple proteins? (y/n, default: n): y
How many proteins? (e.g., 3): 4

[Input for RdRp, GnGc, N_protein, NSL...]

📁 Enter output directory: my_analysis

============================================================
RUNNING EPITOPE MAPPING
============================================================

📊 Processing protein 1/4: RdRp
------------------------------------------------------------
🔄 Processing full search...
✅ Saved: my_analysis/RdRp_fullsearch.xlsx
📊 Saved visualization: my_analysis/RdRp_fullsearch.svg
[... flank results ...]

📊 Processing protein 2/4: GnGc
------------------------------------------------------------
🔄 Processing full search...
✅ Saved: my_analysis/GnGc_fullsearch.xlsx
📊 Saved visualization: my_analysis/GnGc_fullsearch.svg
[... flank results ...]

📊 Processing protein 3/4: N_protein
------------------------------------------------------------
[... results ...]

📊 Processing protein 4/4: NSL
------------------------------------------------------------
[... results ...]

============================================================
✅ ANALYSIS COMPLETE!
============================================================
📁 Results saved to: my_analysis
🧪 Processed 4 protein(s)
```

---

## Troubleshooting

### "File not found" Error
- Check file paths are correct
- Use absolute paths: `/home/user/data/sequences.fasta`

### "Column not found" Error
- Verify CSV headers match across all proteins
- Column names are case-insensitive
- Check spelling: `Sequence`, `Start`, `End`

### Script stops on one protein
- Fix the file or CSV format for that protein
- Error message will tell you which file failed
- Re-run after fixing

### No visualizations
- Install: `pip install seaborn matplotlib`
- Check Excel file has "All_Scores" sheet

---

## Dependencies

```bash
pip install pandas biopython seaborn matplotlib numpy
```

---

## Single vs Multiple Proteins

| Feature | Single | Multiple |
|---------|--------|----------|
| # Proteins | 1 | 2+ |
| Setup time | Quick | 2-3 min |
| Process time | Fast | Takes longer |
| Output files | 2-4 | 6-20+ |
| Flank sizes | Same for all | Different per protein |
| Column names | One set | Must be same |

---

## Output Interpretation

**Excel sheets per protein:**

**Original:**
- Epitope names
- Position (1-based)
- Similarity scores (0-1)
- Matching sequences

**All_Scores:**
- Epitope + scores only
- Easy filtering/sorting

**Multiple_Best_Matches:**
- Epitopes with tied matches
- All positions found

**Heatmaps:**
- Y: Epitope indices
- X: Grouped by genus
- Color: Score (Red=high, Blue=low)
- Gray: Missing data

---

## Quick Workflow

### Process 3 Proteins:

```bash
# 1. Run script
python epitope_analyzer.py

# 2. Answer: y to multiple proteins
# 3. Enter: 3
# 4. For each protein:
#    - Name
#    - Epitope CSV path
#    - FASTA path
#    - Flank sizes
# 5. Common settings (column names, viz, colormap)
# 6. Wait for all proteins to process
# 7. Check output folder
```

**That's it!** All 3 proteins analyzed. 🎯

---

Good luck! Questions? Check the file formats above.
