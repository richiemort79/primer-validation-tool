# Primer Validation Tool

A Python script that validates PCR primer pairs by checking if they bind to the correct mouse genes and generate amplicons in the expected size range (70-150 bp). Automatically checks both genomic DNA and mRNA sequences to handle primers that span exon boundaries.

## Features

- ✅ Validates primers against NCBI mouse gene sequences
- ✅ Checks both mRNA (spliced) and genomic DNA (with introns)
- ✅ Calculates amplicon sizes
- ✅ Validates amplicon size is 70-150 bp
- ✅ Calculates primer statistics (length, GC%, Tm using Wallace rule)
- ✅ Generates detailed Excel report
- ✅ Handles primers designed for RT-PCR that span exons
- ✅ Batch processing of multiple primer pairs

## Requirements

- Python 3.6 or higher
- Internet connection (to access NCBI databases)

## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/primer-validation-tool.git
cd primer-validation-tool
```

2. Create a virtual environment (recommended):
```bash
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### Basic Usage

```bash
python check_primers.py your_primers.xlsx
```

### Interactive Mode

```bash
python check_primers.py
```
The script will prompt you for the input file path.

### Custom Output Name

```bash
python check_primers.py input.xlsx -o custom_results.xlsx
```

### Different Organism

```bash
python check_primers.py primers.xlsx --organism "Rattus norvegicus"
```

## Input File Format

The script expects an Excel file in IDT template format:

- Primers starting from row 5 (row 4 is the header)
- Column A: Primer names (e.g., `Gene_Fwd`, `Gene_Rev`)
- Column C: Primer sequences
- Primer names must contain `_Fwd` or `_Rev` (case-insensitive)

Example:
```
| Oligo Name  | ... | Sequence              | ...
|-------------|-----|-----------------------|-----
| Mitf_Fwd    | ... | CCAGGTGCCGATGGAAGTC  | ...
| Mitf_Rev    | ... | GGTGAGCTCAGGACTTGGC  | ...
| Kit_Fwd     | ... | CTGGCTCTGGACCTGGATG  | ...
| Kit_Rev     | ... | GATGTCTCTGGCTAGCCCG  | ...
```

## Output

The script generates an Excel file with the following information:

| Column | Description |
|--------|-------------|
| **Gene** | Gene name |
| **Status** | PASS / WARN / FAIL / ERROR |
| **Amplicon_Size_bp** | Predicted amplicon size |
| **Sequence_Type** | mRNA / Genomic / mRNA (Genomic also checked) |
| **Details** | Explanation of result |
| **Fwd_Primer / Rev_Primer** | Primer sequences |
| **Fwd_Length / Rev_Length** | Primer lengths |
| **Fwd_GC% / Rev_GC%** | GC content |
| **Fwd_Tm_C / Rev_Tm_C** | Melting temp in °C (Wallace rule: 2(A+T) + 4(G+C)) |

### Status Meanings

- **PASS**: Primers bind correctly and generate 70-150 bp amplicon ✓
- **WARN**: Primers bind but amplicon size is outside target range ⚠️
- **FAIL**: Primers don't bind or generate no amplicon ✗
- **ERROR**: Missing primers or gene not found ⚠️

### Sequence Types

- **mRNA**: Primers bind to spliced mRNA (no introns)
- **Genomic**: Primers bind to genomic DNA (with introns)
- **Genomic (mRNA also checked)**: Both sequences tested

## How It Works

1. Extracts primer pairs from Excel file
2. Fetches both mRNA and genomic DNA sequences from NCBI for each gene
3. Checks mRNA first (most primers are designed for RT-PCR)
4. Falls back to genomic DNA if mRNA doesn't work
5. Calculates amplicon sizes and validates 70-150 bp range
6. Generates comprehensive Excel report

## Example Output

```
======================================================================
Validating 33 primer pairs against Mus musculus genes
======================================================================

[1/33] Checking Mitf...
  Fetching sequences for Mitf... ✓ (mRNA: 1234 bp, Genomic: 5678 bp)
  ✓ Amplicon size: 125 bp ✓ [mRNA]

[2/33] Checking Kit...
  Fetching sequences for Kit... ✓ (mRNA: 2345 bp, Genomic: 6789 bp)
  ✓ Amplicon size: 142 bp ✓ [mRNA]

======================================================================
VALIDATION SUMMARY
======================================================================

Total primer pairs checked: 33
✓ PASS (70-150 bp):         28 (84.8%)
⚠ WARN (wrong size):        3 (9.1%)
✗ FAIL (no amplicon):       2 (6.1%)
⚠ ERROR (missing primers):  0 (0.0%)
```

## Troubleshooting

### "Gene sequence not found in NCBI"
- Verify gene name is correct (use official gene symbols)
- Check gene exists in NCBI database for your organism

### "Primer does not bind to gene"
- Primer may span exon boundary (script checks both mRNA and genomic)
- Verify primer sequence is correct
- Primer may not be specific to this gene

### No internet connection
- Script requires internet to access NCBI databases
- Ensure Python can access `eutils.ncbi.nlm.nih.gov`

## Dependencies

- pandas >= 1.3.0
- openpyxl >= 3.0.0
- requests >= 2.25.0

## License

MIT License - see LICENSE file for details

## Contributing

Contributions welcome! Please feel free to submit a Pull Request.

## Citation

If you use this tool in your research, please cite:
```
Primer Validation Tool
https://github.com/yourusername/primer-validation-tool
```

## Support

For issues or questions, please open an issue on GitHub.
