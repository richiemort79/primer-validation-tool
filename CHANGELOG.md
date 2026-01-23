# Changelog

## Version 1.1.0 (2025-01-22)

### New Features
- **Multiple isoform checking**: Script now fetches and checks up to 5 mRNA isoforms per gene
- **Transcript ID reporting**: New `Transcript_ID` column shows which specific RefSeq transcript matched (e.g., NM_008601.4)
- This solves issues where primers are designed for a specific isoform/variant that differs from the primary RefSeq

### Why This Matters
- Genes often have multiple transcript variants
- Primers designed from CCDS or specific isoforms may not match the "primary" RefSeq mRNA
- Now the script tries all available isoforms and reports which one worked

## Version 1.0.1 (2025-01-22)

### Bug Fixes
- **CRITICAL: Tm calculation now correct for all primers**
  - Previously used incorrect formula for primers ≥14bp
  - Now correctly uses Wallace rule for ALL primers: Tm = 2(A+T) + 4(G+C)
  - All primers now give accurate Tm values

## Version 1.0.0 (2025-01-22)

### Changes
- **Amplicon size range updated**: Changed from 100-150 bp to 70-150 bp
- **Tm calculation fixed**: Now uses Wallace rule correctly (Tm = 2(A+T) + 4(G+C))
- **Column naming simplified**: Single Tm column (Fwd_Tm_C / Rev_Tm_C) instead of two redundant columns

### Bug Fixes
- Fixed incorrect Tm calculation that was giving values lower than expected
- Corrected the "nearest neighbor" formula which was actually wrong

### Technical Details

**Tm Calculation:**
- Old (incorrect): `tm = 64.9 + 41 * (gc_count - 16.4) / length` for primers >14bp
- New (correct): Wallace rule for all primers: `Tm = 2(A+T) + 4(G+C)`
- This gives accurate Tm values that match the simple G/C=4°C, A/T=2°C formula

**Amplicon Range:**
- Changed validation range from 100-150 bp to 70-150 bp per user requirements
- All status messages and documentation updated accordingly

### Notes
- The Wallace rule (2(A+T) + 4(G+C)) and the simple formula (GC*4 + AT*2) are mathematically identical
- For more accurate Tm calculations considering salt concentration, use specialized tools like IDT OligoAnalyzer
