# Changelog

## Version 2.1.0 (2026-10-07)

### Pairing
- Direction word can sit anywhere in the name, e.g. `Mitf_Fwd_IM` (previously these all came out unpaired)
- Unpaired Fwd/Rev primers on neighbouring rows are paired (order sheets list pairs on alternate lines)
- Near-identical names are paired as a probable typo (`Tyrp1_Fwd_IM` + `Tryp1_Rev_IM`) and flagged WARN
- Gene names not found in NCBI are retried with neighbouring letters swapped and NCBI's spelling suggestion (`Endrb` -> Ednrb); only official symbols are accepted, and the pair is flagged WARN

## Version 2.0.0 (2026-10-07)

### Bug fixes
- Gene names starting with "Chr" (Chrm1, Chrna7, Chrd...) were treated as chromosomes; only exact tokens (Chr1, ChrX, ChrMT) count now
- Fwd/Rev detection matched anywhere in the name, so `Rev3l_Fwd` counted as both; it now requires a direction suffix (case-insensitive, with optional pair number: `_F1`/`_R1`)
- Chromosome loci are searched on both strands; the whole-chromosome efetch (which timed out) is replaced by a cached UCSC download, and `chromosome 1[Title]` no longer matches chr10-19
- Genomic sequence now comes from the gene's own coordinates in NCBI Gene (not the first genomic search hit, which could be the wrong record)
- Transcripts come from NCBI Gene links, so only that gene's RefSeq RNAs are checked, and all of them (not 5 of 10)
- NCBI requests are rate-limited and retried; network failures are reported as ERROR rather than "not found"
- Input header row is detected instead of assumed; sheets with fewer than 10 columns no longer crash; CSV supported

### New
- Mismatch tolerance (`--max-mismatches`, default 1, 3' end must match) and degenerate IUPAC bases
- Reports every product: multiple products, single-primer products and swapped Fwd/Rev are flagged
- Genomic DNA product size per pair (shows whether the primers span an intron)
- Which/how many transcripts each pair amplifies
- Products inferred between unpartnered locus primers (e.g. Chr1_L_Flank + Chr1_L_Int)
- `--genome-check`: genome-wide off-target product search
- Nearest-neighbour Tm (Primer3) replaces the Wallace rule
- Output has Pairs, Primers and Settings sheets; size window, product limit, flank etc. are options
- NCBI results and genome files cached in ~/.cache/primer_validation_tool
- Offline unit tests (`pytest tests`)
- `setup.sh` creates the environment outside the repo (~/.venvs); `./check_primers` wrapper runs setup automatically when needed

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
