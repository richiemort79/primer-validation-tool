# Changelog

## Version 2.4.2 (2026-10-15)

### Isoforms
- Removed the `Isoform` column: an isoform is now named only in the primer name (`Mitf-M_Fwd`), so upload sheets need no changes
- The catch-all warning now suggests naming the primers `Gene-<isoform>_Fwd/_Rev` instead of adding an Isoform column

## Version 2.4.1 (2026-10-15)

### Amplicon size
- Default size window widened from 70-150 to 70-175 bp (`--max-size`): products up to ~200 bp amplify efficiently with standard SYBR cycling
- New Primer3 designs still aim for ~100 bp: each bp over that adds to the design penalty

## Version 2.4.0 (2026-10-15)

### Isoforms
- Gene pairs are checked as catch-all by default: a pair that misses any curated (NM_/NR_) transcript now WARNs and names it. Missed predicted (XM_/XR_) transcripts are only noted
- New optional `Isoform` (or `Transcript`) column: list RefSeq accessions to make a pair isoform-specific (e.g. Mitf-M: NM_008601, NM_001410059). It WARNs if it misses a target or amplifies another curated transcript
- A hyphen after the gene in a primer name (`Mitf-M_Fwd`, `MITF-M_For`) marks the pair as isoform-targeted, so Sigma/IDT upload sheets work unchanged: no catch-all WARN, Notes lists the curated transcripts it amplifies, and `--suggest` aims for the same ones
- `Transcripts_Amplified` and `Matched_Transcripts` count curated transcripts only; predicted XM_/XR_ transcripts are ignored (e.g. Glis3 shows 3/3, not 4/7)
- Gene lookup also tries the part of the name before a `-` (`MITF-M` -> MITF)
- New `Isoforms` column on the Pairs and Alternatives sheets

### Suggestions
- Catch-all pairs get extra Primer3 designs with both primers in sequence shared by every curated transcript; isoform-specific pairs get designs with a primer in sequence unique to the targets
- Ranking puts isoform coverage first, then a single product size across isoforms
- Isoform-specific suggestions never amplify other curated transcripts

## Version 2.3.0 (2026-10-13)

### Tm
- Tm is calculated at SYBR mastermix conditions (3 mM Mg2+, 250 nM primer; `--mg`, `--primer-conc`) instead of Primer3's defaults (1.5 mM Mg2+, 50 nM), so values are ~4 °C higher than before and closer to the real reaction
- New `Tm_2+4` columns show the 2+4 (Wallace) rule for comparison
- Pairs WARN if a primer's Tm is below the annealing temperature (`--anneal`, default 60 °C), or the two Tms differ by more than `--max-tm-diff` (default 5 °C)

### G/C clamps
- New `GC_Clamp` columns show which primer ends are G/C, and flag GC-rich 3' ends (>3 G/C in the last 5); a missing 3' clamp is mentioned in Notes
- Primer3 designs require a 3' G/C clamp (with at most 3 G/C in the last 5) and a 5' G/C (`--no-5prime-clamp` to drop the 5' rule); clamped pairs rank higher

### Suggestions
- For primers with a low or mismatched Tm, `--suggest` offers your own pair lengthened at the 5' end (3' end unchanged), ranked first when it passes cleanly
- Primer3 designs aim 2-7 °C above the annealing temperature, at the same reaction conditions
- `--suggest` only looks for replacements for WARN/FAIL pairs; PASS pairs that don't span an intron are no longer included
- With `--suggest`, the Pairs sheet ends with a `Suggested_Alternative` column holding the top-ranked replacement pair

## Version 2.2.0 (2026-10-07)

### New
- `--suggest`: for pairs that don't span an intron or that WARN/FAIL, suggests up to `--suggest-count` (default 3) better pairs on a new Alternatives sheet
  - candidates from Harvard PrimerBank (mouse/human; lab-validated pairs marked) and Primer3 designs around every exon junction (junction-crossing, or flanking introns too long to amplify)
  - exon structure derived from the transcript vs the gene's genomic sequence
  - every candidate gets the same checks as the input primers (all transcripts, genomic region, genome with `--genome-check`); only clean, genomic-DNA-safe pairs are offered, ranked by isoform coverage, lab validation and design quality, with no overlapping primer positions

## Version 2.1.0 (2026-10-07)

### Output
- Results are saved next to the input file by default, not in the current directory
- `Details` split into `Problem` (what is wrong and what to do; blank for PASS) and `Notes` (information only); messages reworded in plain language, and the console shows the problem on its own line

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
