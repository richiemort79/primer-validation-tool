# Primer Validation Tool

A Python script that checks PCR primers in silico before you order them or run them:

- **Gene-targeted pairs** (e.g. `Mitf_Fwd` / `Mitf_Rev`) are tested against **every RefSeq transcript** of the gene and against the gene's **genomic region**, using NCBI. You get the amplicon size, which isoforms it amplifies, and whether it spans an intron (no genomic DNA product).
- **Chromosome-locus primers** (e.g. `Chr1_L_Flank`, `Chr1_R_Int`) are located on that chromosome, on both strands. Every product they can form with each other is reported, so you don't need Fwd/Rev names for genotyping primer sets.
- **Optional genome-wide check** (`--genome-check`) searches the whole genome for off-target products, similar to UCSC In-Silico PCR.

## Requirements

- Linux or macOS with Python 3.9 or higher
- An internet connection (NCBI E-utilities and the UCSC download server)

## Installation

```bash
./setup.sh
```

This creates a Python environment in `~/.venvs/primer_validation_tool` and installs the requirements. The environment is kept **outside the repo** on purpose: a venv synced by Dropbox/OneDrive breaks when it moves to another machine, so each machine gets its own.

You don't usually need to run setup yourself. The `./check_primers` wrapper runs it automatically the first time, and again whenever `requirements.txt` changes.

To put the environment somewhere else, set `PRIMER_TOOL_VENV=/path/to/venv` for both `setup.sh` and `check_primers`.

## Usage

```bash
./check_primers primers.xlsx                  # check pairs
./check_primers                               # prompts for the file
./check_primers primers.xlsx -o results.xlsx  # custom output name
./check_primers primers.xlsx --genome-check   # + genome-wide off-target search
./check_primers primers.xlsx --organism "Rattus norvegicus"
```

The wrapper works from any directory (e.g. `~/Dropbox/git-repos/primer_validation_tool/check_primers primers.xlsx`), and you can symlink it into `~/bin`. You can also call the script directly with `~/.venvs/primer_validation_tool/bin/python check_primers.py ...`.

### Options

| Option | Default | Meaning |
|--------|---------|---------|
| `--min-size` / `--max-size` | 70 / 150 | Target amplicon window for gene (qPCR) pairs |
| `--max-product` | 4000 | Largest product treated as amplifiable (genomic DNA, off-target and locus products) |
| `--max-mismatches` | 1 | Mismatches tolerated per primer. Any mismatch makes the pair WARN |
| `--exact-3prime` | 5 | 3′ bases that must match exactly (gene checks) |
| `--anchor` | 15 | 3′ bases that must match exactly (genome searches, as in UCSC isPCR) |
| `--flank` | 1000 | bp either side of the gene included in its genomic region |
| `--genome-check` | off | Genome-wide off-target search. The first run downloads the genome (~800 MB for mouse) |
| `--genome-db` | from organism | UCSC assembly (`mm39`, `rn7`, `hg38`, `danRer11`, `galGal6`, `dm6`) |
| `--sheet` | first sheet | Worksheet to read |
| `--refresh` | off | Re-download NCBI gene data instead of using the cache |
| `--cache-dir` | `~/.cache/primer_validation_tool` | Where NCBI data and genome files are cached |

Environment variables:

- `NCBI_API_KEY`: raises NCBI's rate limit from 3 to 10 requests/s.
- `NCBI_EMAIL`: contact address NCBI asks scripted clients to send.

## Input file

Use an Excel order template (IDT/Sigma) or a CSV. The script finds the header row by looking for columns containing **"Name"** and **"Sequence"**, so extra rows above the header are fine. If there is no header, it uses column A for names and column C for sequences.

Inline modification codes (`/5Phos/`, `[Btn]`), whitespace and lower case are cleaned up. Degenerate IUPAC bases (N, R, Y…) are supported. Any other character marks the primer as invalid.

### Naming and pairing

Each name is split into words on `_`, `-`, `.` or spaces. A direction word (`Fwd`, `For`, `Forward`, `Rev`, `Reverse`, case-insensitive, optionally numbered like `Fwd2`) can appear anywhere after the gene name. A single letter `F`/`R` only counts at the very end, so `Chr1_R_Flank` isn't read as a reverse primer.

| Name | Read as |
|------|---------|
| `Mitf_Fwd`, `Mitf_For`, `Mitf_F`, `Mitf_Forward` | forward primer for gene *Mitf* |
| `Mitf_Rev`, `Mitf_R`, `Mitf_reverse` | reverse primer for gene *Mitf* |
| `Mitf_Fwd_IM` + `Mitf_Rev_IM` | *Mitf* pair with tag `IM` (tags after the direction are fine) |
| `Gapdh_F1` + `Gapdh_R1`, `Gapdh_F2` + `Gapdh_R2` | two separate *Gapdh* pairs |
| `Tyr_Pro_For` | forward primer, looked up as `Tyr_Pro`, then `Tyr` |
| `Chr1_L_Flank`, `Chr1_R_Int` | chromosome 1 locus primers. Products are inferred between them |

Primers are paired in this order:

1. **Same name** apart from the direction (`Mitf_Fwd_IM` + `Mitf_Rev_IM`).
2. **Neighbouring rows**: an unpaired Fwd directly above or below an unpaired Rev, as order sheets usually list pairs on alternate lines (`Cre_For` + `ERT_Rev`). This is noted in *Notes*.
3. **Near-identical names** anywhere in the sheet (`Tyrp1_Fwd` + `Tryp1_Rev`). This is reported as a probable typo, and the pair is marked WARN so you fix the name before ordering.

If a gene name isn't found, the script tries swapping neighbouring letters and NCBI's spelling suggestion (`Endrb` → *Ednrb*). It only accepts a suggestion that is an official symbol. The pair is then checked against that gene and marked WARN.

Genes whose names merely start with "Chr" (`Chrm1`, `Chrna7`, `Chrd`) are treated as genes; only an exact chromosome token (`Chr1`, `chrX`, `ChrMT`) counts as a chromosome locus. A primer that still has no partner is listed on the **Primers** sheet with where it binds.

## Output

An Excel file with three sheets. By default it is saved next to the input file as `<input>_validation_results.xlsx`, whichever directory you run the script from. `-o` overrides this; a relative `-o` path is relative to your current directory.

### Pairs

| Column | Description |
|--------|-------------|
| **Pair** | `Fwd name + Rev name` |
| **Target** | Gene symbol and description, or chromosome |
| **Status** | PASS / WARN / FAIL / ERROR / INFERRED (see below) |
| **Problem** | What is wrong and what to do about it. Blank for PASS; several problems are numbered |
| **Amplicon_Size_bp** | Product size on the first amplified transcript (NM_ preferred) |
| **Notes** | Information only: product size, transcripts amplified, whether it spans an intron |
| **Matched_Transcripts** | Every RefSeq transcript amplified, with sizes |
| **Transcripts_Amplified** | e.g. `7/9` |
| **Genomic_Amplicon_bp** | Product on genomic DNA (blank = primers span an intron/junction) |
| **Genome_Products** | With `--genome-check`: on-target / off-target products and locations |
| **Fwd/Rev_Primer, _Length, _GC%, _Tm_C** | Primer properties |
| **Fwd/Rev_Mismatches** | Mismatches at the sites used for the product |

### Primers

One row per input primer: its source row, sequence, Tm, GC%, its partner(s), and where it binds. Binding is filled in for unpartnered primers, for chromosome-locus primers, and for every primer with `--genome-check`.

### Settings

The parameters used for the run, for your records.

### Status meanings

- **PASS**: product on the gene's transcript(s) within the size window. No mismatches, no extra products, no off-target products found.
- **WARN**: a product forms, but something needs a look. Possible reasons: size outside the window, a mismatch, several products, product only on genomic DNA, or off-target genome products. The *Problem* column says what to check.
- **FAIL**: no product. *Problem* says whether one primer doesn't bind, or both bind but can't form a product.
- **ERROR**: the gene couldn't be found or fetched, or a sequence has invalid characters.
- **INFERRED**: a product predicted between primers that weren't named as a pair (chromosome-locus sets; also unpartnered primers with `--genome-check`).

For **chromosome-locus** pairs the qPCR size window isn't applied. They PASS if there is exactly one product on that chromosome.

## How it works

1. **Gene lookup**: the gene is resolved in NCBI Gene by its exact symbol. If that fails, it tries an alias, which is reported in *Notes*. All of the gene's RefSeq RNAs (NM_/NR_/XM_/XR_) are fetched, plus its genomic region ± flank, in gene orientation. Results are cached, so re-runs are fast.
2. **Binding**: both strands are searched. The 3′-most bases must match exactly, and up to `--max-mismatches` are tolerated elsewhere.
3. **Products**: any forward-facing site followed by a reverse-facing site within `--max-product` bp. This also catches single-primer products and swapped Fwd/Rev names.
4. **Genome**: chromosomes are downloaded once from UCSC (`hgdownload`) and searched locally. This covers chromosome loci, and every chromosome with `--genome-check`.
5. **Tm**: nearest-neighbour model (SantaLucia 1998) via Primer3, at Primer3 default conditions (50 mM Na⁺, 1.5 mM Mg²⁺, 0.6 mM dNTPs, 50 nM oligo). Expect values a few °C different from IDT OligoAnalyzer, which uses different default conditions.

### Limitations

- The genome search is an exact-anchor search, not an alignment. Primers with a mismatch in their last 15 bases (`--anchor`) won't be found off-target. This is the same assumption UCSC In-Silico PCR makes.
- Transgene sequences (Cre, ERT2, reporters…) aren't in the reference genome or NCBI Gene, so pairs that amplify only a transgene can't be validated here.
- Off-target checking is for genomic DNA. Off-target products on other genes' *transcripts* aren't searched.

## Troubleshooting

- **"No Mus musculus gene named …"**: use the official gene symbol as the first part of the primer name (`Gapdh_Fwd`, not `GAPDH-mouse_Fwd`).
- **Primers listed as "without a Fwd/Rev partner"**: give them matching names with Fwd/Rev, or put each Fwd on the row directly above its Rev.
- **"NCBI request failed"**: a network or NCBI outage. Re-run; successfully fetched genes are cached.
- **Slow first `--genome-check`**: the genome is downloading into the cache directory. Later runs take about a minute.

## Tests

```bash
./setup.sh --dev                                   # installs pytest too
~/.venvs/primer_validation_tool/bin/python -m pytest tests
```

## License

MIT License - see LICENSE file for details
