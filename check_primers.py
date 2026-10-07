#!/usr/bin/env python3
"""
Primer Validation Tool

Checks PCR primer pairs in silico:
  * gene-targeted pairs (e.g. Mitf_Fwd / Mitf_Rev) are tested against every
    RefSeq transcript of the gene and against the gene's genomic region (NCBI);
  * chromosome-locus primers (e.g. Chr1_L_Flank) are located on that chromosome
    (UCSC genome download, both strands) and the products they form are reported;
  * optionally (--genome-check) every primer is searched genome-wide to find
    off-target products.

Usage:
    python check_primers.py                          # Interactive file selection
    python check_primers.py primers.xlsx             # Specify file directly
    python check_primers.py primers.xlsx -o out.xlsx # Custom output name
    python check_primers.py primers.xlsx --genome-check
"""

import argparse
import gzip
import json
import os
import re
import sys
import time
from bisect import bisect_left
from dataclasses import dataclass, field, replace
from datetime import datetime
from pathlib import Path

import pandas as pd
import requests

try:
    import primer3
except ImportError:  # Tm is reported as blank if primer3-py is missing
    primer3 = None

try:
    import ahocorasick
except ImportError:  # only needed for genome searches
    ahocorasick = None


# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------

IUPAC = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT',
}
COMPLEMENT = str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN')

# Organism -> UCSC assembly. These match the assemblies NCBI Gene currently
# annotates, so NCBI gene coordinates and UCSC chromosome positions agree.
UCSC_DBS = {
    'mus musculus': 'mm39',
    'rattus norvegicus': 'rn7',
    'homo sapiens': 'hg38',
    'danio rerio': 'danRer11',
    'gallus gallus': 'galGal6',
    'drosophila melanogaster': 'dm6',
}
UCSC_DOWNLOAD = 'https://hgdownload.soe.ucsc.edu/goldenPath'

# Primer names are split into tokens on _ - . or space. A direction word
# (Fwd/For/Forward/Rev/Reverse, optionally numbered: Fwd2) may be any token
# after the first, so "Mitf_Fwd_IM" works; a single letter F/R (optionally
# numbered: F1) only counts as the last token, so "Chr1_R_Flank" is not
# read as a reverse primer.
NAME_SEP_RE = re.compile(r'([_\-\s.]+)')
DIRECTION_WORD_RE = re.compile(r'(?P<dir>fwd|for|forward|rev|reverse)(?P<num>\d*)', re.IGNORECASE)
DIRECTION_LETTER_RE = re.compile(r'(?P<dir>[fr])(?P<num>\d*)', re.IGNORECASE)
# Only an exact chromosome token counts, so genes such as Chrm1, Chrna7 or
# Chrd are still treated as genes.
CHROM_RE = re.compile(r'^chr(?P<chrom>\d+|X|Y|M|MT)$', re.IGNORECASE)

STATUS_ORDER = {'FAIL': 0, 'ERROR': 1, 'WARN': 2, 'PASS': 3, 'INFERRED': 4}

DEFAULT_CACHE = Path(os.environ.get(
    'PRIMER_TOOL_CACHE', Path.home() / '.cache' / 'primer_validation_tool'))


# --------------------------------------------------------------------------
# Sequence utilities
# --------------------------------------------------------------------------

def reverse_complement(seq):
    """Return reverse complement of a DNA sequence (IUPAC-aware)"""
    return seq.upper().translate(COMPLEMENT)[::-1]


def clean_sequence(raw):
    """Strip whitespace and inline modification codes (/5Phos/, [Phos]); RNA U -> T"""
    seq = re.sub(r'/[^/]*/', '', str(raw))
    seq = re.sub(r'\[[^\]]*\]', '', seq)
    return re.sub(r'\s+', '', seq).upper().replace('U', 'T')


def invalid_bases(seq):
    return sorted(set(seq) - set(IUPAC))


def calculate_tm(seq):
    """
    Nearest-neighbour Tm (SantaLucia 1998) via Primer3, using Primer3's default
    conditions: 50 mM monovalent, 1.5 mM Mg2+, 0.6 mM dNTPs, 50 nM oligo.
    Returns None for degenerate primers or if primer3-py is not installed.
    """
    if primer3 is None or not seq or set(seq) - set('ACGT'):
        return None
    return round(primer3.calc_tm(seq), 1)


def calculate_primer_stats(seq):
    """Length, GC content and Tm for a primer"""
    if not seq:
        return {'length': None, 'gc_content': None, 'tm': None}
    gc = seq.count('G') + seq.count('C') + seq.count('S')
    return {
        'length': len(seq),
        'gc_content': round(gc / len(seq) * 100, 1),
        'tm': calculate_tm(seq),
    }


@dataclass
class Site:
    """A primer binding site. strand '+' = primer sequence reads left-to-right
    in the target (extends rightwards); '-' = its reverse complement does."""
    start: int
    end: int
    strand: str
    mismatches: int
    primer: object = None
    chrom: str = None


@dataclass
class Product:
    left: Site   # '+' site
    right: Site  # '-' site
    size: int

    @property
    def chrom(self):
        return self.left.chrom

    def location(self):
        prefix = f'{self.chrom}:' if self.chrom else ''
        return f'{prefix}{self.left.start + 1}-{self.right.end}'


def _expand_degenerate(pattern, cap=64):
    """All concrete ACGT sequences a degenerate pattern stands for (None if > cap)"""
    combos = ['']
    for base in pattern:
        combos = [c + b for c in combos for b in IUPAC[base]]
        if len(combos) > cap:
            return None
    return combos


def _count_mismatches(allowed, target, limit):
    mm = 0
    for ok, base in zip(allowed, target):
        if base not in ok:
            mm += 1
            if mm > limit:
                break
    return mm


def _find_anchor(anchor, seq):
    """Start positions of an (optionally degenerate) anchor in seq, overlapping"""
    if not set(anchor) - set('ACGT'):
        pos = seq.find(anchor)
        while pos != -1:
            yield pos
            pos = seq.find(anchor, pos + 1)
    else:
        pattern = ''.join(f'[{IUPAC[b]}]' for b in anchor)
        for m in re.finditer(f'(?=({pattern}))', seq):
            yield m.start()


def find_primer_sites(primer, seq, max_mismatches=0, exact_3prime=5):
    """
    Find where a primer anneals in seq, on both strands.
    The 3'-most `exact_3prime` bases must match exactly (a 3' mismatch blocks
    extension); up to `max_mismatches` are tolerated elsewhere. Degenerate
    bases in the primer match any base they stand for.
    """
    L = len(primer)
    k = min(exact_3prime, L)
    sites = {}
    for strand, probe in (('+', primer), ('-', reverse_complement(primer))):
        allowed = [IUPAC[b] for b in probe]
        # The primer's 3' end is at the right of the probe on '+', left on '-'
        offset = L - k if strand == '+' else 0
        anchor = probe[offset:offset + k]
        for pos in _find_anchor(anchor, seq):
            start = pos - offset
            if start < 0 or start + L > len(seq):
                continue
            mm = _count_mismatches(allowed, seq[start:start + L], max_mismatches)
            if mm <= max_mismatches:
                sites[(start, strand)] = Site(start, start + L, strand, mm)
    return sorted(sites.values(), key=lambda s: (s.start, s.strand))


def find_products(sites, max_size):
    """
    All PCR products formed by a set of sites on one sequence: a '+' site
    followed, within max_size, by a '-' site. Sites may come from any primer,
    so primer-pair, self-primed and cross-pair products are all found.
    """
    plus = sorted((s for s in sites if s.strand == '+'), key=lambda s: s.start)
    minus = sorted((s for s in sites if s.strand == '-'), key=lambda s: s.start)
    minus_starts = [s.start for s in minus]
    products = []
    for a in plus:
        i = bisect_left(minus_starts, a.start)
        while i < len(minus) and minus[i].start < a.start + max_size:
            b = minus[i]
            size = b.end - a.start
            if b.end > a.end and size <= max_size:
                products.append(Product(a, b, size))
            i += 1
    return products


@dataclass
class SequenceResult:
    fwd_sites: list
    rev_sites: list
    products: list   # fwd/rev products (either orientation)
    extra: list      # products primed by a single primer (fwd+fwd or rev+rev)


def evaluate_on_sequence(fwd, rev, seq, max_mismatches, exact_3prime, max_size):
    """Locate both primers on seq and work out the products they form"""
    f_sites = find_primer_sites(fwd, seq, max_mismatches, exact_3prime)
    r_sites = find_primer_sites(rev, seq, max_mismatches, exact_3prime)
    for s in f_sites:
        s.primer = 'F'
    for s in r_sites:
        s.primer = 'R'
    products, extra = [], []
    for p in find_products(f_sites + r_sites, max_size):
        (products if p.left.primer != p.right.primer else extra).append(p)
    products.sort(key=lambda p: p.size)
    return SequenceResult(f_sites, r_sites, products, extra)


def product_mismatches(product):
    """(fwd mismatches, rev mismatches) for a fwd/rev product"""
    if product.left.primer == 'F':
        return product.left.mismatches, product.right.mismatches
    return product.right.mismatches, product.left.mismatches


# --------------------------------------------------------------------------
# Input
# --------------------------------------------------------------------------

@dataclass
class Primer:
    name: str
    seq: str
    row: int
    base: str = ''
    direction: str = None   # 'F', 'R' or None
    pair_num: str = ''
    tag: str = ''           # text after the direction, e.g. 'IM' in Mitf_Fwd_IM
    chrom: str = None       # e.g. 'chr1' for Chr1_* locus primers
    stats: dict = field(default_factory=dict)
    error: str = None
    paired_with: list = field(default_factory=list)
    binding: str = ''

    @property
    def pair_key(self):
        return (self.base.lower(), self.pair_num, self.tag.lower())

    @property
    def label(self):
        return self.base + self.pair_num


def parse_primer_name(name):
    """
    Split a primer name into (base, direction, pair number, tag).
    'Gapdh_Fwd' -> ('Gapdh', 'F', '', ''); 'Gapdh_R2' -> ('Gapdh', 'R', '2', '');
    'Mitf_Fwd_IM' -> ('Mitf', 'F', '', 'IM');
    'Chr1_L_Flank' -> ('Chr1_L_Flank', None, '', '')
    """
    name = str(name).strip()
    parts = NAME_SEP_RE.split(name)   # tokens at even indexes, separators between
    tokens = parts[0::2]
    found = None
    for k in range(len(tokens) - 1, 0, -1):
        m = DIRECTION_WORD_RE.fullmatch(tokens[k])
        if m:
            found = (k, m)
            break
    if found is None and len(tokens) > 1:
        m = DIRECTION_LETTER_RE.fullmatch(tokens[-1])
        if m:
            found = (len(tokens) - 1, m)
    if found is None:
        return name, None, '', ''
    k, m = found
    base = ''.join(parts[:2 * k - 1])
    tag = ''.join(parts[2 * k + 2:])
    direction = 'F' if m['dir'].lower().startswith('f') else 'R'
    return base, direction, m['num'], tag


def name_distance(a, b):
    """Edit distance counting a swap of neighbouring letters as one edit"""
    a, b = a.lower(), b.lower()
    d = [[max(i, j) if min(i, j) == 0 else 0 for j in range(len(b) + 1)] for i in range(len(a) + 1)]
    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            cost = a[i - 1] != b[j - 1]
            d[i][j] = min(d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + cost)
            if i > 1 and j > 1 and a[i - 1] == b[j - 2] and a[i - 2] == b[j - 1]:
                d[i][j] = min(d[i][j], d[i - 2][j - 2] + 1)
    return d[len(a)][len(b)]


def probable_typo(a, b):
    """Two different base names that are probably the same name mistyped"""
    return a.lower() != b.lower() and min(len(a), len(b)) >= 4 and name_distance(a, b) <= 2


def chromosome_of(base):
    """'Chr1_L_Flank' -> 'chr1'; 'ChrX' -> 'chrX'; 'Chrm1' -> None"""
    m = CHROM_RE.match(str(base).split('_')[0])
    if not m:
        return None
    chrom = m['chrom'].upper()
    return 'chr' + ('M' if chrom == 'MT' else chrom)


def gene_symbol_candidates(base):
    """Symbols to try in NCBI Gene: the full base name, then its first token"""
    candidates = [base]
    first = base.split('_')[0]
    if first and first != base:
        candidates.append(first)
    return candidates


def _find_columns(raw):
    """Locate the header row and the name/sequence columns"""
    for i in range(min(30, len(raw))):
        cells = [str(c).lower() if pd.notna(c) else '' for c in raw.iloc[i]]
        seq_col = next((j for j, c in enumerate(cells) if 'sequence' in c), None)
        name_col = next((j for j, c in enumerate(cells) if 'name' in c), None)
        if seq_col is not None and name_col is not None:
            return i, name_col, seq_col
    return None, 0, 2  # no header found: IDT layout, column A / column C


def load_primers(filepath, sheet=None):
    """
    Load primers from an Excel (IDT/Sigma order template) or CSV file.
    The header row is found by looking for 'Name' and 'Sequence' columns;
    without one, column A is taken as the name and column C as the sequence.
    """
    print(f"\nLoading primers from: {filepath}")
    if str(filepath).lower().endswith('.csv'):
        raw = pd.read_csv(filepath, header=None, dtype=str)
    else:
        raw = pd.read_excel(filepath, sheet_name=sheet if sheet is not None else 0,
                            header=None, dtype=str)

    header_row, name_col, seq_col = _find_columns(raw)
    first_row = 0 if header_row is None else header_row + 1

    primers = []
    for idx in range(first_row, len(raw)):
        raw_seq = raw.iat[idx, seq_col] if seq_col < raw.shape[1] else None
        if pd.isna(raw_seq) or not str(raw_seq).strip():
            continue
        seq = clean_sequence(raw_seq)
        if header_row is None and (len(seq) < 10 or invalid_bases(seq)):
            continue  # no header to anchor on: skip rows that aren't DNA
        raw_name = raw.iat[idx, name_col]
        name = str(raw_name).strip() if pd.notna(raw_name) and str(raw_name).strip() else f'Row{idx + 1}'
        base, direction, num, tag = parse_primer_name(name)
        primer = Primer(name=name, seq=seq, row=idx + 1, base=base, direction=direction,
                        pair_num=num, tag=tag, chrom=chromosome_of(base))
        bad = invalid_bases(seq)
        if bad:
            primer.error = f"Invalid characters in sequence: {''.join(bad)}"
        primer.stats = calculate_primer_stats(seq if not bad else '')
        primers.append(primer)

    print(f"✓ Found {len(primers)} primers")
    return primers


@dataclass
class PrimerPair:
    fwd: Primer
    rev: Primer
    how: str = 'name'   # 'name', 'adjacent' (neighbouring rows) or 'typo'
    note: str = ''


def pair_primers(primers):
    """
    Pair forward and reverse primers. Returns (pairs, unpaired).
      1. Same name apart from the direction (Mitf_Fwd_IM + Mitf_Rev_IM).
         Groups with several fwd or rev primers give every combination.
      2. Neighbouring rows holding one Fwd and one Rev (Cre_For + ERT_Rev),
         as order sheets usually list pairs on alternate lines.
      3. Names that differ by a likely typo (Tyrp1_Fwd + Tryp1_Rev).
    Only primers whose name has a direction are paired.
    """
    pairs = []

    def link(f, r, how='name', note=''):
        pairs.append(PrimerPair(f, r, how, note))
        f.paired_with.append(r.name)
        r.paired_with.append(f.name)

    groups = {}
    for p in primers:
        if p.direction:
            groups.setdefault(p.pair_key, []).append(p)
    for members in groups.values():
        for f in (p for p in members if p.direction == 'F'):
            for r in (p for p in members if p.direction == 'R'):
                link(f, r)

    def free(p):
        return p.direction and not p.paired_with

    i = 0
    while i < len(primers) - 1:
        a, b = primers[i], primers[i + 1]
        if free(a) and free(b) and {a.direction, b.direction} == {'F', 'R'}:
            f, r = (a, b) if a.direction == 'F' else (b, a)
            if probable_typo(f.base, r.base):
                link(f, r, 'typo', f"Primer names don't match: '{f.name}' vs '{r.name}'. "
                                   f"Probably a typo; fix the name before ordering")
            elif f.base.lower() == r.base.lower():
                link(f, r, 'adjacent', 'Paired because they are on neighbouring rows (name endings differ)')
            else:
                link(f, r, 'adjacent', f"Paired because they are on neighbouring rows "
                                       f"('{f.base}' and '{r.base}' differ)")
            i += 2
        else:
            i += 1

    for f in [p for p in primers if free(p) and p.direction == 'F']:
        candidates = [r for r in primers if free(r) and r.direction == 'R'
                      and (r.pair_num, r.tag.lower()) == (f.pair_num, f.tag.lower())
                      and probable_typo(f.base, r.base)]
        if candidates:
            r = min(candidates, key=lambda r: name_distance(f.base, r.base))
            link(f, r, 'typo', f"Primer names don't match: '{f.name}' vs '{r.name}'. "
                               f"Probably a typo; fix the name before ordering")

    unpaired = [p for p in primers if not p.paired_with]
    return pairs, unpaired


# --------------------------------------------------------------------------
# NCBI
# --------------------------------------------------------------------------

class NCBIError(Exception):
    pass


class NCBI:
    """Minimal E-utilities client with rate limiting and retries.
    Set NCBI_API_KEY to raise the limit from 3 to 10 requests/second."""
    BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'

    def __init__(self, api_key=None, email=None):
        self.session = requests.Session()
        self.api_key = api_key
        self.email = email
        self.interval = 0.11 if api_key else 0.34
        self._last = 0.0

    def request(self, endpoint, **params):
        params['tool'] = 'primer_validation_tool'
        if self.api_key:
            params['api_key'] = self.api_key
        if self.email:
            params['email'] = self.email
        error = None
        for attempt in range(5):
            wait = self.interval - (time.monotonic() - self._last)
            if wait > 0:
                time.sleep(wait)
            self._last = time.monotonic()
            try:
                # POST so long ID lists don't overflow the URL
                r = self.session.post(self.BASE + endpoint, data=params, timeout=120)
            except requests.RequestException as e:
                error = str(e)
            else:
                if r.status_code == 200:
                    return r
                error = f'HTTP {r.status_code}'
                if r.status_code not in (429, 500, 502, 503, 504):
                    break
            time.sleep(2 ** attempt)
        raise NCBIError(f'{endpoint} failed: {error}')

    def json(self, endpoint, **params):
        data = self.request(endpoint, retmode='json', **params).json()
        err = data.get('error') or data.get('esearchresult', {}).get('ERROR')
        if err:
            raise NCBIError(f'{endpoint}: {err}')
        return data


def parse_fasta(text):
    records, header, chunks = [], None, []
    for line in text.splitlines():
        if line.startswith('>'):
            if header is not None:
                records.append((header, ''.join(chunks).upper()))
            header, chunks = line[1:].strip(), []
        elif header is not None:
            chunks.append(line.strip())
    if header is not None:
        records.append((header, ''.join(chunks).upper()))
    return records


@dataclass
class GeneRecord:
    symbol: str
    gene_id: str
    description: str
    note: str                 # e.g. "matched alias of Xyz"
    chrom: str                # UCSC-style, e.g. 'chr6' (None if unplaced)
    chrom_start: int          # 0-based span of the gene on the chromosome
    chrom_end: int
    transcripts: list         # [(accession, sequence)], NM_ first
    genomic: str              # gene region +/- flank, in gene orientation
    genomic_flank: int

    def overlaps(self, chrom, start, end):
        return (self.chrom == chrom and self.chrom_start is not None
                and start < self.chrom_end + self.genomic_flank
                and end > self.chrom_start - self.genomic_flank)


def _accession_rank(acc):
    order = {'NM': 0, 'NR': 1, 'XM': 2, 'XR': 3}
    return (order.get(acc[:2], 4), acc)


def fetch_gene(ncbi, symbol, organism, flank=1000, cache_dir=None, refresh=False):
    """
    Resolve a gene symbol in NCBI Gene and fetch all its RefSeq transcripts
    plus its genomic region (+/- flank bp, gene orientation).
    Returns a GeneRecord, or None if no gene with that symbol (or alias) exists.
    Raises NCBIError on network/API failure.
    """
    cache_file = None
    if cache_dir:
        slug = re.sub(r'\W+', '_', organism.lower())
        cache_file = Path(cache_dir) / 'ncbi' / slug / f'{symbol.lower()}_f{flank}.json'
        if cache_file.exists() and not refresh:
            return GeneRecord(**json.loads(cache_file.read_text()))

    ids = ncbi.json('esearch.fcgi', db='gene', retmax=20,
                    term=f'"{symbol}"[Gene Name] AND "{organism}"[Organism] AND alive[prop]'
                    )['esearchresult'].get('idlist', [])
    if not ids:
        return None
    summary = ncbi.json('esummary.fcgi', db='gene', id=','.join(ids))['result']
    docs = [summary[i] for i in summary.get('uids', [])]

    note = ''
    match = next((d for d in docs if d.get('name', '').lower() == symbol.lower()), None)
    if match is None:
        alias_hits = [d for d in docs
                      if symbol.lower() in [a.strip().lower() for a in d.get('otheraliases', '').split(',')]]
        if not alias_hits:
            return None
        match = alias_hits[0]
        note = f"'{symbol}' is an alias of {match['name']}"
        if len(alias_hits) > 1:
            note += f" (also an alias of {', '.join(d['name'] for d in alias_hits[1:])})"

    gene_id = str(match['uid'])

    # Transcripts: exactly this gene's RefSeq RNAs
    links = ncbi.json('elink.fcgi', dbfrom='gene', db='nuccore', id=gene_id,
                      linkname='gene_nuccore_refseqrna')
    tx_ids = []
    for linkset in links.get('linksets', []):
        for lsdb in linkset.get('linksetdbs', []):
            tx_ids.extend(lsdb.get('links', []))
    transcripts = []
    if tx_ids:
        fasta = ncbi.request('efetch.fcgi', db='nuccore', id=','.join(map(str, tx_ids)),
                             rettype='fasta', retmode='text').text
        transcripts = sorted(((h.split()[0], s) for h, s in parse_fasta(fasta)),
                             key=lambda t: _accession_rank(t[0]))

    # Genomic region on the chromosome RefSeq
    chrom = chrom_start = chrom_end = None
    genomic = None
    info = (match.get('genomicinfo') or [None])[0]
    if info and info.get('chraccver'):
        a, b = int(info['chrstart']), int(info['chrstop'])
        chrom_start, chrom_end = min(a, b), max(a, b) + 1
        loc = str(info.get('chrloc', ''))
        chrom = 'chr' + ('M' if loc == 'MT' else loc) if loc else None
        fasta = ncbi.request('efetch.fcgi', db='nuccore', id=info['chraccver'],
                             rettype='fasta', retmode='text',
                             seq_start=max(1, chrom_start + 1 - flank),
                             seq_stop=chrom_end + flank,
                             strand=2 if a > b else 1).text
        records = parse_fasta(fasta)
        genomic = records[0][1] if records else None

    record = GeneRecord(symbol=match['name'], gene_id=gene_id,
                        description=match.get('description', ''), note=note,
                        chrom=chrom, chrom_start=chrom_start, chrom_end=chrom_end,
                        transcripts=transcripts, genomic=genomic, genomic_flank=flank)
    if cache_file:
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        cache_file.write_text(json.dumps(record.__dict__))
    return record


def spelling_suggestions(ncbi, symbol, organism):
    """
    Likely intended gene symbols for a misspelt one: any real gene that is
    the name with two neighbouring letters swapped (Endrb -> Ednrb), then
    NCBI's spelling suggestion (which is sometimes flaky or far off).
    """
    swaps = _transposed_symbols(ncbi, symbol, organism)
    try:
        text = ncbi.request('espell.fcgi', db='gene', term=symbol).text
        m = re.search(r'<CorrectedQuery>([^<]*)</CorrectedQuery>', text)
        suggestion = m.group(1).strip() if m else ''
        if suggestion and suggestion.lower() != symbol.lower() and re.fullmatch(r'[\w\-.]+', suggestion):
            return swaps + [suggestion]
    except NCBIError:
        pass
    return swaps


def _transposed_symbols(ncbi, symbol, organism):
    """Real gene symbols equal to `symbol` with two neighbouring letters swapped"""
    variants = {symbol[:i] + symbol[i + 1] + symbol[i] + symbol[i + 2:]
                for i in range(len(symbol) - 1)} - {symbol}
    if not variants:
        return []
    term = ('(' + ' OR '.join(f'"{v}"[Gene Name]' for v in sorted(variants))
            + f') AND "{organism}"[Organism] AND alive[prop]')
    try:
        ids = ncbi.json('esearch.fcgi', db='gene', retmax=20, term=term)['esearchresult'].get('idlist', [])
        if not ids:
            return []
        summary = ncbi.json('esummary.fcgi', db='gene', id=','.join(ids))['result']
    except NCBIError:
        return []
    wanted = {v.lower() for v in variants}
    names = [summary[i].get('name', '') for i in summary.get('uids', [])]
    return [n for n in names if n.lower() in wanted]


def resolve_gene(ncbi, bases, organism, gene_cache, opts):
    """
    Look up the gene for a primer pair's base name(s), trying each name, then
    its first token, then NCBI's spelling suggestion; memoised.
    Returns (GeneRecord or None, error message or None, warning or None).
    A warning is returned when the gene was only found via a spelling suggestion.
    """
    def lookup(symbol):
        key = symbol.lower()
        if key not in gene_cache:
            try:
                gene_cache[key] = (fetch_gene(ncbi, symbol, organism, opts.flank,
                                              opts.cache_dir, opts.refresh), None)
            except NCBIError as e:
                gene_cache[key] = (None, f'NCBI request failed: {e}')
        return gene_cache[key]

    errors = []
    for base in bases:
        for symbol in gene_symbol_candidates(base):
            gene, err = lookup(symbol)
            if gene:
                if symbol != base:
                    note = f"looked up as '{symbol}'"
                    gene = replace(gene, note=f'{gene.note}; {note}' if gene.note else note)
                return gene, None, None
            if err:
                errors.append(err)
    if errors:
        return None, errors[0], None

    for base in bases:
        for symbol in gene_symbol_candidates(base):
            for suggestion in spelling_suggestions(ncbi, symbol, organism):
                gene, _ = lookup(suggestion)
                # Only an official symbol: a suggestion that is merely an alias
                # (tryp -> Prss2) is too likely to be the wrong gene
                if gene and gene.symbol.lower() == suggestion.lower():
                    return gene, None, (f"There is no gene called '{symbol}', so the primers were "
                                        f"checked against {gene.symbol} instead. Probably a typo; "
                                        f"rename to {gene.symbol} if that is the intended gene")

    tried = "' / '".join(c for b in bases for c in gene_symbol_candidates(b))
    return None, (f"No {organism} gene called '{tried}' in NCBI Gene: check the name "
                  f"(transgenes such as Cre can't be checked)"), None


# --------------------------------------------------------------------------
# Genome (UCSC downloads)
# --------------------------------------------------------------------------

def _download(url, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + '.part')
    with requests.get(url, stream=True, timeout=120) as r:
        if r.status_code != 200:
            raise RuntimeError(f'download failed ({r.status_code}): {url}')
        total = int(r.headers.get('Content-Length', 0))
        done = 0
        with open(tmp, 'wb') as fh:
            for chunk in r.iter_content(1 << 20):
                fh.write(chunk)
                done += len(chunk)
                if total:
                    print(f'\r  Downloading {dest.name}: {done / 1e6:.0f}/{total / 1e6:.0f} MB',
                          end='', flush=True)
    print()
    tmp.rename(dest)


def genome_chromosomes(db, cache_dir):
    """Primary chromosome names for a UCSC assembly (no random/Un/alt contigs)"""
    path = Path(cache_dir) / 'genomes' / db / f'{db}.chrom.sizes'
    if not path.exists():
        _download(f'{UCSC_DOWNLOAD}/{db}/bigZips/{db}.chrom.sizes', path)
    names = [line.split('\t')[0] for line in path.read_text().splitlines() if line.strip()]
    return [n for n in names if re.fullmatch(r'chr(\d+|X|Y|M)', n)]


def load_chromosome(db, chrom, cache_dir):
    """Chromosome sequence (uppercase), downloaded once and cached on disk"""
    path = Path(cache_dir) / 'genomes' / db / f'{chrom}.fa.gz'
    if not path.exists():
        _download(f'{UCSC_DOWNLOAD}/{db}/chromosomes/{chrom}.fa.gz', path)
    with gzip.open(path, 'rt') as fh:
        fh.readline()  # FASTA header
        return fh.read().replace('\n', '').upper()


def scan_genome(primers, chroms, db, cache_dir, max_mismatches, anchor_len=15,
                max_sites=2000):
    """
    Find binding sites of every primer on the given chromosomes, both strands.
    Like UCSC In-Silico PCR, the 3'-most `anchor_len` bases must match exactly;
    up to `max_mismatches` are allowed in the rest of the primer.
    Returns {primer index: [Site]} and a set of primer indexes that hit
    max_sites (repetitive primers).
    """
    if ahocorasick is None:
        raise RuntimeError('pyahocorasick is required for genome searches: pip install pyahocorasick')

    automaton = ahocorasick.Automaton()
    keys = {}
    for i, p in enumerate(primers):
        L = len(p.seq)
        k = min(anchor_len, L)
        for strand, probe in (('+', p.seq), ('-', reverse_complement(p.seq))):
            offset = L - k if strand == '+' else 0
            allowed = [IUPAC[b] for b in probe]
            for anchor in _expand_degenerate(probe[offset:offset + k]) or []:
                keys.setdefault(anchor, []).append((i, strand, offset, L, allowed))
    if not keys:
        return {}, set()
    for anchor, value in keys.items():
        automaton.add_word(anchor, (len(anchor), value))
    automaton.make_automaton()

    sites = {i: [] for i in range(len(primers))}
    saturated = set()
    for chrom in chroms:
        seq = load_chromosome(db, chrom, cache_dir)
        print(f'  Scanning {db} {chrom}...', end=' ', flush=True)
        t0 = time.time()
        for end, (k, entries) in automaton.iter(seq):
            anchor_start = end - k + 1
            for i, strand, offset, L, allowed in entries:
                if i in saturated:
                    continue
                start = anchor_start - offset
                if start < 0 or start + L > len(seq):
                    continue
                mm = _count_mismatches(allowed, seq[start:start + L], max_mismatches)
                if mm <= max_mismatches:
                    sites[i].append(Site(start, start + L, strand, mm, primer=i, chrom=chrom))
                    if len(sites[i]) >= max_sites:
                        saturated.add(i)
        del seq
        print(f'{time.time() - t0:.0f}s')
    return sites, saturated


def genome_products(primer_ids, genome_sites, max_size):
    """Every product formed by any combination of the given primers, genome-wide"""
    by_chrom = {}
    for i in primer_ids:
        for s in genome_sites.get(i, []):
            by_chrom.setdefault(s.chrom, []).append(s)
    products = []
    for chrom_sites in by_chrom.values():
        products.extend(find_products(chrom_sites, max_size))
    return sorted(products, key=lambda p: p.size)


def describe_sites(sites, limit=5):
    if not sites:
        return ''
    shown = '; '.join(f"{s.chrom + ':' if s.chrom else ''}{s.start + 1}({s.strand})"
                      + (f' {s.mismatches}mm' if s.mismatches else '') for s in sites[:limit])
    more = f' ... +{len(sites) - limit} more' if len(sites) > limit else ''
    return f'{len(sites)} site(s): {shown}{more}'


def describe_products(products, names, limit=5):
    """'chr1:100-250 (151 bp, A+B); ...'"""
    parts = [f'{p.location()} ({p.size} bp, {names[p.left.primer]}+{names[p.right.primer]})'
             for p in products[:limit]]
    more = f' ... +{len(products) - limit} more' if len(products) > limit else ''
    return '; '.join(parts) + more


# --------------------------------------------------------------------------
# Evaluation
# --------------------------------------------------------------------------

def blank_result(fwd, rev, target):
    return {
        'Pair': f'{fwd.name} + {rev.name}',
        'Target': target,
        'Status': 'FAIL',
        'Problem': [],          # list while evaluating; joined by format_problems for output
        'Amplicon_Size_bp': None,
        'Notes': '',
        'Matched_Transcripts': '',
        'Transcripts_Amplified': '',
        'Genomic_Amplicon_bp': None,
        'Genome_Products': '',
        'Fwd_Name': fwd.name,
        'Rev_Name': rev.name,
        'Fwd_Primer': fwd.seq,
        'Rev_Primer': rev.seq,
        'Fwd_Length': fwd.stats.get('length'),
        'Rev_Length': rev.stats.get('length'),
        'Fwd_GC%': fwd.stats.get('gc_content'),
        'Rev_GC%': rev.stats.get('gc_content'),
        'Fwd_Tm_C': fwd.stats.get('tm'),
        'Rev_Tm_C': rev.stats.get('tm'),
        'Fwd_Mismatches': None,
        'Rev_Mismatches': None,
    }


def evaluate_gene_pair(fwd, rev, gene, opts, genome=None):
    """Check a fwd/rev pair against a gene's transcripts and genomic region"""
    result = blank_result(fwd, rev, gene.symbol)
    if gene.description:
        result['Target'] = f'{gene.symbol} ({gene.description})'
    notes, warnings = [], []
    if gene.note:
        notes.append(gene.note)

    def run(seq):
        return evaluate_on_sequence(fwd.seq, rev.seq, seq, opts.max_mismatches,
                                    opts.exact_3prime, opts.max_product)

    tx_results = [(acc, run(seq)) for acc, seq in gene.transcripts]
    gen = run(gene.genomic) if gene.genomic else None
    amplified = [(acc, r) for acc, r in tx_results if r.products]
    result['Transcripts_Amplified'] = f'{len(amplified)}/{len(tx_results)}'
    result['Matched_Transcripts'] = '; '.join(
        f"{acc} ({', '.join(str(p.size) for p in r.products)} bp)" for acc, r in amplified)
    if gen and gen.products:
        result['Genomic_Amplicon_bp'] = gen.products[0].size

    if amplified:
        acc, primary = amplified[0]
        product = primary.products[0]
        size = product.size
        result['Amplicon_Size_bp'] = size
        f_mm, r_mm = product_mismatches(product)
        result['Fwd_Mismatches'], result['Rev_Mismatches'] = f_mm, r_mm

        notes.insert(0, f'{size} bp product on {acc}')
        if not opts.min_size <= size <= opts.max_size:
            warnings.append(f'Product is {size} bp, outside the {opts.min_size}-{opts.max_size} bp target')
        for which, mm in (('Forward', f_mm), ('Reverse', r_mm)):
            if mm:
                warnings.append(f"{which} primer has {mm} mismatch{'es' if mm > 1 else ''} with {acc}: "
                                f'check the sequence (typo or strain difference?)')
        if len(primary.products) > 1:
            warnings.append(f'Makes {len(primary.products)} different products on {acc} ('
                            + ', '.join(f'{p.size}' for p in primary.products)
                            + ' bp): expect extra bands / melt peaks')
        if primary.extra:
            warnings.append(f'One primer on its own also amplifies {acc} ('
                            + ', '.join(f'{p.size}' for p in primary.extra)
                            + ' bp): may give extra products')
        if product.left.primer == 'R':
            notes.append('Fwd/Rev names are swapped relative to the transcript (works the same)')
        sizes = {p.products[0].size for _, p in amplified}
        if len(sizes) > 1:
            notes.append(f"Product size differs between isoforms ({', '.join(map(str, sorted(sizes)))} bp)")
        notes.append(f'Amplifies all {len(tx_results)} {gene.symbol} transcripts'
                     if len(amplified) == len(tx_results) and len(tx_results) > 1 else
                     f'Amplifies {len(amplified)} of {len(tx_results)} {gene.symbol} transcripts'
                     if len(tx_results) > 1 else f'{gene.symbol} has one RefSeq transcript')
        if gen is not None:
            if gen.products:
                notes.append(f'Does not span an intron: genomic DNA also gives a {gen.products[0].size} bp '
                             f'product (DNase-treat RNA)')
            else:
                notes.append('Spans an intron/exon junction (no product from genomic DNA)')
        result['Status'] = 'WARN' if warnings else 'PASS'
    elif gen and gen.products:
        product = gen.products[0]
        result['Amplicon_Size_bp'] = product.size
        result['Fwd_Mismatches'], result['Rev_Mismatches'] = product_mismatches(product)
        warnings.append(f"Doesn't amplify any {gene.symbol} mRNA ({len(tx_results)} RefSeq transcripts "
                        f'checked), only genomic DNA ({product.size} bp): it will not measure expression')
        if not opts.min_size <= product.size <= opts.max_size:
            warnings.append(f'Product is {product.size} bp, outside the {opts.min_size}-{opts.max_size} bp target')
        result['Status'] = 'WARN'
    else:
        all_results = [r for _, r in tx_results] + ([gen] if gen else [])
        f_binds = any(r.fwd_sites for r in all_results)
        r_binds = any(r.rev_sites for r in all_results)
        checked = (f'{len(tx_results)} {gene.symbol} transcripts and the {gene.symbol} genomic region '
                   f'±{opts.flank} bp')
        if not f_binds and not r_binds:
            warnings.append(f"Neither primer matches {gene.symbol} (checked {checked})")
        elif not f_binds:
            warnings.append(f"Forward primer doesn't match {gene.symbol} (checked {checked})")
        elif not r_binds:
            warnings.append(f"Reverse primer doesn't match {gene.symbol} (checked {checked})")
        else:
            warnings.append(f'Both primers match {gene.symbol} but cannot make a product together '
                            f'(facing the wrong way or more than {opts.max_product} bp apart)')
        if not gene.transcripts and not gene.genomic:
            result['Status'] = 'ERROR'
            warnings[:] = [f'NCBI has no sequences for {gene.symbol}']

    if genome is not None:
        products, names, saturated = genome
        on = [p for p in products if gene.overlaps(p.chrom, p.left.start, p.right.end)]
        off = [p for p in products if p not in on]
        result['Genome_Products'] = (f'{len(on)} on-target, {len(off)} off-target'
                                     + (f': {describe_products(off, names)}' if off else ''))
        if saturated:
            warnings.append('A primer matches too many places in the genome (repetitive sequence)')
        if off:
            same_size = result['Amplicon_Size_bp'] in {p.size for p in off}
            sites = ', '.join(f'{p.location()}, {p.size} bp' for p in off[:3])
            if len(off) > 3:
                sites += ' ...'
            warnings.append(f"Also amplifies {len(off)} other genomic site{'s' if len(off) > 1 else ''} "
                            f'({sites}): genomic DNA contamination could give a false signal'
                            + ('. Same size as the cDNA product, so probably a processed pseudogene'
                               if same_size else ''))
            if result['Status'] == 'PASS':
                result['Status'] = 'WARN'

    result['Problem'] = warnings
    result['Notes'] = '; '.join(notes)
    return result


def format_problems(problems):
    """One problem as-is; several numbered so each stands out"""
    if len(problems) <= 1:
        return ''.join(problems)
    return ' '.join(f'({i}) {p}.' for i, p in enumerate(problems, 1))


def evaluate_chrom_pair(fwd, rev, f_idx, r_idx, genome_sites, saturated, names, opts):
    """Check an explicit fwd/rev pair whose name refers to a chromosome locus"""
    result = blank_result(fwd, rev, fwd.chrom)
    products = genome_products([f_idx, r_idx], genome_sites, opts.max_product)
    pair_products = [p for p in products if p.left.primer != p.right.primer]
    on = [p for p in pair_products if p.chrom == fwd.chrom]
    off = [p for p in products if p not in on]
    if off:
        result['Genome_Products'] = f'{len(off)} other: {describe_products(off, names)}'

    if not on:
        f_here = [s for s in genome_sites.get(f_idx, []) if s.chrom == fwd.chrom]
        r_here = [s for s in genome_sites.get(r_idx, []) if s.chrom == fwd.chrom]
        if not f_here and not r_here:
            result['Problem'] = [f'Neither primer matches {fwd.chrom}']
        elif not f_here:
            result['Problem'] = [f"Forward primer doesn't match {fwd.chrom}"]
        elif not r_here:
            result['Problem'] = [f"Reverse primer doesn't match {fwd.chrom}"]
        else:
            result['Problem'] = [f'Both primers match {fwd.chrom} but cannot make a product together '
                                 f'(facing the wrong way or more than {opts.max_product} bp apart)']
        return result

    product = on[0]
    left_is_fwd = product.left.primer == f_idx
    f_mm = product.left.mismatches if left_is_fwd else product.right.mismatches
    r_mm = product.right.mismatches if left_is_fwd else product.left.mismatches
    result.update({'Amplicon_Size_bp': product.size, 'Fwd_Mismatches': f_mm, 'Rev_Mismatches': r_mm})
    warnings = []
    if len(on) > 1:
        warnings.append(f'Makes {len(on)} different products on {fwd.chrom} '
                        f'({describe_products(on, names)}): expect extra bands')
    for which, mm in (('Forward', f_mm), ('Reverse', r_mm)):
        if mm:
            warnings.append(f"{which} primer has {mm} mismatch{'es' if mm > 1 else ''} with the genome: "
                            f'check the sequence')
    if off:
        warnings.append(f'Also makes {len(off)} other product(s) elsewhere or from one primer alone '
                        f'({describe_products(off, names, 3)})')
    if {f_idx, r_idx} & saturated:
        warnings.append('A primer matches too many places in the genome (repetitive sequence)')
    result['Status'] = 'WARN' if warnings else 'PASS'
    result['Problem'] = warnings
    result['Notes'] = (f'{product.size} bp product at {product.location()} '
                       f'(qPCR size window not applied to chromosome loci)')
    return result


def inferred_rows(group_ids, primers, genome_sites, names, opts, existing_pairs, chrom=None):
    """Products formed by primers that have no named partner (e.g. Chr1_L_Flank + Chr1_L_Int)"""
    rows = []
    products = genome_products(group_ids, genome_sites, opts.max_product)
    for p in products:
        if chrom and p.chrom != chrom:
            continue
        a, b = primers[p.left.primer], primers[p.right.primer]
        if frozenset((a.name, b.name)) in existing_pairs:
            continue
        row = blank_result(a, b, p.chrom)
        row.update({
            'Pair': f'{a.name} + {b.name}' if a is not b else f'{a.name} (alone)',
            'Status': 'INFERRED',
            'Amplicon_Size_bp': p.size,
            'Fwd_Mismatches': p.left.mismatches,
            'Rev_Mismatches': p.right.mismatches,
            'Notes': (f'{p.size} bp product at {p.location()}, predicted from where the primers bind'
                      + ('; this one primer binds both ends' if a is b else '')),
        })
        rows.append(row)
    return rows


# --------------------------------------------------------------------------
# Main workflow
# --------------------------------------------------------------------------

def validate(primers, opts):
    """Run all checks; returns (pairs DataFrame, primers DataFrame)"""
    pairs, unpaired = pair_primers(primers)
    index = {id(p): i for i, p in enumerate(primers)}
    names = {i: p.name for i, p in enumerate(primers)}
    rows = []

    # ---- Genome scan (chromosome loci, and everything with --genome-check)
    genome_sites, saturated = {}, set()
    chrom_needed = {p.chrom for p in primers if p.chrom and not p.error}
    db = opts.genome_db or UCSC_DBS.get(opts.organism.lower())
    genome_error = None
    if chrom_needed or opts.genome_check:
        if not db:
            genome_error = (f"No UCSC assembly known for '{opts.organism}'; pass --genome-db "
                            f"(e.g. {', '.join(sorted(UCSC_DBS.values()))})")
        else:
            try:
                if opts.genome_check:
                    chroms = genome_chromosomes(db, opts.cache_dir)
                    scan_ids = [i for i, p in enumerate(primers) if not p.error]
                else:
                    chroms = sorted(chrom_needed, key=lambda c: (len(c), c))
                    scan_ids = [i for i, p in enumerate(primers) if p.chrom and not p.error]
                print(f"\n{'=' * 70}\nSearching {len(scan_ids)} primers on {db} "
                      f"({len(chroms)} chromosome{'s' if len(chroms) > 1 else ''})\n{'=' * 70}")
                sub = [primers[i] for i in scan_ids]
                sites, sat = scan_genome(sub, chroms, db, opts.cache_dir, opts.max_mismatches,
                                         opts.anchor)
                # Re-key by index in the full primer list
                for j, i in enumerate(scan_ids):
                    genome_sites[i] = sites.get(j, [])
                    for s in genome_sites[i]:
                        s.primer = i
                saturated = {scan_ids[j] for j in sat}
            except Exception as e:
                genome_error = f'Genome search failed: {e}'
                print(f'  ❌ {genome_error}')

    # ---- Pairs
    ncbi = NCBI(os.environ.get('NCBI_API_KEY'), os.environ.get('NCBI_EMAIL'))
    gene_cache = {}
    print(f"\n{'=' * 70}\nValidating {len(pairs)} primer pairs against {opts.organism}\n{'=' * 70}\n")
    for n, pair in enumerate(pairs, 1):
        fwd, rev = pair.fwd, pair.rev
        f_idx, r_idx = index[id(fwd)], index[id(rev)]
        print(f'[{n}/{len(pairs)}] {fwd.name} + {rev.name}')
        warning = None
        if fwd.error or rev.error:
            row = blank_result(fwd, rev, fwd.label)
            row.update(Status='ERROR', Problem=[e for e in (fwd.error, rev.error) if e])
        elif fwd.chrom and rev.chrom:
            if genome_error or db is None:
                row = blank_result(fwd, rev, fwd.chrom)
                row.update(Status='ERROR', Problem=[genome_error or 'Genome not searched'])
            else:
                row = evaluate_chrom_pair(fwd, rev, f_idx, r_idx, genome_sites, saturated, names, opts)
        else:
            # Only fall back to the reverse primer's name when it looks like a typo
            bases = [fwd.base] + ([rev.base] if pair.how == 'typo' else [])
            gene, err, warning = resolve_gene(ncbi, bases, opts.organism, gene_cache, opts)
            if gene is None:
                row = blank_result(fwd, rev, fwd.label)
                row.update(Status='ERROR', Problem=[err])
            else:
                genome = None
                if opts.genome_check and not genome_error:
                    genome = (genome_products([f_idx, r_idx], genome_sites, opts.max_product),
                              names, {f_idx, r_idx} & saturated)
                row = evaluate_gene_pair(fwd, rev, gene, opts, genome)
        # Naming problems go first: they are the thing to fix before ordering
        name_problems = [t for t in (warning, pair.note if pair.how == 'typo' else None) if t]
        if name_problems:
            row['Problem'] = name_problems + row['Problem']
            if row['Status'] == 'PASS':
                row['Status'] = 'WARN'
        if pair.how == 'adjacent':
            row['Notes'] = '; '.join(t for t in (pair.note, row['Notes']) if t)
        print_result(row)
        rows.append(row)

    # ---- Primers without a named partner
    existing = {frozenset((pair.fwd.name, pair.rev.name)) for pair in pairs}
    chrom_groups = {}
    for p in primers:
        if p.chrom and not p.error:
            chrom_groups.setdefault(p.chrom, []).append(index[id(p)])
    if not genome_error:
        for chrom, ids in chrom_groups.items():
            rows.extend(inferred_rows(ids, primers, genome_sites, names, opts, existing, chrom))
        loose = [index[id(p)] for p in unpaired if not p.chrom and not p.error]
        if opts.genome_check and loose:
            rows.extend(inferred_rows(loose, primers, genome_sites, names, opts, existing))

    for p in unpaired:
        if p.error:
            p.binding = p.error
        elif p.chrom or opts.genome_check:
            continue  # described from genome sites below
        else:
            gene, err, warning = resolve_gene(ncbi, [p.base], opts.organism, gene_cache, opts)
            if gene is None:
                p.binding = err
                continue
            if warning:
                gene = replace(gene, note=warning)
            parts = []
            for acc, seq in gene.transcripts:
                sites = find_primer_sites(p.seq, seq, opts.max_mismatches, opts.exact_3prime)
                if sites:
                    parts.append(f'{acc} {describe_sites(sites, 2)}')
                    break
            if gene.genomic:
                sites = find_primer_sites(p.seq, gene.genomic, opts.max_mismatches, opts.exact_3prime)
                parts.append(f'{gene.symbol} genomic region: {describe_sites(sites, 3)}' if sites
                             else f'does not bind {gene.symbol} genomic region (±{opts.flank} bp)')
            p.binding = '; '.join(([gene.note] if gene.note else []) + parts) or f'does not bind {gene.symbol}'
    for i, p in enumerate(primers):
        if genome_error and not p.binding and (p.chrom or opts.genome_check):
            p.binding = genome_error
        if i in genome_sites and not p.binding:
            p.binding = (describe_sites(genome_sites[i]) or f'no binding site in {db}'
                         if i not in saturated else 'repetitive: too many genome sites')

    for row in rows:
        row['Problem'] = format_problems(row['Problem'])
    pairs_df = pd.DataFrame(rows, columns=list(blank_result(Primer('', '', 0), Primer('', '', 0), '').keys()))
    if not pairs_df.empty:
        pairs_df['_order'] = pairs_df['Status'].map(STATUS_ORDER)
        pairs_df = pairs_df.sort_values('_order', kind='stable').drop(columns='_order')

    primers_df = pd.DataFrame([{
        'Row': p.row,
        'Name': p.name,
        'Sequence': p.seq,
        'Direction': {'F': 'Fwd', 'R': 'Rev'}.get(p.direction, ''),
        'Paired_With': ', '.join(p.paired_with) if p.paired_with else '(no partner)',
        'Length': p.stats.get('length'),
        'GC%': p.stats.get('gc_content'),
        'Tm_C': p.stats.get('tm'),
        'Binding': p.binding,
    } for p in primers])
    return pairs_df, primers_df, unpaired


def print_result(row):
    symbol = {'PASS': '✓', 'WARN': '⚠', 'FAIL': '✗', 'ERROR': '✗'}[row['Status']]
    print(f"  {symbol} {row['Status']}")
    if row['Problem']:
        print(f"      Problem: {format_problems(row['Problem'])}")
    if row['Notes']:
        print(f"      Notes:   {row['Notes']}")
    print()


def print_summary(pairs_df, unpaired, opts):
    print(f"\n{'=' * 70}\nVALIDATION SUMMARY\n{'=' * 70}\n")
    named = pairs_df[pairs_df['Status'] != 'INFERRED']
    total = len(named)
    counts = named['Status'].value_counts()
    labels = [
        ('PASS', '✓', f'PASS ({opts.min_size}-{opts.max_size} bp, clean)'),
        ('WARN', '⚠', 'WARN (product forms; see Problem)'),
        ('FAIL', '✗', 'FAIL (no product)'),
        ('ERROR', '✗', 'ERROR (lookup/input problem)'),
    ]
    print(f'Primer pairs checked: {total}')
    for status, sym, label in labels:
        c = int(counts.get(status, 0))
        pct = f' ({c / total * 100:.1f}%)' if total else ''
        print(f'{sym} {label:<36} {c}{pct}')
    inferred = int((pairs_df['Status'] == 'INFERRED').sum())
    if inferred:
        print(f'ℹ Inferred products from unpartnered primers: {inferred}')
    if unpaired:
        print(f"\nPrimers without a Fwd/Rev partner ({len(unpaired)}): "
              + ', '.join(p.name for p in unpaired))
        print('  (name pairs <Gene>_Fwd / <Gene>_Rev, or put each Fwd directly above its Rev)')

    for status, title in (('WARN', 'PAIRS WITH WARNINGS'), ('FAIL', 'PAIRS THAT FAILED'),
                          ('ERROR', 'PAIRS WITH ERRORS')):
        subset = named[named['Status'] == status]
        if len(subset):
            print(f"\n{'=' * 70}\n{title}:\n{'=' * 70}")
            for _, row in subset.iterrows():
                print(f"  {row['Pair']}")
                print(f"      {row['Problem']}")


def write_report(pairs_df, primers_df, output_file, opts):
    settings = pd.DataFrame([
        ('Run', datetime.now().strftime('%Y-%m-%d %H:%M')),
        ('Organism', opts.organism),
        ('Target amplicon size', f'{opts.min_size}-{opts.max_size} bp'),
        ('Max product size considered', f'{opts.max_product} bp'),
        ('Max mismatches per primer', opts.max_mismatches),
        ("3' bases that must match (gene check)", opts.exact_3prime),
        ("3' bases that must match (genome search)", opts.anchor),
        ('Genome-wide off-target check', 'yes' if opts.genome_check else 'no'),
        ('Genome assembly', opts.genome_db or UCSC_DBS.get(opts.organism.lower(), '')),
        ('Tm method', 'Primer3 nearest-neighbour (SantaLucia 1998); '
                      '50 mM Na+, 1.5 mM Mg2+, 0.6 mM dNTP, 50 nM oligo'),
    ], columns=['Setting', 'Value'])

    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        for name, df in (('Pairs', pairs_df), ('Primers', primers_df), ('Settings', settings)):
            df.to_excel(writer, index=False, sheet_name=name)
            ws = writer.sheets[name]
            ws.freeze_panes = 'A2'
            for col in ws.columns:
                width = max((len(str(c.value)) for c in col if c.value is not None), default=8)
                ws.column_dimensions[col[0].column_letter].width = min(max(width + 2, 8), 60)


def select_input_file():
    """Interactive file selection"""
    print("\n" + "=" * 70)
    print("PRIMER VALIDATION TOOL")
    print("=" * 70)
    print("\nPlease enter the path to your Excel file:")
    print("(You can drag and drop the file here, or type the path)")
    filepath = input("\nFile path: ").strip().strip('"').strip("'")
    if not os.path.exists(filepath):
        print(f"\n❌ Error: File not found: {filepath}")
        sys.exit(1)
    return filepath


def build_parser():
    parser = argparse.ArgumentParser(
        description='Validate PCR primer pairs against NCBI genes and the genome',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python check_primers.py                              # Interactive mode
  python check_primers.py primers.xlsx                 # Specify input file
  python check_primers.py primers.xlsx -o results.xlsx # Custom output name
  python check_primers.py primers.xlsx --genome-check  # Also look for off-target products

Environment:
  NCBI_API_KEY   raises NCBI's rate limit from 3 to 10 requests/second
  NCBI_EMAIL     contact address NCBI asks scripted clients to send
  PRIMER_TOOL_CACHE  cache directory (default ~/.cache/primer_validation_tool)
        """)
    parser.add_argument('input_file', nargs='?', help='Excel (.xlsx) or CSV file with primers')
    parser.add_argument('-o', '--output', help='Output Excel file (default: <input>_validation_results.xlsx in the input file\'s folder)')
    parser.add_argument('--sheet', help='Worksheet name (default: first sheet)')
    parser.add_argument('--organism', default='Mus musculus', help='Organism (default: Mus musculus)')
    parser.add_argument('--min-size', type=int, default=70, help='Minimum target amplicon size (default 70)')
    parser.add_argument('--max-size', type=int, default=150, help='Maximum target amplicon size (default 150)')
    parser.add_argument('--max-product', type=int, default=4000,
                        help='Largest product considered amplifiable, for gDNA/off-target/locus '
                             'products (default 4000)')
    parser.add_argument('--max-mismatches', type=int, default=1,
                        help="Mismatches tolerated per primer outside its 3' end (default 1); "
                             'any mismatch makes the pair WARN')
    parser.add_argument('--exact-3prime', type=int, default=5,
                        help="3' bases that must match exactly in gene checks (default 5)")
    parser.add_argument('--anchor', type=int, default=15,
                        help="3' bases that must match exactly in genome searches (default 15)")
    parser.add_argument('--flank', type=int, default=1000,
                        help='bp either side of the gene included in its genomic region (default 1000)')
    parser.add_argument('--genome-check', action='store_true',
                        help='Search the whole genome for off-target products '
                             '(first run downloads the genome, ~800 MB for mouse)')
    parser.add_argument('--genome-db', help='UCSC assembly to use (default: from --organism, e.g. mm39)')
    parser.add_argument('--cache-dir', default=str(DEFAULT_CACHE), help=f'Cache directory (default {DEFAULT_CACHE})')
    parser.add_argument('--refresh', action='store_true', help='Ignore cached NCBI gene data')
    return parser


def main():
    opts = build_parser().parse_args()

    input_file = opts.input_file or select_input_file()
    if not os.path.exists(input_file):
        print(f"\n❌ Error: File not found: {input_file}")
        sys.exit(1)
    # Default: next to the input file, wherever the script is run from
    output_file = opts.output or f"{os.path.splitext(os.path.abspath(input_file))[0]}_validation_results.xlsx"

    try:
        primers = load_primers(input_file, opts.sheet)
    except Exception as e:
        print(f"\n❌ Error loading primer file: {e}")
        sys.exit(1)
    if not primers:
        print("\n❌ No primer sequences found (expected 'Name' and 'Sequence' columns)")
        sys.exit(1)

    pairs_df, primers_df, unpaired = validate(primers, opts)
    print_summary(pairs_df, unpaired, opts)

    print(f"\n{'=' * 70}\nSaving results to: {output_file}")
    write_report(pairs_df, primers_df, output_file, opts)
    print(f"✓ Results saved (sheets: Pairs, Primers, Settings)\n{'=' * 70}\n")


if __name__ == "__main__":
    main()
