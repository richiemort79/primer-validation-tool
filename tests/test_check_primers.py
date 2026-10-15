"""Offline tests for the sequence, naming and input logic (no network)."""
import random
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_primers as cp  # noqa: E402


@pytest.mark.parametrize('name, expected', [
    ('Mitf_Fwd', ('Mitf', 'F', '', '')),
    ('Mitf_Rev', ('Mitf', 'R', '', '')),
    ('Tyr_EnH_rev', ('Tyr_EnH', 'R', '', '')),
    ('Cre_For', ('Cre', 'F', '', '')),
    ('Gapdh_F1', ('Gapdh', 'F', '1', '')),
    ('Gapdh-R2', ('Gapdh', 'R', '2', '')),
    ('Actb_FORWARD', ('Actb', 'F', '', '')),
    ('Rev3l_Fwd', ('Rev3l', 'F', '', '')),       # gene name contains "Rev"
    ('Forl_Rev', ('Forl', 'R', '', '')),
    ('Mitf_Fwd_IM', ('Mitf', 'F', '', 'IM')),    # direction mid-name
    ('Tryp1_Rev_IM', ('Tryp1', 'R', '', 'IM')),
    ('Sox10_Fwd2_qPCR_v1', ('Sox10', 'F', '2', 'qPCR_v1')),
    ('Chr1_L_Flank', ('Chr1_L_Flank', None, '', '')),
    ('Chr1_R_Flank', ('Chr1_R_Flank', None, '', '')),   # mid-name single R is not a direction
    ('Chr1_L_F', ('Chr1_L', 'F', '', '')),
    ('Cdkn2aR', ('Cdkn2aR', None, '', '')),      # no separator: not a direction
    ('Fwd', ('Fwd', None, '', '')),
])
def test_parse_primer_name(name, expected):
    assert cp.parse_primer_name(name) == expected


def test_name_distance_and_typo():
    assert cp.name_distance('Tyrp1', 'Tryp1') == 1   # transposition
    assert cp.probable_typo('Tyrp1', 'Tryp1')
    assert not cp.probable_typo('Mitf', 'mitf')
    assert not cp.probable_typo('Kit', 'Fit')        # too short to judge
    assert not cp.probable_typo('Snai1', 'Sox10')


@pytest.mark.parametrize('base, chrom', [
    ('Chr1_L_Flank', 'chr1'),
    ('chr19_site', 'chr19'),
    ('ChrX_del', 'chrX'),
    ('ChrMT_x', 'chrM'),
    ('Chrm1', None),      # cholinergic receptor genes are not chromosomes
    ('Chrna7', None),
    ('Chrd', None),
    ('Mitf', None),
])
def test_chromosome_of(base, chrom):
    assert cp.chromosome_of(base) == chrom


def test_clean_sequence():
    assert cp.clean_sequence(' /5Phos/acg tuA[Btn]\n') == 'ACGTTA'
    assert cp.invalid_bases('ACGTNRX') == ['X']


def test_reverse_complement_iupac():
    assert cp.reverse_complement('AACGTR') == 'YACGTT'


def test_tm_nearest_neighbour():
    tm = cp.calculate_tm('GGCTGTATTCCCCTCCATCG')
    assert 60 < tm < 67                      # mastermix conditions; 2+4 rule gives 64
    assert cp.calculate_tm('ACGTNACGTACGTACGTAC') is None


def test_tm_conditions_follow_mg_and_primer_conc():
    seq = 'ATGGCTGAAAGCAAGGAGTG'
    default = cp.calculate_tm(seq)
    try:
        cp.set_tm_conditions(1.5, 50)            # Primer3's own defaults
        assert cp.calculate_tm(seq) < default - 2
    finally:
        cp.set_tm_conditions(3.0, 250)


def test_wallace_tm():
    assert cp.wallace_tm('GGCTGTATTCCCCTCCATCG') == 64
    assert cp.wallace_tm('ACGTSW') == 2 * 3 + 4 * 3
    assert cp.wallace_tm('ACGTN') is None


@pytest.mark.parametrize('seq, expected', [
    ('GATTACAATTAC', 'both ends'),
    ('GATTACAATTAT', "5' only"),
    ('TATTACAATTAG', "3' only"),
    ('TATTACAATTAT', 'neither end'),
    ('GATTACAGCGCG', "both ends (3' end GC-rich: 5/5)"),
])
def test_clamp_status(seq, expected):
    assert cp.clamp_status(seq) == expected


def _random_seq(n, seed=1):
    rng = random.Random(seed)
    return ''.join(rng.choice('ACGT') for _ in range(n))


def test_sites_both_strands_and_product():
    seq = _random_seq(1000)
    fwd = seq[100:120]
    rev = cp.reverse_complement(seq[200:220])
    r = cp.evaluate_on_sequence(fwd, rev, seq, 0, 5, 4000)
    assert [(s.start, s.strand) for s in r.fwd_sites] == [(100, '+')]
    assert [(s.start, s.strand) for s in r.rev_sites] == [(200, '-')]
    assert [p.size for p in r.products] == [120]


def test_swapped_primers_still_amplify():
    seq = _random_seq(1000)
    fwd = seq[100:120]
    rev = cp.reverse_complement(seq[200:220])
    r = cp.evaluate_on_sequence(rev, fwd, seq, 0, 5, 4000)
    assert [p.size for p in r.products] == [120]
    assert r.products[0].left.primer == 'R'


def test_wrong_orientation_no_product():
    seq = _random_seq(1000)
    fwd = cp.reverse_complement(seq[100:120])
    rev = seq[200:220]
    r = cp.evaluate_on_sequence(fwd, rev, seq, 0, 5, 4000)
    assert r.fwd_sites and r.rev_sites and not r.products


def test_mismatch_tolerance_and_3prime_rule():
    seq = _random_seq(1000)
    site = seq[300:320]
    swap = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}
    five_prime_mm = swap[site[0]] + site[1:]
    three_prime_mm = site[:-1] + swap[site[-1]]
    assert cp.find_primer_sites(five_prime_mm, seq, 0, 5) == []
    hits = cp.find_primer_sites(five_prime_mm, seq, 1, 5)
    assert [(h.start, h.mismatches) for h in hits] == [(300, 1)]
    assert cp.find_primer_sites(three_prime_mm, seq, 1, 5) == []


def test_degenerate_primer():
    seq = _random_seq(1000)
    site = seq[400:420]
    degenerate = 'N' + site[1:10] + 'N' + site[11:]
    hits = cp.find_primer_sites(degenerate, seq, 0, 5)
    assert [(h.start, h.strand) for h in hits] == [(400, '+')]


def test_products_respect_max_size():
    seq = _random_seq(10000)
    fwd = seq[100:120]
    rev = cp.reverse_complement(seq[5000:5020])
    assert cp.evaluate_on_sequence(fwd, rev, seq, 0, 5, 4000).products == []
    assert [p.size for p in cp.evaluate_on_sequence(fwd, rev, seq, 0, 5, 6000).products] == [4920]


def _primers(names):
    return [cp.Primer(n, 'ACGT' * 5, i, *cp.parse_primer_name(n)) for i, n in enumerate(names)]


def _pairs(pairs):
    return {(p.fwd.name, p.rev.name, p.how) for p in pairs}


def test_pairing_by_name():
    primers = _primers(['Mitf_Fwd', 'Rev3l_Fwd', 'Mitf_Rev', 'Rev3l_Rev', 'Gapdh_F1', 'Gapdh_F2',
                        'Gapdh_R1', 'Gapdh_R2', 'Chr1_L_Flank', 'Chr1_R_Int'])
    pairs, unpaired = cp.pair_primers(primers)
    assert _pairs(pairs) == {
        ('Mitf_Fwd', 'Mitf_Rev', 'name'), ('Rev3l_Fwd', 'Rev3l_Rev', 'name'),
        ('Gapdh_F1', 'Gapdh_R1', 'name'), ('Gapdh_F2', 'Gapdh_R2', 'name')}
    assert [p.name for p in unpaired] == ['Chr1_L_Flank', 'Chr1_R_Int']


def test_pairing_melanoblast_sheet():
    primers = _primers(['Mitf_Fwd_IM', 'Mitf_Rev_IM', 'Tyrp1_Fwd_IM', 'Tryp1_Rev_IM',
                        'Endrb_Fwd_IM', 'Endrb_Rev_IM'])
    pairs, unpaired = cp.pair_primers(primers)
    assert _pairs(pairs) == {('Mitf_Fwd_IM', 'Mitf_Rev_IM', 'name'),
                             ('Tyrp1_Fwd_IM', 'Tryp1_Rev_IM', 'typo'),
                             ('Endrb_Fwd_IM', 'Endrb_Rev_IM', 'name')}
    assert unpaired == []


def test_pairing_adjacent_rows():
    primers = _primers(['Cre_For', 'ERT_Rev', 'Tyr_Pro_For', 'Tyr_EnH_rev',
                        'Chr1_L_Flank', 'Chr1_L_Int', 'Lone_Fwd', 'Other_Fwd'])
    pairs, unpaired = cp.pair_primers(primers)
    assert _pairs(pairs) == {('Cre_For', 'ERT_Rev', 'adjacent'),
                             ('Tyr_Pro_For', 'Tyr_EnH_rev', 'adjacent')}
    assert [p.name for p in unpaired] == ['Chr1_L_Flank', 'Chr1_L_Int', 'Lone_Fwd', 'Other_Fwd']


def test_typo_pairing_not_adjacent():
    primers = _primers(['Tyrp1_Fwd', 'Mitf_Fwd', 'Mitf_Rev', 'Tryp1_Rev'])
    pairs, _ = cp.pair_primers(primers)
    assert ('Tyrp1_Fwd', 'Tryp1_Rev', 'typo') in _pairs(pairs)


def test_load_primers_finds_header(tmp_path):
    rows = [[None, None, None],
            ['WHEN COPYING...', None, None],
            ['Oligo Name', "5' Mod", 'Sequence'],
            ['Mitf_Fwd', None, 'ccaggtgccgatggaagtc'],
            ['Mitf_Rev', None, 'GGTGAGCTCAGGACTTGGC '],
            [None, None, None],
            ['Bad_Fwd', None, 'ACGTXACGT']]
    path = tmp_path / 'p.xlsx'
    pd.DataFrame(rows).to_excel(path, header=False, index=False)
    primers = cp.load_primers(path)
    assert [(p.name, p.row, p.seq) for p in primers][:2] == [
        ('Mitf_Fwd', 4, 'CCAGGTGCCGATGGAAGTC'), ('Mitf_Rev', 5, 'GGTGAGCTCAGGACTTGGC')]
    assert primers[2].error and 'X' in primers[2].error


def test_load_primers_few_columns(tmp_path):
    path = tmp_path / 'p.csv'
    path.write_text('Name,Sequence\nA_Fwd,ACGTACGTACGT\nA_Rev,TTTTACGTACGT\n')
    assert [p.name for p in cp.load_primers(path)] == ['A_Fwd', 'A_Rev']


def test_load_primers_isoform_column(tmp_path):
    path = tmp_path / 'p.csv'
    path.write_text('Name,Sequence,Isoform\n'
                    'Mitf_Fwd,ACGTACGTACGT,"NM_008601.4, nm_001410059"\n'
                    'Mitf_Rev,TTTTACGTACGT,\n'
                    'Tyr_Fwd,ACGTACGTACGT,\n')
    primers = cp.load_primers(path)
    assert [p.isoforms for p in primers] == [['NM_008601.4', 'NM_001410059'], [], []]
    pairs, _ = cp.pair_primers(primers)
    assert pairs[0].isoforms == ['NM_008601.4', 'NM_001410059']   # either primer's cell counts


# --- Isoforms: catch-all and isoform-specific -----------------------------------

def _isoform_gene():
    """Synthetic gene: exon A or B, then shared exons C and D. Two curated
    transcripts use exon A (NM_1, NM_2), one exon B (NM_3); XM_4 skips exon C."""
    a, b, c, d = (_random_seq(n, seed) for n, seed in ((150, 21), (150, 22), (200, 23), (200, 24)))
    tx = [('NM_1.1', a + c + d), ('NM_2.1', a + c + d + 'ACGT'), ('NM_3.1', b + c + d),
          ('XM_4.1', a + d)]
    gene = cp.GeneRecord('Gx', '1', '', '', None, None, None, tx, None, 1000)
    return gene, (a, b, c, d)


def test_isoform_coverage_catch_all():
    gene, _ = _isoform_gene()
    problems, notes, score, _ = cp.isoform_coverage(gene, ['NM_1.1', 'NM_2.1', 'NM_3.1'])
    assert problems == [] and score == (0, 0)
    assert notes == ['Amplifies all 3 curated Gx transcripts']   # predicted XM_4 ignored
    problems, _, score, _ = cp.isoform_coverage(gene, ['NM_1.1', 'NM_2.1'])
    assert score == (1, 0) and 'Misses 1 of 3 curated Gx transcripts (NM_3.1)' in problems[0]


def test_isoform_coverage_specific():
    gene, _ = _isoform_gene()
    targets, unknown = cp.resolve_isoforms(gene, ['NM_1', 'NM_2.5', 'NM_9'])
    assert targets == ['NM_1.1', 'NM_2.1'] and unknown == ['NM_9']   # version ignored
    problems, notes, score, _ = cp.isoform_coverage(gene, ['NM_1.1', 'NM_2.1', 'XM_4.1'], targets)
    assert problems == [] and score == (0, 0) and notes == ['Specific to NM_1.1, NM_2.1']
    problems, _, score, _ = cp.isoform_coverage(gene, ['NM_1.1', 'NM_3.1'], targets)
    assert score == (1, 1)
    assert "Doesn't amplify target isoform NM_2.1" in problems[0]
    assert 'also amplifies NM_3.1' in problems[1]


def test_shared_and_unique_regions():
    gene, (a, b, c, d) = _isoform_gene()
    seqs = dict(gene.transcripts)
    template = seqs['NM_1.1']
    # Catch-all: only exons C and D are in every curated transcript
    excluded = cp.mask_regions(cp.shared_mask(template, [seqs['NM_2.1'], seqs['NM_3.1']]), False)
    # (the boundary can move a few bases where exons A and B happen to end alike)
    [[start, length]] = excluded
    assert start == 0 and len(a) - 5 <= length <= len(a)
    # NM_1 vs NM_3: exon A is what tells them apart
    [[start, length]] = cp.unique_regions(template, [seqs['NM_3.1']])
    assert start == 0 and len(a) - 5 <= length <= len(a)


def test_evaluate_gene_pair_isoforms():
    gene, (a, b, c, d) = _isoform_gene()
    f = cp.Primer('Gx_Fwd', a[100:120], 0)                    # exon A: NM_1, NM_2, XM_4
    r = cp.Primer('Gx_Rev', cp.reverse_complement(d[20:40]), 0)
    for p in (f, r):
        p.stats = cp.calculate_primer_stats(p.seq)
    res = cp.evaluate_gene_pair(f, r, gene, _Opts)
    assert res['Isoforms'] == 'All (catch-all)' and res['Status'] == 'WARN'
    # Counts and lists are curated only: XM_4 is amplified but not shown
    assert res['Transcripts_Amplified'] == '2/3'
    assert 'XM_4' not in res['Matched_Transcripts']
    assert any('Misses 1 of 3' in t for t in res['Problem'])
    res = cp.evaluate_gene_pair(f, r, gene, _Opts, isoforms=['NM_1.1', 'NM_2.1'])
    assert res['Isoforms'] == 'NM_1.1, NM_2.1'
    assert not any('isoform' in t.lower() for t in res['Problem'])


def test_isoform_label_from_name():
    gene, _ = _isoform_gene()
    assert cp.isoform_label('Gx-A', gene) == 'A'
    assert cp.isoform_label('GX-a_Pro', gene) == 'a'
    for base in ('Gx', 'Gx_A', 'Gx_Pro', 'Other-A'):     # only a hyphen after the gene counts
        assert cp.isoform_label(base, gene) == ''
    hyphenated = cp.GeneRecord('H2-K1', '1', '', '', None, None, None, [], None, 1000)
    assert cp.isoform_label('H2-K1', hyphenated) == ''


def test_isoform_coverage_named_in_primer_name():
    gene, _ = _isoform_gene()
    problems, notes, score, _ = cp.isoform_coverage(gene, ['NM_1.1', 'NM_2.1', 'XM_4.1'], label='A')
    assert problems == [] and score == (0, 0)                 # no catch-all WARN
    assert 'Targets isoform Gx-A (from the name): amplifies 2 of 3 curated transcripts (NM_1.1, NM_2.1)' in notes[0]
    # Accessions from an Isoform column still take precedence
    problems, _, _, _ = cp.isoform_coverage(gene, ['NM_1.1'], ['NM_1.1', 'NM_2.1'], label='A')
    assert problems


def test_gene_symbol_candidates_hyphen():
    assert cp.gene_symbol_candidates('MITF-M') == ['MITF-M', 'MITF']
    assert cp.gene_symbol_candidates('Tyr_Pro') == ['Tyr_Pro', 'Tyr']


# --- Alternative primers (--suggest) -------------------------------------------

def _gene_model(seed=7):
    """Synthetic gene: three exons separated by two 1 kb introns"""
    exons = [_random_seq(200, seed), _random_seq(180, seed + 1), _random_seq(220, seed + 2)]
    introns = ['GT' + _random_seq(996, seed + 3) + 'AG', 'GT' + _random_seq(996, seed + 4) + 'AG']
    genomic = _random_seq(300, seed + 5) + exons[0] + introns[0] + exons[1] + introns[1] + exons[2]
    return ''.join(exons), genomic


def test_map_exons_and_junctions():
    transcript, genomic = _gene_model()
    exons = cp.map_exons(transcript, genomic)
    assert [(a, b) for a, b, _ in exons] == [(0, 200), (200, 380), (380, 600)]
    assert [j for j, _ in cp.exon_junctions(exons)] == [200, 380]
    assert all(abs(intron - 1000) <= 2 for _, intron in cp.exon_junctions(exons))


def test_map_exons_tolerates_snp():
    transcript, genomic = _gene_model()
    pos = genomic.index(transcript[50:70]) + 10                 # a base inside exon 1
    genomic = genomic[:pos] + ('A' if genomic[pos] != 'A' else 'C') + genomic[pos + 1:]
    assert [(a, b) for a, b, _ in cp.map_exons(transcript, genomic)] == [(0, 200), (200, 380), (380, 600)]


_Opts = cp.build_parser().parse_args([])   # the command-line defaults


def test_design_primer3_avoids_genomic_dna():
    transcript, genomic = _gene_model()
    exons = cp.map_exons(transcript, genomic)
    designs = cp.design_primer3(transcript, cp.exon_junctions(exons), _Opts)
    assert designs
    lo, _, hi = cp.design_tm_range(_Opts)
    for c in designs:
        for seq in (c.fwd, c.rev):
            assert cp.has_clamp(seq, 5) and cp.has_clamp(seq, 3), seq
            assert cp.end_gc(seq) <= cp.MAX_END_GC
            assert lo - 0.5 <= cp.calculate_tm(seq) <= hi + 0.5
        kind, how, _ = cp.intron_spanning(c.fwd, c.rev, transcript, exons, _Opts)
        assert kind in (0, 1), how
        assert not cp.evaluate_on_sequence(c.fwd, c.rev, genomic, 1, 5, 4000).products


def test_lengthen_pair_raises_low_tm_keeping_3prime_end():
    transcript, _ = _gene_model()
    # Short primers either side of the exon 1/2 junction: Tm well under 60 °C
    fwd, rev = transcript[150:164], cp.reverse_complement(transcript[230:244])
    assert cp.calculate_tm(fwd) < 60 and cp.calculate_tm(rev) < 60
    c = cp.lengthen_pair(fwd, rev, transcript, _Opts)
    assert c.source == 'Lengthened'
    assert c.fwd.endswith(fwd) and c.rev.endswith(rev)        # only the 5' ends grow
    assert transcript.find(c.fwd) == 150 - (len(c.fwd) - len(fwd))
    assert cp.reverse_complement(c.rev) in transcript
    min_tm = cp.design_tm_range(_Opts)[0]
    assert cp.calculate_tm(c.fwd) >= min_tm and cp.calculate_tm(c.rev) >= min_tm
    # Nothing to do for primers that are already fine
    assert cp.lengthen_pair(c.fwd, c.rev, transcript, _Opts) is None


def test_tm_and_clamp_checks():
    low = cp.Primer('A_Fwd', 'ATTATAATTTAC', 0)
    ok = cp.Primer('A_Rev', 'GCTCACCTCATTCTCCGAAGT', 0)
    for p in (low, ok):
        p.stats = cp.calculate_primer_stats(p.seq)
    problems, notes = cp.tm_and_clamp_checks(low, ok, _Opts)
    assert any('Forward primer Tm' in t for t in problems)
    assert any('Tms differ' in t for t in problems)
    problems, notes = cp.tm_and_clamp_checks(ok, ok, _Opts)
    assert problems == [] and notes == ["Fwd primer has no G/C at its 3' end",
                                        "Rev primer has no G/C at its 3' end"]


def test_intron_spanning_within_exon():
    transcript, genomic = _gene_model()
    exons = cp.map_exons(transcript, genomic)
    fwd, rev = transcript[20:40], cp.reverse_complement(transcript[120:140])
    assert cp.intron_spanning(fwd, rev, transcript, exons, _Opts)[0] == 2


PRIMERBANK_PAGE = """<html><body>
<p>Primer Pair 1 (Click here for cDNA and amplicon sequence) :</p>
<table><tr><td>PrimerBank ID</td><td>142353045c1</td></tr>
<tr><td>Amplicon Size</td><td>88</td></tr>
<tr><td>Forward Primer</td><td>GGTGTCCCAAGACAACTTGTAA</td><td>22</td></tr>
<tr><td>Reverse Primer</td><td>CTCTCCAGCAGTTAGACCCCT</td><td>21</td></tr></table>
<p>Primer Pair 2 (Click here for cDNA and amplicon sequence) :</p>
<p>Validation Results&nbsp;&nbsp; (Click here to view experimental validation data)</p>
<table><tr><td>PrimerBank ID</td><td>31981217a1</td></tr>
<tr><td>Forward Primer</td><td>GAGCTTCCTTCCCGTGCTT</td></tr>
<tr><td>Reverse Primer</td><td>TGCCTGTTCCAGGTTTTAGTTAC</td></tr></table>
</body></html>"""


def test_fetch_primerbank_parses_page(monkeypatch):
    class Response:
        text = PRIMERBANK_PAGE

        def raise_for_status(self):
            pass
    monkeypatch.setattr(cp.requests, 'post', lambda *a, **k: Response())
    pairs, err = cp.fetch_primerbank('Pmel', 'Mus musculus')
    assert err is None
    assert [(c.source_id, c.fwd, c.rev, c.validated) for c in pairs] == [
        ('142353045c1', 'GGTGTCCCAAGACAACTTGTAA', 'CTCTCCAGCAGTTAGACCCCT', False),
        ('31981217a1', 'GAGCTTCCTTCCCGTGCTT', 'TGCCTGTTCCAGGTTTTAGTTAC', True)]
    assert cp.fetch_primerbank('Pmel', 'Danio rerio') == ([], None)   # mouse/human only
