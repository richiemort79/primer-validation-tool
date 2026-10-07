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
    assert 55 < tm < 65                      # Wallace rule would give 64
    assert cp.calculate_tm('ACGTNACGTACGTACGTAC') is None


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


class _Opts:
    min_size, max_size, max_product = 70, 150, 4000
    max_mismatches, exact_3prime = 1, 5


def test_design_primer3_avoids_genomic_dna():
    transcript, genomic = _gene_model()
    exons = cp.map_exons(transcript, genomic)
    designs = cp.design_primer3(transcript, cp.exon_junctions(exons), _Opts)
    assert designs
    for c in designs:
        kind, how, _ = cp.intron_spanning(c.fwd, c.rev, transcript, exons, _Opts)
        assert kind in (0, 1), how
        assert not cp.evaluate_on_sequence(c.fwd, c.rev, genomic, 1, 5, 4000).products


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
