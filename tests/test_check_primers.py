"""Offline tests for the sequence, naming and input logic (no network)."""
import random
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_primers as cp  # noqa: E402


@pytest.mark.parametrize('name, expected', [
    ('Mitf_Fwd', ('Mitf', 'F', '')),
    ('Mitf_Rev', ('Mitf', 'R', '')),
    ('Tyr_EnH_rev', ('Tyr_EnH', 'R', '')),
    ('Cre_For', ('Cre', 'F', '')),
    ('Gapdh_F1', ('Gapdh', 'F', '1')),
    ('Gapdh-R2', ('Gapdh', 'R', '2')),
    ('Actb_FORWARD', ('Actb', 'F', '')),
    ('Rev3l_Fwd', ('Rev3l', 'F', '')),       # gene name contains "Rev"
    ('Forl_Rev', ('Forl', 'R', '')),
    ('Chr1_L_Flank', ('Chr1_L_Flank', None, '')),
    ('Cdkn2aR', ('Cdkn2aR', None, '')),      # no separator: not a direction
])
def test_parse_primer_name(name, expected):
    assert cp.parse_primer_name(name) == expected


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


def test_pairing():
    names = ['Mitf_Fwd', 'Mitf_Rev', 'Rev3l_Fwd', 'Rev3l_Rev', 'Gapdh_F1', 'Gapdh_R1',
             'Gapdh_F2', 'Gapdh_R2', 'Cre_For', 'ERT_Rev', 'Chr1_L_Flank']
    primers = [cp.Primer(n, 'ACGT' * 5, i, *cp.parse_primer_name(n)) for i, n in enumerate(names)]
    pairs, unpaired = cp.pair_primers(primers)
    assert {(f.name, r.name) for f, r in pairs} == {
        ('Mitf_Fwd', 'Mitf_Rev'), ('Rev3l_Fwd', 'Rev3l_Rev'),
        ('Gapdh_F1', 'Gapdh_R1'), ('Gapdh_F2', 'Gapdh_R2')}
    assert [p.name for p in unpaired] == ['Cre_For', 'ERT_Rev', 'Chr1_L_Flank']


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
