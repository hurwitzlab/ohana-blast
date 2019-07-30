import io

import extractseqs


blast_output_text = """\
A\tHOT234_1_0200m_rep_c55158_2\t100.00\t147\t0\t0\t1\t147\t400\t546 2e-70\t272
A\tHOT234_1_0200m_c10096_4\t100.00\t147\t0\t0\t1\t147\t400\t546\t2e-70\t272
A\tHOT238_1c_0200m_c3_1\t100.00\t147\t0\t0\t1\t147\t1\t147\t2e-70\t272
A\tHOT238_1c_0200m_rep_c260499_1\t100.00\t147\t0\t0\t1\t147\t400\t546\t2e-70\t272
B\tHOT234_1_0200m_rep_c55158_2\t100.00\t147\t0\t0\t1\t147\t400\t546 2e-70\t272
B\tHOT234_1_0200m_c10096_4\t100.00\t147\t0\t0\t1\t147\t400\t546\t2e-70\t272
C\tHOT238_1c_0200m_c3_1\t100.00\t147\t0\t0\t1\t147\t1\t147\t2e-70\t272"""

fasta_text = """\
>A
GTGC
>B
ATGC
>C
ATGG"""


def test_get_blast_hits():
    blast_output_file = io.StringIO(blast_output_text)
    blast_hits, seqid_counts = extractseqs.get_blast_reference_hits(blast_output_file=blast_output_file, blast_output_row_limit=None)

    assert len(blast_hits) == 2
    assert len(blast_hits['HOT234_1_0200m']) == 2
    assert 'HOT234_1_0200m_rep_c55158_2' in blast_hits['HOT234_1_0200m']
    assert 'HOT234_1_0200m_c10096_4' in blast_hits['HOT234_1_0200m']
    assert len(blast_hits['HOT238_1c_0200m']) == 2
    assert 'HOT238_1c_0200m_c3_1' in blast_hits['HOT238_1c_0200m']
    assert 'HOT238_1c_0200m_rep_c260499_1' in blast_hits['HOT238_1c_0200m']

    assert len(seqid_counts) == 4
    assert seqid_counts['HOT234_1_0200m_rep_c55158_2'] == 2
    assert seqid_counts['HOT234_1_0200m_c10096_4'] == 2
    assert seqid_counts['HOT238_1c_0200m_c3_1'] == 2
    assert seqid_counts['HOT238_1c_0200m_rep_c260499_1'] == 1


def test_get_blast_hits__limit():
    blast_output_file = io.StringIO(blast_output_text)
    blast_hits, seqid_counts = extractseqs.get_blast_reference_hits(blast_output_file=blast_output_file, blast_output_row_limit=3)

    assert len(blast_hits) == 2
    assert len(blast_hits['HOT234_1_0200m']) == 2
    assert 'HOT234_1_0200m_rep_c55158_2' in blast_hits['HOT234_1_0200m']
    assert 'HOT234_1_0200m_c10096_4' in blast_hits['HOT234_1_0200m']
    assert len(blast_hits['HOT238_1c_0200m']) == 1
    assert 'HOT238_1c_0200m_c3_1' in blast_hits['HOT238_1c_0200m']

    assert len(seqid_counts) == 3
    assert seqid_counts['HOT234_1_0200m_rep_c55158_2'] == 1
    assert seqid_counts['HOT234_1_0200m_c10096_4'] == 1
    assert seqid_counts['HOT238_1c_0200m_c3_1'] == 1


def test_find_sequences():
    fasta_file = io.StringIO(fasta_text)
    search_results = list(extractseqs.find_sequences(['C', 'B', 'A'], fasta_file=fasta_file))
    assert len(search_results) == 3
    assert search_results[0].id == 'A'
    assert search_results[1].id == 'B'
    assert search_results[2].id == 'C'


def test_find_sequences__first():
    fasta_file = io.StringIO(fasta_text)
    search_results = list(extractseqs.find_sequences(['A'], fasta_file=fasta_file))
    assert len(search_results) == 1
    assert search_results[0].id == 'A'


def test_find_sequences__middle():
    fasta_file = io.StringIO(fasta_text)
    search_results = list(extractseqs.find_sequences(['B'], fasta_file=fasta_file))
    assert len(search_results) == 1
    assert search_results[0].id == 'B'


def test_find_sequences__last():
    fasta_file = io.StringIO(fasta_text)
    search_results = list(extractseqs.find_sequences(['C'], fasta_file=fasta_file))
    assert len(search_results) == 1
    assert search_results[0].id == 'C'


def test_parse_blast_output_filename__contigs():
    blast_input_file_name, seq_type = extractseqs.parse_muscope_blast_output_filename('/some/dir/test_HOT224_1_0025m.fa-contigs.tab')
    assert blast_input_file_name == 'test_HOT224_1_0025m.fa'
    assert seq_type == 'contigs'


def test_parse_blast_output_filename__genes():
    blast_input_file_name, seq_type = extractseqs.parse_muscope_blast_output_filename('/some/dir/test_HOT224_1_0025m.fa-genes.tab')
    assert blast_input_file_name == 'test_HOT224_1_0025m.fa'
    assert seq_type == 'genes'


def test_parse_blast_output_filename__proteins():
    blast_input_file_name, seq_type = extractseqs.parse_muscope_blast_output_filename('/some/dir/test_HOT224_1_0025m.fa-proteins.tab')
    assert blast_input_file_name == 'test_HOT224_1_0025m.fa'
    assert seq_type == 'proteins'
