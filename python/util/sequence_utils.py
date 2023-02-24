import numpy as np


# The following annotations are present in called/aligned sequences
# N: No-call
# C: cytosine
# A: adenine
# G: guanine
# T: thymine
# I: insertion in an aligned sequence
# D: deletion in a called sequence
# V: masked variant in a called sequence

NUCLEOTIDE2IDX = {'N' : 0, 'C' : 1, 'A' : 2, 'G' : 3, 'T' : 4, 'I': 5, 'D': 6, 'V': 7, 'S': 8}
NUC2COMPLEMENT = dict(zip("ACTGNIDVS", "TGACNIDVS"))
NUCS = 'CAGT'

def reverse_complement(input_sequence):
    """Generate reverse complement of input sequence

    Args:
        input_sequence (string): Input nucleotide sequence (A,C,G,T or N)

    Returns:
        string: Reverse-complement of input sequence
    """
    return "".join(NUC2COMPLEMENT[base] for base in reversed(input_sequence))
assert reverse_complement("SVIDNGCAT") == "ATGCNDIVS"


def process_sequence_records(sequence_iter, max_read_length, num_reads, right_shift = False):
    """Process sequence iterable into numpy basecalls byte array. See NUCLEOTIDE2IDX for mapping of
    nucleotide to byte.

    Args:
        sequence_iter (iter<string>): Iterable of input sequences
        max_read_length (int): Length of prefix to extract from input sequences
        num_reads (int): Expected number of records in iterable
        right_shift (bool): if true, right justify the final matrix of reads

    Returns:
        np.array<uint8>: max_read_length x num_reads
    """

    result = np.zeros((num_reads, max_read_length), dtype=np.uint8)

    for read_idx, record in enumerate(sequence_iter):
        if record is None:
            continue
        offset = 0
        
        if right_shift:
            offset = max_read_length - len(record)
        for cycle_idx, nucleotide in enumerate(record):
            if cycle_idx < max_read_length:
                result[read_idx, offset + cycle_idx] = NUCLEOTIDE2IDX[nucleotide]
    return result

def byte_encode_seq(seq, max_read_length=150, right_shift=False):
    """Process a single sequenceinto numpy basecalls byte array. See NUCLEOTIDE2IDX for mapping of
    nucleotide to byte.

    Args:
        seq (string): Input sequence
        max_read_length (int): Length of prefix to extract from input sequences
        right_shift (bool): if true, right justify the final matrix of reads

    Returns:
        np.array<uint8>: max_read_length
    """
    
    read_length = len(seq) if right_shift is False else max_read_length

    offset = 0
    if right_shift:
        offset = read_length - len(seq)

    byte_seq_array = np.zeros(read_length, dtype=np.uint8)

    for cycle_idx, nucleotide in enumerate(seq):
        if cycle_idx < read_length:
            byte_seq_array[offset + cycle_idx] = NUCLEOTIDE2IDX[nucleotide]

    return byte_seq_array