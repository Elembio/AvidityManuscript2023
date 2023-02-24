import pysam

from .genome_mask import GENOME2DATA
from .sequence_utils import reverse_complement


def get_aligned_segment(called_seq, md_tag):
    """Use information in MD tag to convert from called sequence of segment to aligned sequence

    Args:
        called_seq (string): Called sequence of segment
        md_tag (string): MD tag (no insertions or deletions allowed)

    Returns:
        string: Called sequence, string: Aligned sequence
    """
    md_tag = md_tag.upper()
    called_list = list(called_seq)
    aligned_list = list(called_seq)
    digit_start = 0
    current_position = 0
    is_del = False

    for md_index in range(len(md_tag)):
        
        if is_del:
            if md_tag[md_index].isdigit():
                is_del = False
                digit_start = md_index
            else:
                called_list.insert(current_position, 'D')
                aligned_list.insert(current_position, md_tag[md_index])
                current_position += 1
        elif md_tag[md_index] == '^':
            if md_index > digit_start:
                current_position += int(md_tag[digit_start:md_index])
            is_del = True
        elif not md_tag[md_index].isdigit():
            nucleotide = md_tag[md_index]
            assert nucleotide in ["A", "C", "G", "T", "N", ], f'nucleotide {nucleotide} is unexpected'

            if md_index > digit_start:
                current_position += int(md_tag[digit_start:md_index])

            assert aligned_list[current_position] != nucleotide, f"{nucleotide} {current_position} {md_tag} {called_seq}"
            aligned_list[current_position] = nucleotide

            digit_start = md_index + 1
            current_position += 1

    return "".join(called_list), "".join(aligned_list)

assert get_aligned_segment("ACAT", "3G") == ("ACAT", "ACAG"), get_aligned_segment("ACAT", "3G")
assert get_aligned_segment("ACAT", "G3") == ("ACAT", "GCAT"), get_aligned_segment("ACAT", "G3")
assert get_aligned_segment("ACAT", "4") == ("ACAT", "ACAT")
assert get_aligned_segment("ACAT", "1GG1") == ("ACAT", "AGGT")
assert get_aligned_segment("ACAT", "1GGG") == ("ACAT", "AGGG")
assert get_aligned_segment("ACACA", "1G1G1") == ("ACACA", "AGAGA"), get_aligned_segment("ACACA", "1G1G1")
assert get_aligned_segment("ACAT", "2G^GC1") == ("ACADDT", "ACGGCT"), get_aligned_segment("ACAT", "2G^GC1")
assert get_aligned_segment("ACAT", "^A4") == ("DACAT", "AACAT")
assert get_aligned_segment("ACAT", "2^T1G") == ("ACDAT", "ACTAG")


def get_aligned_sequence(read, genome_mask = None, skip_ins=True):
    """Use MD tag and CIGAR string to build aligned sequence for read. Only handle S and M operations in CIGAR string (no indel or hard clip)

    Args:
        read (pysam.AlignedSegment): The input aligned read

    Returns:
        string: The aligned sequence
    """
    if not read.cigartuples:
        return None, None

    called_seq = read.seq

    current_start = 0
    preclip = 0
    for cigartuple in read.cigartuples:
        # skip, do not handle CINS
        if cigartuple[0] == pysam.CINS and not skip_ins:
            called_seq = called_seq[:current_start] + called_seq[current_start+cigartuple[1]:]
        elif cigartuple[0] == pysam.CMATCH:
            current_start += cigartuple[1]
        elif cigartuple[0] == pysam.CSOFT_CLIP:
            if current_start == 0:
                preclip = cigartuple[1]
            called_seq = called_seq[:current_start] + ("N" * cigartuple[1]) + called_seq[current_start+cigartuple[1]:]
            current_start += cigartuple[1]
        elif cigartuple[0] == pysam.CDEL:
            pass
        else:
            return None, None

    called_seq, aligned_seq = get_aligned_segment(called_seq[preclip:], read.get_tag("MD"))
    
    if genome_mask is not None:
        reference_positions = read.get_reference_positions(full_length=True)
        reference_name = read.reference_name

        # remove any chr outside of 1-22
        if reference_name not in GENOME2DATA["hg38"]['CONTIG2POSITION']:
            print(f"OOB chr name, {reference_name} not in {GENOME2DATA['hg38']['CONTIG2POSITION']}")
            return None, None
        else:
            masked_sequence = ''

            # fill in gaps (del or softclip) missing in in reference_positions
            padded_ref_positions = []
            last_pos = None
            for pos in reference_positions:
                if pos is not None and last_pos is not None:
                    temp_pos = last_pos + 1
                    while temp_pos != pos:
                        padded_ref_positions.append(temp_pos)
                        temp_pos += 1
                padded_ref_positions.append(pos)
                last_pos = pos

            # adjust for soft-clipping
            padded_ref_positions = padded_ref_positions[preclip:]

            assert len(called_seq) == len(padded_ref_positions), f'called-len [{len(called_seq)}] != padded-ref-len [{len(padded_ref_positions)}]'

            for i, (ref_pos, called_base) in enumerate(zip(padded_ref_positions, called_seq)):
                
                if ref_pos is None or genome_mask.slice(reference_name, ref_pos, 1)[0]:
                    masked_sequence += 'V'
                else:
                    masked_sequence += called_base

        called_seq = masked_sequence

    return called_seq, aligned_seq


def get_called_and_aligned(read, genome_mask: str = None, skip_ins=False) -> tuple:
    """Retrieve records from BAM file. Sequences are returned on original read strand. 

    Args:
        read (pysam.AlignedSegment) : BAM/SAM record
        genome_mask : SNP variant mask

    Returns:
        (list<string>,list<string>): tuple of (list of called sequences, list of aligned sequences)
    """

    called_seq, aligned_seq = get_aligned_sequence(read, genome_mask, skip_ins=skip_ins)
    if not called_seq:
        return None, None

    if read.is_reverse:
        called_seq = reverse_complement(called_seq)
        aligned_seq = reverse_complement(aligned_seq)

    assert len(called_seq) == len(aligned_seq), f'{len(called_seq)} != {len(aligned_seq)}\n{called_seq}\n{aligned_seq}'

    return called_seq, aligned_seq


