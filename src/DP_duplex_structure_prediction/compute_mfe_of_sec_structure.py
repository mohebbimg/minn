import os

mfes_dict = {}


def read_mfe_file_to_dict(fname):
    with open(fname, 'r') as infile:
        for line in infile.readlines():
            parts = line.split(":")
            mfes_dict[parts[0]] = float(parts[1])


def find_base_pairs(dot_bracket, sequence):
    """
    Find base pairs in the dot-bracket notation and return their indices and nucleotide types.
    """
    stack = []
    pairs = []

    # Iterate through each character in the dot-bracket string
    for i, char in enumerate(dot_bracket):
        if char == '(':  # Opening bracket
            stack.append(i)
        elif char == ')':  # Closing bracket
            if stack:
                j = stack.pop()  # Find the matching opening bracket
                pairs.append((j, i, sequence[j], sequence[i]))  # Store indices and the nucleotide types

    return pairs[::-1]


mfe_gen_rules = {
    "IMI": 4.09,  # intermolecular initiation
    "AU_p": 0.45,  # AU end penalty
    "SYM": 0.43,  # symmetry mfe, only if AU end penalty applied
    "Special_C_b":-0.9  # special bulge len 1, and it is C
}

curr_file_path = os.path.join(os.path.dirname(__file__), '')
read_mfe_file_to_dict(f'{curr_file_path}/mfe_files/stack.converted.txt')
read_mfe_file_to_dict(f'{curr_file_path}/mfe_files/int11.converted.txt')
read_mfe_file_to_dict(f'{curr_file_path}/mfe_files/int21.converted.txt')
read_mfe_file_to_dict(f'{curr_file_path}/mfe_files/int22.converted.txt')
read_mfe_file_to_dict(f'{curr_file_path}/mfe_files/loop_mfes.converted.txt')


def compute_MFE(seq1_5p3p_and_seq2_5p3p, sec_struct):
    sequences = seq1_5p3p_and_seq2_5p3p.replace('&', '')
    structure = sec_struct.replace('&', '')
    all_pairs = find_base_pairs(structure, sequences)

    # print("Base pairs:")
    bp_stacks = []
    bulges = []
    internal_loops = []
    computed_mfe = mfe_gen_rules["IMI"]  # initial mfe value of duplex sequences

    # for k in range(len(all_pairs)):
    #     i, j, nt1, nt2 = all_pairs[k]
    #     print(f"({i}, {j}) - {nt1}{nt2}")

    for k in range(len(all_pairs) - 1):
        tmp_tp1 = all_pairs[k]
        tmp_tp2 = all_pairs[k + 1]
        if abs(tmp_tp1[0] - tmp_tp2[0]) == 1 and abs(tmp_tp1[1] - tmp_tp2[1]) == 1:
            bp_stacks.append([tmp_tp1[2] + tmp_tp1[3] + "-" + tmp_tp2[2] + tmp_tp2[3], tmp_tp1, tmp_tp2])
        elif abs(tmp_tp1[0] - tmp_tp2[0]) == 1 or abs(tmp_tp1[1] - tmp_tp2[1]) == 1:
            bulges.append([tmp_tp1[2] + tmp_tp1[3] + "-" + tmp_tp2[2] + tmp_tp2[3], tmp_tp1, tmp_tp2])
        else:
            internal_loops.append([tmp_tp1[2] + tmp_tp1[3] + "-" + tmp_tp2[2] + tmp_tp2[3], tmp_tp1, tmp_tp2])

    for st in bp_stacks:
        computed_mfe += mfes_dict[st[0]]
        # print(st)

    # print(f"computed mfe: {computed_mfe}")

    # print("bulges:")
    for bl in bulges:
        surrounding_bps, bp1, bp2 = bl
        if abs(bp1[0] - bp2[0]) > 1:
            sub_seq = sequences[min(bp1[0], bp2[0]) + 1: max(bp1[0], bp2[0])]
        else:
            sub_seq = sequences[min(bp1[1], bp2[1]) + 1: max(bp1[1], bp2[1])]
        # print(bl, sub_seq)
        computed_mfe += mfes_dict[f"BULGE-{len(sub_seq)}"]
        if sub_seq == "C":
            computed_mfe += mfe_gen_rules["Special_C_b"]
        # print(f"computed mfe: {computed_mfe}")

    # print("internal loops:")
    for il in internal_loops:
        surrounding_bps, bp1, bp2 = il
        sub_seq1 = sequences[min(bp1[0], bp2[0]) + 1: max(bp1[0], bp2[0])]
        sub_seq2 = sequences[min(bp1[1], bp2[1]) + 1: max(bp1[1], bp2[1])]
        # print(il, sub_seq1, sub_seq2)

        if len(sub_seq1) == len(sub_seq2) == 2:
            big_X = sub_seq1[0] + sub_seq2[0]
            big_Y = sub_seq1[1] + sub_seq2[1]
            computed_mfe += mfes_dict[f"{surrounding_bps}-{big_X}-{big_Y}"]
        elif len(sub_seq1) == 1 and len(sub_seq2) == 2:
            big_X = sub_seq1
            YN = sub_seq2
            computed_mfe += mfes_dict[f"{surrounding_bps}-{big_X}-{YN}"]

        # print(f"computed mfe: {computed_mfe}")
    return computed_mfe


if __name__ == "__main__":
    # Input RNA sequences and their structures
    sequences = "UAAAGUAAAUAUGCACCAAAA&CGUUUUUAGCUUUGUGUUUACUUUU"
    structure = ".(((((((((((((.......&........))...)))))))))))."  # RNAcofold mfe (-10.10)

    # structure = ".((((..(((((((.......&........))...)))))..))))."  # changed

    computed_mfe = compute_MFE(sequences, structure)
    print(f"computed mfe: {computed_mfe}")
