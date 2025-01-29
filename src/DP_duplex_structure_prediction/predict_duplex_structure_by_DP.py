import csv
import os

import numpy as np

from src.DP_duplex_structure_prediction import compute_mfe_of_sec_structure


def load_base_pairing_scores(score_file):
    """Load the base pairing scores from a CSV file."""
    score_dict = {}
    all_weights = []
    valid_pairs = ['AU', 'UA', 'CG', 'GC', 'GU', 'UG']  # Define valid base pairs for single and combined pairs

    with open(score_file, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip header

        for row in reader:
            nt1_index, nt2_index, bps, score = int(row[0]), int(row[1]), row[2], float(row[3])
            all_weights.append(float(row[3]))

            # Check for single, double, or triple base pairs
            base_pairs = bps.split('-')

            # Store score in a nested dictionary structure for easier access
            if len(base_pairs) == 1:
                score_dict[(nt1_index, nt2_index, bps)] = score
            elif len(base_pairs) == 2:
                key = (nt1_index, nt2_index, '2bp', bps)
                score_dict[key] = score
            elif len(base_pairs) == 3:
                key = (nt1_index, nt2_index, '3bp', bps)
                score_dict[key] = score

    return score_dict, all_weights


def initialize_matrix(m, n):
    """Initialize the DP and traceback matrices for two sequences."""
    dp = [[0 for _ in range(n)] for _ in range(m)]
    traceback = [[None for _ in range(n)] for _ in range(m)]
    mx_mfe = [[0 for _ in range(n)] for _ in range(m)]
    return dp, traceback, mx_mfe


def rna_secondary_structure_two_sequences(microRNA_5p3p, cts_3p5p):
    score_dict = score_dict_g
    """Predict the secondary structure between two RNA sequences using dynamic programming and score maximization."""
    m, n = len(microRNA_5p3p), len(cts_3p5p)
    dp, traceback, mx_mfe = initialize_matrix(m, n)

    for i in range(1, m):
        for j in range(1, n):
            # Option 1: Ignore the i-th base of microRNA_5p3p or the j-th base of cts_3p5p
            if i > 0 and dp[i - 1][j] > dp[i][j]:
                dp[i][j] = dp[i - 1][j]
                traceback[i][j] = (i - 1, j)
                mx_mfe[i][j] = mx_mfe[i - 1][j]

            if j > 0 and dp[i][j - 1] > dp[i][j]:
                dp[i][j] = dp[i][j - 1]
                traceback[i][j] = (i, j - 1)
                mx_mfe[i][j] = mx_mfe[i][j - 1]

            # Handle single base pair
            try:
                single_key = (i, j, microRNA_5p3p[i] + cts_3p5p[j])
            except:
                print(microRNA_5p3p[i], cts_3p5p[j], microRNA_5p3p, cts_3p5p)
                raise

            if single_key in score_dict:
                score = dp[i - 1][j - 1] + score_dict[single_key]
                if score > dp[i][j]:
                    dp[i][j] = score
                    traceback[i][j] = (i - 1, j - 1)
                    two_seqs, dot_bracket = traceback_return_dot_bracket(traceback, microRNA_5p3p[:i + 1],
                                                                         cts_3p5p[:j + 1])
                    mx_mfe[i][j] = compute_mfe_of_sec_structure.compute_MFE(two_seqs, dot_bracket)

                    # Handle double and triple base pairs
            for k in range(2, 4):  # k = 2 for double, k = 3 for triple base pairs
                if i - k >= 0 and j - k >= 0:
                    tmp_pairs_arr = []
                    for r in range(k - 1, -1, -1):
                        tmp_pairs_arr.append(f'{microRNA_5p3p[i - r]}{cts_3p5p[j - r]}')
                    pair_key = '-'.join(tmp_pairs_arr)

                    key = (i - k, j - k, f'{k}bp', pair_key)
                    if key in score_dict:
                        score = dp[i - k][j - k] + score_dict[key]
                        if score > dp[i][j]:
                            dp[i][j] = score
                            traceback[i][j] = (i - k, j - k)
                            two_seqs, dot_bracket = traceback_return_dot_bracket(traceback, microRNA_5p3p[:i + 1],
                                                                                 cts_3p5p[:j + 1])
                            mx_mfe[i][j] = compute_mfe_of_sec_structure.compute_MFE(two_seqs, dot_bracket)

    return dp, traceback, mx_mfe


def traceback_print_bps(traceback, microRNA_5p3p, cts_3p5p):
    """Reconstruct and print the predicted RNA secondary structure based on traceback matrix."""
    i, j = len(microRNA_5p3p) - 1, len(cts_3p5p) - 1
    structure = []

    while i > 0 and j > 0:
        if traceback[i][j] is None:
            break

        prev_i, prev_j = traceback[i][j]
        if prev_i is None or prev_j is None:
            break

        if i - prev_i == 1 and j - prev_j == 1:
            # Single base pair match
            structure.append(f"({microRNA_5p3p[i]}-{cts_3p5p[j]})")
        elif i - prev_i == 2 and j - prev_j == 2:
            # Double base pair match
            structure.append(f"({microRNA_5p3p[prev_i + 1:i + 1]}-{cts_3p5p[prev_j + 1:j + 1]})")
        elif i - prev_i == 3 and j - prev_j == 3:
            # Triple base pair match
            structure.append(f"({microRNA_5p3p[prev_i + 1:i + 1]}-{cts_3p5p[prev_j + 1:j + 1]})")

        i, j = prev_i, prev_j

    structure.reverse()
    print("Predicted RNA Secondary Structure:")
    print(" -> ".join(structure))


def traceback_return_dot_bracket(traceback, microRNA_5p3p, cts_3p5p):
    """Generate dot-bracket notation for the predicted RNA secondary structure."""
    m, n = len(microRNA_5p3p), len(cts_3p5p)
    structure1 = ['.'] * m
    structure2 = ['.'] * n

    i, j = m - 1, n - 1

    while i > 0 and j > 0:
        if traceback[i][j] is None:
            break

        prev_i, prev_j = traceback[i][j]
        if prev_i is None or prev_j is None:
            break

        if i - prev_i == 1 and j - prev_j == 1:
            # Single base pair match
            structure1[i] = '('
            structure2[j] = ')'
        elif i - prev_i == 2 and j - prev_j == 2:
            # Double base pair match
            structure1[prev_i + 1:i + 1] = ['('] * (i - prev_i)
            structure2[prev_j + 1:j + 1] = [')'] * (j - prev_j)
        elif i - prev_i == 3 and j - prev_j == 3:
            # Triple base pair match
            structure1[prev_i + 1:i + 1] = ['('] * (i - prev_i)
            structure2[prev_j + 1:j + 1] = [')'] * (j - prev_j)

        i, j = prev_i, prev_j

    # Combine the two structures
    two_seqs = microRNA_5p3p + '&' + cts_3p5p[::-1]  # both seqs are 5p to 3p
    dot_bracket = ''.join(structure1) + '&' + ''.join(reversed(structure2))
    return two_seqs, dot_bracket


def find_min_MFE_indices(mfe_mx):
    n, m = len(mfe_mx), len(mfe_mx[0])
    min_mfe = None
    min_i, min_j = -1, -1
    for i in range(n):
        for j in range(m):
            if min_mfe is None or min_mfe >= mfe_mx[i][j]:
                min_mfe = mfe_mx[i][j]
                min_i, min_j = i, j
    return min_i, min_j


def compute_sample_score(microRNA_seq, cts_seq, all_weights):
    valid_pairs = sorted(['AU', 'UA', 'CG', 'GC', 'GU', 'UG'])
    max_miRNA_len, max_target_len = 25, 25

    miRNA = microRNA_seq
    target_candidate = cts_seq

    num_valid_pairs = len(valid_pairs)
    # MX1; one nt
    mx_of_one_nt = np.zeros((max_miRNA_len, max_target_len, num_valid_pairs))
    for i in range(min(len(miRNA), max_miRNA_len)):
        for j in range(min(len(target_candidate), max_target_len)):
            nts_pair = miRNA[i] + target_candidate[j]
            if nts_pair in valid_pairs:
                one_hot = valid_pairs.index(nts_pair)
                mx_of_one_nt[i, j, one_hot] = 1

    flattened_one_nt = mx_of_one_nt.reshape(-1)

    # MX2; two nts; each nt to interact with max 5 nts
    mx_of_two_nts = np.zeros((max_miRNA_len, 5, num_valid_pairs * num_valid_pairs))
    for i in range(min(len(miRNA), max_miRNA_len)):
        tmp_mx_col = 0
        for j in range(max(0, i - 2), min(len(target_candidate), i + 3)):  # +3 because last not included.
            if i + 1 >= len(miRNA) or j + 1 >= len(target_candidate):
                continue

            nts_pair1 = miRNA[i] + target_candidate[j]
            nts_pair2 = miRNA[i + 1] + target_candidate[j + 1]
            if nts_pair1 in valid_pairs and nts_pair2 in valid_pairs:
                one_hot1 = valid_pairs.index(nts_pair1)
                one_hot2 = valid_pairs.index(nts_pair2)
                mixed_onehot = one_hot1 * num_valid_pairs + one_hot2
                mx_of_two_nts[i][tmp_mx_col][mixed_onehot] = 1

            tmp_mx_col += 1

    flattened_two_nts = mx_of_two_nts.reshape(-1)

    # MX3; three nts; each nt to interact with max 3 nts
    mx_of_three_nts = np.zeros((max_miRNA_len, 3, num_valid_pairs * num_valid_pairs * num_valid_pairs))
    for i in range(min(len(miRNA), max_miRNA_len)):
        tmp_mx_col = 0
        for j in range(max(0, i - 1), min(len(target_candidate), i + 2)):  # +2 because last not included.
            if i + 2 >= len(miRNA) or j + 2 >= len(target_candidate):
                continue

            nts_pair1 = miRNA[i] + target_candidate[j]
            nts_pair2 = miRNA[i + 1] + target_candidate[j + 1]
            nts_pair3 = miRNA[i + 2] + target_candidate[j + 2]
            if nts_pair1 in valid_pairs and nts_pair2 in valid_pairs and nts_pair3 in valid_pairs:
                one_hot1 = valid_pairs.index(nts_pair1)
                one_hot2 = valid_pairs.index(nts_pair2)
                one_hot3 = valid_pairs.index(nts_pair3)
                mixed_onehot = one_hot1 * num_valid_pairs * num_valid_pairs + one_hot2 * num_valid_pairs + one_hot3
                mx_of_three_nts[i][tmp_mx_col][mixed_onehot] = 1

            tmp_mx_col += 1

    flattened_three_nts = mx_of_three_nts.reshape(-1)

    # MX4; three nts; each nt to interact with max 3 nts
    mx_of_four_nts = np.zeros((max_miRNA_len, 1, num_valid_pairs * num_valid_pairs * num_valid_pairs * num_valid_pairs))
    for i in range(min(len(miRNA), max_miRNA_len)):
        j = i
        if i + 3 >= len(miRNA) or j + 3 >= len(target_candidate):
            continue

        nts_pair1 = miRNA[i] + target_candidate[j]
        nts_pair2 = miRNA[i + 1] + target_candidate[j + 1]
        nts_pair3 = miRNA[i + 2] + target_candidate[j + 2]
        nts_pair4 = miRNA[i + 3] + target_candidate[j + 3]

        if nts_pair1 in valid_pairs and nts_pair2 in valid_pairs and nts_pair3 in valid_pairs and nts_pair4 in valid_pairs:
            one_hot1 = valid_pairs.index(nts_pair1)
            one_hot2 = valid_pairs.index(nts_pair2)
            one_hot3 = valid_pairs.index(nts_pair3)
            one_hot4 = valid_pairs.index(nts_pair4)

            mixed_onehot = one_hot1 * num_valid_pairs * num_valid_pairs * num_valid_pairs + \
                           one_hot2 * num_valid_pairs * num_valid_pairs + \
                           one_hot3 * num_valid_pairs + one_hot4

            mx_of_four_nts[i][0][mixed_onehot] = 1

    flattened_four_nts = mx_of_four_nts.reshape(-1)

    # sample_vector = np.concatenate((flattened_one_nt, flattened_two_nts, flattened_three_nts, flattened_four_nts))
    sample_vector = np.concatenate((flattened_one_nt, flattened_two_nts, flattened_three_nts))

    result = sum(x * y for x, y in zip(sample_vector, all_weights))
    sigmoid_value = 1 / (1 + np.exp(-result))
    return sigmoid_value


def print_matrix(matrix, name="Matrix"):
    """Prints the matrix in a readable format, handling None values."""
    print(f"\n{name}:")
    for row in matrix:
        print(" ".join(f"{val:.2f}" if isinstance(val, (int, float)) else str(val) for val in row))


def pred_structure_get_total_score(microRNA_5p3p, cts_5p3p):
    cts_3p5p = cts_5p3p[::-1]
    dp_matrix, traceback_matrix, mfe_mx = rna_secondary_structure_two_sequences(microRNA_5p3p, cts_3p5p)
    # test the whole structure
    two_seqs, struct_dot_bracket = traceback_return_dot_bracket(traceback_matrix, microRNA_5p3p, cts_3p5p)

    computed_mfe = compute_mfe_of_sec_structure.compute_MFE(two_seqs, struct_dot_bracket)

    # structure_score = dp_matrix[end_i][end_j]
    structure_score = dp_matrix[-1][-1]
    s_score_sigmoid = 1 / (1 + np.exp(-structure_score))

    ML_model_score = compute_sample_score(microRNA_5p3p, cts_3p5p, all_weights_g)

    # Print DP matrix for debugging
    # print_matrix(dp_matrix, name="DP Matrix")
    # print_matrix(traceback_matrix, name="Traceback Matrix")
    # print_matrix(mfe_mx, "MFE mx")

    return two_seqs, struct_dot_bracket, computed_mfe, s_score_sigmoid, ML_model_score, traceback_matrix, microRNA_5p3p, cts_3p5p

# Load the base pairing scores from the CSV file
file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'weights_for_DP_algorithm.csv')

score_dict_g, all_weights_g = load_base_pairing_scores(file_path)

#
# if __name__ == "__main__":
#     # Test with two sample RNA sequences
#     rna_sequence_1g = "CUAUCUCACGUCUGGUCCCAGA"[::-1]
#     rna_sequence_2g = "UUAAUUUAAUUAAAGAGUAGGGUUU"
#
#     two_seqs0, struct_dot_bracket0, mfe, s_score_sigmoid0, ML_model_score0, traceback_matrix0, rna_sequence_1_0, rna_sequence_2_0 = \
#         pred_structure_get_total_score(rna_sequence_1g, rna_sequence_2g)
#
#     # traceback_print_bps(traceback_matrix0, rna_sequence_1_0, rna_sequence_2_0)
#     print(f'Secondary Struct score: {s_score_sigmoid0} \t ML model score: {ML_model_score0}')
#     print("Dot-Bracket Notation:")
#     print(f'{two_seqs0}\n{struct_dot_bracket0}')
#     print(f"MFE: {mfe}")

