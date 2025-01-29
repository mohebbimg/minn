# Load the CSV data
import csv
from src.DP_duplex_structure_prediction.predict_duplex_structure_by_DP import pred_structure_get_total_score


def pred_structures(seq_file):
    sequences_1 = []
    sequences_2 = []
    labels = []

    with open(seq_file, newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        for row in reader:
            sequences_1.append(row[0])
            sequences_2.append(row[1])
            try:
                labels.append(int(row[2]))
            except:
                labels.append(-1)   # when the label not provided, fill it with -1 for unkown label

    sec_struct_fname = f'{seq_file.replace(".csv", "").replace(".CSV", "")}.dp_sec_struct.out'
    structure_scores = []
    with open(sec_struct_fname, 'w') as output_file:
        for idx, (seq1, seq2, label) in enumerate(zip(sequences_1, sequences_2, labels)):
            two_seqs, struct_dot_bracket, mfe_val, s_score_sigmoid, ML_model_score, traceback_matrix, rna_sequence_1, rna_sequence_2\
                = pred_structure_get_total_score(seq1, seq2)

            output_file.write(f">seq{idx+1} ML_score={ML_model_score:.3f} sec_score={s_score_sigmoid:.3f} label={label}\n")
            output_file.write(f"{two_seqs}\n")
            output_file.write(f"{struct_dot_bracket} ({mfe_val:.3f})\n")

            structure_scores.append(s_score_sigmoid)

    return sec_struct_fname

