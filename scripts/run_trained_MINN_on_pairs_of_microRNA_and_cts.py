import sys
import os

# Add the project root directory to sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.DP_duplex_structure_prediction import run_DP_alg_on_a_set
from src.MINN_main_algorithm import train_and_test, load_and_predict


def run_on_pairs_of_seqs(microRNA_seq, cts_seq):
    tmp_fname = "temp_rundata.csv"
    with open(tmp_fname, 'w') as output_file:
        output_file.write(f"microRNA_5p3p,target_site_5p3p,label\n")
        output_file.write(f"{microRNA_seq},{cts_seq},{-1}\n")

    # predict secondary structure of input pair, with our DP algorithm:
    tmp_sec_struct_fname = run_DP_alg_on_a_set.pred_structures(tmp_fname)

    model_fullpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../model_of_MINN.keras")

    predictions, _ = load_and_predict(model_fullpath, tmp_sec_struct_fname, True)

    # delete temp files:
    os.remove(tmp_fname)
    os.remove(tmp_sec_struct_fname)


    return predictions[0][0]


def main():
    if len(sys.argv) < 2:
        print("Usage: python run_MINN_on_pairs_of_microRNA_and_cts.py <microRNA seq> <Candidate Target Site (CTS) seq>")
        return

    microRNA_seq = sys.argv[1]
    cts_seq = sys.argv[2]
    pred = run_on_pairs_of_seqs(microRNA_seq, cts_seq)
    print(f"MicroRNA seq: {microRNA_seq}")
    print(f"Candidate Target Site: {cts_seq}")
    print(f"Prediction score: {pred}")


if __name__ == "__main__":
    main()
