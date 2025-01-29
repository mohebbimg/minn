import sys
import os


# Add the project root directory to sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.run_trained_MINN_on_pairs_of_microRNA_and_cts import run_on_pairs_of_seqs


def main():
    microRNA_seq = "CAAAGUGCUUACAGUGCAGGUAG"
    cts_seq = "AUUCGCUACCGCCUGCAGCACUACU"
    pred = run_on_pairs_of_seqs(microRNA_seq, cts_seq)

    print(f"MicroRNA seq: {microRNA_seq}")
    print(f"Candidate Target Site: {cts_seq}")
    print(f"Prediction score: {pred}")


if __name__ == "__main__":
    main()
