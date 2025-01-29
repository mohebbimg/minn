import sys
import os


# Add the project root directory to sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.MINN_main_algorithm import load_and_predict
from src.DP_duplex_structure_prediction import run_DP_alg_on_a_set


def main():
    # Check if any arguments are passed
    if len(sys.argv) < 2:
        print("Usage: python run_trained_MINN_on_testset_file.py <test_set.csv>")
        return

    test_set_fname = sys.argv[1]
    # for the input test file, first we need to predict secondary structures of samples with our DP algorithm:
    test_sec_struct_fname = test_set_fname.replace('.csv', '.dp_sec_struct.out')
    if not os.path.exists(test_sec_struct_fname):
        test_sec_struct_fname = run_DP_alg_on_a_set.pred_structures(test_set_fname)

    model_fullpath = os.path.normpath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "../model_of_MINN.keras"))

    load_and_predict(model_fullpath, test_sec_struct_fname)
    # delete temp files:
    os.remove(test_sec_struct_fname)


if __name__ == "__main__":
    main()
