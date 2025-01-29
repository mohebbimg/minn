import sys
import os

# Add the project root directory to sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.DP_duplex_structure_prediction import run_DP_alg_on_a_set
from src.MINN_main_algorithm import train_and_test, load_and_predict


def main():
    # Check if any arguments are passed
    if len(sys.argv) < 3:
        print("Usage: python train_test.py <training_set.csv> <validation_set.csv> <test_set.csv>")
        return

    training_set_fname = sys.argv[1]
    validation_set_fname = sys.argv[2]
    test_set_fname = sys.argv[3]

    training_sec_struct_fname = run_DP_alg_on_a_set.pred_structures(training_set_fname)
    validation_sec_struct_fname = run_DP_alg_on_a_set.pred_structures(validation_set_fname)
    test_sec_struct_fname = run_DP_alg_on_a_set.pred_structures(test_set_fname)

    predictions_file, saved_model_file = \
        train_and_test(training_sec_struct_fname, validation_sec_struct_fname, test_sec_struct_fname)

    print(f"Prediction results were saved in {predictions_file}")
    print(f"Trained model was saved in {saved_model_file}")

    # model_fullpath = os.path.normpath(
    #     os.path.join(os.path.dirname(os.path.abspath(__file__)), "../model_of_MINN.keras"))
    #
    # load_and_predict(model_fullpath, test_sec_struct_fname)

if __name__ == "__main__":
    main()
