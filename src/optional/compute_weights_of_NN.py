import csv
import multiprocessing
import os
from itertools import product

import pandas as pd
# import seaborn as sns
from datetime import datetime
import numpy as np
from keras.models import Sequential
from keras.src.metrics import Precision, Recall, AUC
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import precision_recall_curve, roc_curve, roc_auc_score, average_precision_score, auc
from keras.layers import Dense
import keras
from keras.models import load_model
from tensorflow.keras.backend import clear_session
import tensorflow as tf
import random
from collections import defaultdict
from sklearn.model_selection import train_test_split


def load_file_create_duplexes_labels(dataset_file):
    duplexes, labels, lines = [], [], []
    with open(dataset_file, newline='', encoding='utf-8-sig') as csvfile:
        tmp_reader = csv.reader(csvfile, delimiter=',')
        next(tmp_reader)  # Skip the header line
        for row in tmp_reader:
            tmp_microRNA, tmp_target, tmp_label = row[0].strip(), row[1].strip(), row[2].strip()
            duplexes.append([tmp_microRNA, tmp_target])
            labels.append(int(tmp_label))
            lines.append(row)

    return duplexes, labels, lines


def create_vectors(duplexes, labels):
    valid_pairs = sorted(['AU', 'UA', 'CG', 'GC', 'GU', 'UG'])
    nts_str = 'ACGU'
    max_miRNA_len, max_target_len = 25, 25
    dataset_X = []
    dataset_y = labels
    k = 0
    for sample_duplex in duplexes:
        miRNA = sample_duplex[0]  # in input file microRNA is 5p3p, no need to change
        target_candidate = sample_duplex[1][::-1]  # cts should be 3p5p

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
        mx_of_four_nts = np.zeros(
            (max_miRNA_len, 1, num_valid_pairs * num_valid_pairs * num_valid_pairs * num_valid_pairs))
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

        # Append to dataset
        dataset_X.append(np.concatenate((flattened_one_nt, flattened_two_nts, flattened_three_nts)))
        # dataset_X.append(np.concatenate((flattened_one_nt, flattened_two_nts)))
        # dataset_X.append(flattened_one_nt)

        k += 1

    dataset_X = np.array(dataset_X)
    dataset_y = np.array(dataset_y)

    print('Dataset shape:', dataset_X.shape, dataset_y.shape, f"Num. of positives: {sum(dataset_y)}")
    return dataset_X, dataset_y


def create_model(features_len):
    # Define the DNN model
    seed = 7
    np.random.seed(seed)
    model = Sequential()
    model.add(keras.Input(shape=(features_len,)))

    model.add(Dense(1, kernel_initializer='uniform', activation='sigmoid'))

    # Compile model with precision, recall, and AUPRC as metrics
    model.compile(optimizer='adam',
                  loss='binary_crossentropy',
                  metrics=[Precision(), Recall(), AUC(curve='PR')])  # AUC with PR curve

    return model


def visualize_weights(model, weights_fname):
    nts_str = 'ACGU'
    valid_pairs = sorted(['AU', 'UA', 'CG', 'GC', 'GU', 'UG'])

    # Extract the weights of the first dense layer
    weights = model.get_layer(index=0).get_weights()[0]  # Get weights only (not biases)
    data = []
    mx1_size = 25 * 25 * 6
    mx2_size = 25 * 5 * 6 * 6
    mx3_size = 25 * 3 * 6 * 6 * 6

    for i in range(len(weights)):
        if True or weights[i][0] > 0:  # all weights included.
            if i < mx1_size:
                # Convert number back to indices for the first matrix size
                mi_index = i // (25 * 6)
                cts_index = (i // 6) % 25
                bp_index = i % 6
                bp_name = valid_pairs[bp_index]

            elif (i - mx1_size) < mx2_size:
                # Convert number back to indices for the second matrix size
                num = i - mx1_size
                bp2_index = num % 6
                bp1_index = (num // 6) % 6
                cts_offset = (num // 36) % 5
                mi_index = num // (5 * 6 * 6)
                cts_index = max(0, mi_index - 2 + cts_offset)
                bp_name = valid_pairs[bp1_index] + '-' + valid_pairs[bp2_index]

            elif (i - mx1_size - mx2_size) < mx3_size:
                # Convert number back to indices for the third matrix size
                num = i - mx1_size - mx2_size
                bp3_index = num % 6
                bp2_index = (num // 6) % 6
                bp1_index = (num // 36) % 6
                cts_offset = (num // 216) % 3
                mi_index = num // (3 * 6 * 6 * 6)
                cts_index = max(0, mi_index - 1 + cts_offset)
                bp_name = valid_pairs[bp1_index] + '-' + valid_pairs[bp2_index] + '-' + valid_pairs[bp3_index]

            # Append indices and base-pair names to the data list
            data.append([mi_index, cts_index, bp_name, weights[i][0]])

    # Convert data to a DataFrame
    df = pd.DataFrame(data, columns=["mi_index", "cts_index", "bp_name", "weight"])

    # Sort the DataFrame by the weight column (last column)
    sorted_df = df.sort_values(by="weight", ascending=False)

    # Get the top 20 rows with the highest weights
    top_recs = sorted_df.head(50)

    pd.set_option('display.max_rows', None)  # Display all rows
    pd.set_option('display.max_columns', None)  # Display all columns

    print(top_recs)

    # Open the CSV file in write mode
    with open(weights_fname, 'w', newline='') as csvfile:
        # Create a CSV writer
        writer = csv.writer(csvfile)

        # Write the header row (column names)
        writer.writerow(["mi_index", "cts_index", "bp_name", "weight"])

        # Write the data rows
        writer.writerows(data)


def evaluate_model(labels, predictions):
    precision, recall, thresholds = precision_recall_curve(labels, predictions)

    plt.title(f'DNN Prediction, test size: {len(labels)}')
    plt.plot(recall, precision, 'b')
    aa = roc_curve(labels, predictions)
    plt.plot(aa[0], aa[1], 'g')

    print('Model performance:')
    print('AUC', roc_auc_score(labels, predictions))
    print('AP', average_precision_score(labels, predictions))
    print('AUCPR', auc(recall, precision))

    # plt.show()


def train_model_on_params(params):
    (ep, bs, X_train, y_train, val_X, val_Y, features_len) = params
    print(f"Training with epochs={ep}, batch_size={bs} ...")

    clear_session()
    model = create_model(features_len)

    # Train the model
    # class_weight = {0: 1.0, 1: 10.0}  # Adjust based on your imbalance ratio
    # class_weight = {0: 1.0, 1: (1 / 0.3496)}
    class_weight = {0: 1.0, 1: 1.0}
    try:
        model.fit(X_train, y_train, epochs=ep, batch_size=bs, validation_data=(val_X, val_Y),
                  verbose=0, class_weight=class_weight)
    except Exception as e:
        print(f"Error training model with epochs={ep}, batch_size={bs}: {e}")
        return -1, None, (ep, bs)  # Return -1 AUPRC if training fails

    # Evaluate on validation set
    val_loss, val_precision, val_recall, val_auprc = model.evaluate(val_X, val_Y, verbose=0)

    return val_auprc, model, (ep, bs)


def train_and_evaluate(X_train, y_train, val_X, val_y):
    # Convert features to NumPy arrays for DNN
    X_train = np.array(X_train)
    y_train = np.array(y_train)
    val_X = np.array(val_X)
    val_y = np.array(val_y)

    features_len = len(X_train[0])
    # epochs_arr = [100, 500, 1000, 2000, 5000, 10000]
    # batch_size_arr = [100, 50, 200, 25]

    epochs_arr = [5, 30, 40]
    batch_size_arr = [32, 64, 128, 256]

    best_auprc = -1
    best_model = None
    best_params = None

    for (ep, bs) in product(epochs_arr, batch_size_arr):
        val_auprc, model, params = train_model_on_params((ep, bs, X_train, y_train, val_X, val_y, features_len))
        if val_auprc > best_auprc:
            print(f"New best model found with AUPRC={val_auprc}")
            best_auprc = val_auprc
            best_model = model
            best_params = params

    print(f"Best parameters: Epochs={best_params[0]}, Batch size={best_params[1]} with AUPRC={best_auprc}")
    id = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    best_model.save(f'model_3nts_M_D_H_{id}.keras')
    return best_model


if __name__ == "__main__":
    training_set_file = "../../sample_data/training_set.csv"
    test_set_file = "../../sample_data/test_set.csv"
    validation_set_file = "../../sample_data/validation_set.csv"

    # Read data from file
    tmp_duplexes, tmp_labels, _ = load_file_create_duplexes_labels(test_set_file)
    test_all_X, test_all_y = create_vectors(tmp_duplexes, tmp_labels)

    tmp_duplexes, tmp_labels, _ = load_file_create_duplexes_labels(training_set_file)
    train_all_X, train_all_y = create_vectors(tmp_duplexes, tmp_labels)

    tmp_duplexes, tmp_labels, _ = load_file_create_duplexes_labels(validation_set_file)
    valid_X, valid_y = create_vectors(tmp_duplexes, tmp_labels)

    # Train model and save weights in file
    model = train_and_evaluate(train_all_X, train_all_y, valid_X, valid_y)
    visualize_weights(model, 'weights_for_DP_algorithm.csv')
    predictions = model.predict(test_all_X).flatten()
    evaluate_model(test_all_y, predictions)

