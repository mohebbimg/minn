# Importing Keras Sequential Model
import csv
from collections import defaultdict
from datetime import datetime
from itertools import product
from keras.src.metrics import Precision, Recall, AUC
from tensorflow.keras.models import load_model

import keras
from keras.models import Sequential, Model
from tensorflow.python.keras.layers import concatenate
from tensorflow.python.keras.models import Input, save_model
from keras import callbacks
from tensorflow.keras.layers import Input, Conv2D, MaxPooling2D, Dropout, Flatten, Dense, concatenate

import numpy as np
from scipy.stats._mstats_basic import trim
from sklearn.preprocessing import StandardScaler
from tensorflow.keras.backend import clear_session

from sklearn.metrics import precision_recall_curve
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, average_precision_score, auc, roc_curve
from sklearn.model_selection import train_test_split

from src.DP_duplex_structure_prediction.predict_duplex_structure_by_DP import rna_secondary_structure_two_sequences

# probabilities of base pairing between nucleotides A, C, G, and U (4 by 4 possibilities).
base_pair_probs = [[0.0519250936329588, 0.0869812734082397, 0.15661423220973783, 0.4965393258426966],
                   [0.0869812734082397, 0.018858436007040482, 0.6978878551672114, 0.04730324365099321],
                   [0.15661423220973783, 0.6978878551672114, 0.030295899758809343, 0.12103653156068005],
                   [0.4965393258426966, 0.04730324365099321, 0.12103653156068005, 0.045476469372286915]]


def create_multi_input_model(microRNA_max_len, target_max_len, kernel_dim=3):
    kernel_size_var = (kernel_dim, kernel_dim)

    # Input layers for microRNA and target sequences
    input1 = Input(shape=(microRNA_max_len, target_max_len, 1), name='ch1')
    input2 = Input(shape=(microRNA_max_len, target_max_len, 1), name='ch2')
    input3 = Input(shape=(microRNA_max_len, target_max_len, 1), name='ch3')
    input4 = Input(shape=(microRNA_max_len, target_max_len, 1), name='ch4')

    # CNN branch for microRNA input
    x1 = Conv2D(32, kernel_size=kernel_size_var, activation='relu', padding='same')(input1)
    x1 = MaxPooling2D((2, 2), padding='same')(x1)
    x1 = Dropout(0.25)(x1)

    x1 = Conv2D(64, kernel_size=kernel_size_var, activation='relu', padding='same')(x1)
    x1 = MaxPooling2D((2, 2), padding='same')(x1)
    x1 = Dropout(0.25)(x1)

    x1 = Conv2D(128, kernel_size=kernel_size_var, activation='relu', padding='same')(x1)
    x1 = MaxPooling2D((2, 2), padding='same')(x1)
    x1 = Dropout(0.25)(x1)

    x1 = Flatten()(x1)  # Flatten for Dense layers

    # CNN branch for target input
    x2 = Conv2D(32, kernel_size=kernel_size_var, activation='relu', padding='same')(input2)
    x2 = MaxPooling2D((2, 2), padding='same')(x2)
    x2 = Dropout(0.25)(x2)

    x2 = Conv2D(64, kernel_size=kernel_size_var, activation='relu', padding='same')(x2)
    x2 = MaxPooling2D((2, 2), padding='same')(x2)
    x2 = Dropout(0.25)(x2)

    x2 = Conv2D(128, kernel_size=kernel_size_var, activation='relu', padding='same')(x2)
    x2 = MaxPooling2D((2, 2), padding='same')(x2)
    x2 = Dropout(0.25)(x2)

    x2 = Flatten()(x2)  # Flatten for Dense layers

    # CNN branch for target input
    x3 = Conv2D(32, kernel_size=kernel_size_var, activation='relu', padding='same')(input3)
    x3 = MaxPooling2D((2, 2), padding='same')(x3)
    x3 = Dropout(0.25)(x3)

    x3 = Conv2D(64, kernel_size=kernel_size_var, activation='relu', padding='same')(x3)
    x3 = MaxPooling2D((2, 2), padding='same')(x3)
    x3 = Dropout(0.25)(x3)

    x3 = Conv2D(128, kernel_size=kernel_size_var, activation='relu', padding='same')(x3)
    x3 = MaxPooling2D((2, 2), padding='same')(x3)
    x3 = Dropout(0.25)(x3)

    x3 = Flatten()(x3)  # Flatten for Dense layers

    # CNN branch for target input
    x4 = Conv2D(32, kernel_size=kernel_size_var, activation='relu', padding='same')(input3)
    x4 = MaxPooling2D((2, 2), padding='same')(x4)
    x4 = Dropout(0.25)(x4)

    x4 = Conv2D(64, kernel_size=kernel_size_var, activation='relu', padding='same')(x4)
    x4 = MaxPooling2D((2, 2), padding='same')(x4)
    x4 = Dropout(0.25)(x4)

    x4 = Conv2D(128, kernel_size=kernel_size_var, activation='relu', padding='same')(x4)
    x4 = MaxPooling2D((2, 2), padding='same')(x4)
    x4 = Dropout(0.25)(x4)

    x4 = Flatten()(x4)  # Flatten for Dense layers

    # Concatenate the outputs of both branches
    # concatenated = concatenate([x1, x2, x3, x4])
    concatenated = x1

    # Dense layers
    dense = Dense(128, activation='relu')(concatenated)
    dense = Dense(64, activation='relu')(dense)
    dense = Dropout(0.25)(dense)

    # Output layer
    output = Dense(1, activation='sigmoid')(dense)

    # Create the Model
    # multi_input_model = Model(inputs=[input1, input2, input3, input4], outputs=output)
    multi_input_model = Model(inputs=input1, outputs=output)

    # Compile the Model
    multi_input_model.compile(optimizer='adam',
                              loss='binary_crossentropy',
                              metrics=[Precision(), Recall(), AUC(curve='PR')])  # AUC with PR curve

    return multi_input_model


def load_and_create_matrices1(duplexes, labels, microRNA_max_len, target_max_len):
    dataset_X = []
    dataset_y = labels
    num_of_channels = 4
    k = 0
    for sample_duplex in duplexes:
        sample_mx = np.zeros((microRNA_max_len, target_max_len, num_of_channels))
        microRNA_5p3p = sample_duplex[0]
        micro_sec_struct = sample_duplex[2]
        cts_3p5p = sample_duplex[1][::-1]
        cts_sec_struct = sample_duplex[3][::-1]

        # get DP table from sec sturct pred module
        dp_table, _, mfe_mx = rna_secondary_structure_two_sequences(microRNA_5p3p, cts_3p5p)

        # print(miRNA, target_candidate)
        max_i = min(len(microRNA_5p3p), microRNA_max_len)
        max_j = min(len(cts_3p5p), target_max_len)
        # add dp table content to a channel
        for i in range(max_i):
            for j in range(max_j):
                sample_mx[i][j][2] = dp_table[i][j]
                sample_mx[i][j][3] = mfe_mx[i][j]

        # add bp probabilities to a channel
        for i in range(max_i):
            for j in range(max_j):
                try:
                    i_nt_index = 'ACGU'.index(microRNA_5p3p[i])
                    j_nt_index = 'ACGU'.index(cts_3p5p[j])
                except:
                    continue

                i_j_bp_prob = base_pair_probs[i_nt_index][j_nt_index]
                sample_mx[i][j][0] = i_j_bp_prob

        # add sec struct to a channel
        last_j = -1
        for i in range(max_i):
            if micro_sec_struct[i] == '(':
                for j in range(last_j + 1, max_j):
                    if cts_sec_struct[j] == ')':
                        last_j = j
                        break

                # sample_mx[i][last_j][1] = bp_score
                i_nt_index = 'ACGU'.index(microRNA_5p3p[i])
                try:
                    j_nt_index = 'ACGU'.index(cts_3p5p[last_j])
                except:
                    continue

                sample_mx[i][last_j][1] = base_pair_probs[i_nt_index][j_nt_index]

        dataset_X.append(np.array(sample_mx))
        k += 1

    dataset_X = np.array(dataset_X)
    dataset_y = np.array(dataset_y)

    # reshape for CNN:
    dataset_X = dataset_X.reshape(-1, microRNA_max_len, target_max_len, num_of_channels)
    dataset_X = dataset_X.astype('float32')

    # print('microRNA len:', microRNA_max_len, 'target len:', target_max_len)
    # print(dataset_X.shape, dataset_y.shape, 'Num. of positives:', sum(value > 0 for value in dataset_y))
    return dataset_X, dataset_y, num_of_channels


def evaluate_model(labels, predictions):
    test_type = 'hsa'
    precision, recall, thresholds = precision_recall_curve(labels, predictions)

    plt.title('Test size: ' + str(len(labels)))
    plt.plot(recall, precision, 'b')
    aa = roc_curve(labels, predictions)
    plt.plot(aa[0], aa[1], 'g')

    # write predictions in file
    with open('customizedRNAfold_CNN_' + test_type + '_prediction' + '.txt', 'w') as outfile:
        for i in range(0, len(predictions)):
            outfile.write(str(predictions[i][0]) + '\t' + str(labels[i]) + '\n')

    print('Model performance:')
    print('AUC', roc_auc_score(labels, predictions))
    print('AP', average_precision_score(labels, predictions))
    print('AUCPR', auc(recall, precision))

    plt.show()


def read_dp_outfile_and_extract_features(dataset_file):
    duplexes, labels = [[] for _ in range(1)], [[] for _ in range(1)]
    seed_len = 6
    link_len, link_seq = 1, "&"
    with open(dataset_file, newline='', encoding='utf-8-sig') as input_file:
        all_lines = input_file.readlines()
        i = 0
        while i < len(all_lines):
            if all_lines[i].startswith('>seq'):
                tmp_label = all_lines[i].split('label=')[1].strip()
                tmp_mi_link_target = all_lines[i + 1].strip()

                microRNA_5p3p = tmp_mi_link_target[:tmp_mi_link_target.index(link_seq)]
                cts_5p3p = tmp_mi_link_target[tmp_mi_link_target.index(link_seq) + link_len:]

                tmp_sec_struct = all_lines[i + 2].split()[0]
                mi_sec, cts_sec_struct = tmp_sec_struct[:tmp_mi_link_target.index(link_seq)], \
                    tmp_sec_struct[tmp_mi_link_target.index(link_seq) + link_len:]

                tmp_mfe = float(all_lines[i + 2].split()[-1].replace('(', '').replace(')', ''))
                i += 3

                seed = microRNA_5p3p[1:1 + seed_len]
                under_seed = cts_5p3p[::-1][1:1 + seed_len]
                # cnt = 0
                # for k in range(seed_len):
                #     if seed[k] + under_seed[k] in ['AU', 'UA', 'CG', 'GC', 'GU', 'UG']:
                #         cnt += 1

                duplexes[0].append([microRNA_5p3p, cts_5p3p, mi_sec, cts_sec_struct, tmp_mfe])
                labels[0].append(int(tmp_label))

    return duplexes, labels



def split_duplexes_by_miRNA(duplexes, labels, test_size_ratio=0.2):
    # Create a dictionary where the key is the microRNA and the value is a list of tuples (duplex, label)
    miRNA_dict = defaultdict(list)
    for duplex, label in zip(duplexes, labels):
        microRNA = duplex[0]  # The first sequence is the microRNA
        miRNA_dict[microRNA].append((duplex, label))

    # Split microRNAs into training and validation with a 4:1 ratio
    miRNAs = list(miRNA_dict.keys())
    train_miRNAs, val_miRNAs = train_test_split(miRNAs, test_size=test_size_ratio, random_state=42)

    # Initialize lists for training and validation duplexes and labels
    training_duplexes = []
    training_labels = []
    validation_duplexes = []
    validation_labels = []

    # Assign the duplexes and labels to the training and validation sets based on the split microRNAs
    for miRNA in train_miRNAs:
        for duplex, label in miRNA_dict[miRNA]:
            training_duplexes.append(duplex)
            training_labels.append(label)

    for miRNA in val_miRNAs:
        for duplex, label in miRNA_dict[miRNA]:
            validation_duplexes.append(duplex)
            validation_labels.append(label)

    return training_duplexes, training_labels, validation_duplexes, validation_labels



def train_and_test(training_sec_struct_file, validation_sec_struct_file, test_sec_struct_file):
    microRNA_max_len, target_max_len = 25, 25
    training_group_duplexes, training_group_labels = read_dp_outfile_and_extract_features(training_sec_struct_file)
    test_group_duplexes, test_group_labels = read_dp_outfile_and_extract_features(validation_sec_struct_file)
    validation_duplexes, validation_labels = read_dp_outfile_and_extract_features(test_sec_struct_file)

    # Read data from file
    print(f'Size of group: training = {len(training_group_duplexes[0])}, test = {len(test_group_duplexes[0])}')
    load_and_create_matrices = load_and_create_matrices1

    train_X, train_Y, num_of_channels = load_and_create_matrices(
        training_group_duplexes[0], training_group_labels[0], microRNA_max_len, target_max_len)
    test_X, test_Y, _ = \
        load_and_create_matrices(test_group_duplexes[0], test_group_labels[0], microRNA_max_len, target_max_len)
    # validation set:
    val_X, val_Y, _ = \
        load_and_create_matrices(validation_duplexes[0], validation_labels[0], microRNA_max_len, target_max_len)

    epochs_arr = [10]   # [20, 30, 40]
    batch_size_arr = [100]  # [128, 256, 512]

    best_auprc = -1
    predictions_file = f"{test_sec_struct_file.replace('.dp_sec_struct.out', '.predictions.txt')}"

    saved_model_file = f'model_of_MINN.tmp.keras'

    # separate channels:
    train_X_ch1 = train_X[:, :, :, 0:1]
    train_X_ch2 = train_X[:, :, :, 1:2]
    train_X_ch3 = train_X[:, :, :, 2:3]
    train_X_ch4 = train_X[:, :, :, 3:4]
    train_X_ch0 = train_X_ch1

    val_X_ch1 = val_X[:, :, :, 0:1]
    val_X_ch2 = val_X[:, :, :, 1:2]
    val_X_ch3 = val_X[:, :, :, 2:3]
    val_X_ch4 = val_X[:, :, :, 3:4]
    val_X_ch0 = val_X_ch1

    test_X_ch1 = test_X[:, :, :, 0:1]
    test_X_ch2 = test_X[:, :, :, 1:2]
    test_X_ch3 = test_X[:, :, :, 2:3]
    test_X_ch4 = test_X[:, :, :, 3:4]
    test_X_ch0 = test_X_ch1

    for (ep, bs) in product(epochs_arr, batch_size_arr):
        print(f"epochs={ep}, batch-size={bs} ...")
        clear_session()
        model = create_multi_input_model(microRNA_max_len, target_max_len)
        try:
            # model.fit([train_X_ch1, train_X_ch2, train_X_ch3, train_X_ch4], train_Y, epochs=ep, batch_size=bs,
            #           validation_data=([val_X_ch1, val_X_ch2, val_X_ch3, val_X_ch4], val_Y), verbose=0)

            model.fit(train_X_ch0, train_Y, epochs=ep, batch_size=bs,
                      validation_data=(val_X_ch0, val_Y), verbose=0)

        except Exception as e:
            print(f"Error training model with epochs={ep}, batch_size={bs}: {e}")
            continue

        # Evaluate on validation set
        # val_loss, val_precision, val_recall, val_auprc = model.evaluate([val_X_ch1, val_X_ch2, val_X_ch3, val_X_ch4], val_Y, verbose=0)
        val_loss, val_precision, val_recall, val_auprc = model.evaluate(val_X_ch0, val_Y, verbose=0)

        if val_auprc > best_auprc:
            print(f"New best model found with AUPRC={val_auprc}")
            best_auprc = val_auprc
            best_model = model
            best_params = (ep, bs)

            # test_loss, test_precision, test_recall, test_auprc = model.evaluate([test_X_ch1, test_X_ch2, test_X_ch3, test_X_ch4], test_Y, verbose=0)
            test_loss, test_precision, test_recall, test_auprc = model.evaluate([test_X_ch0], test_Y, verbose=0)


            print(f"Best results so far: AUPRC={test_auprc}, best params:{best_params}")
            best_model.save(saved_model_file)

            # Get predictions for the test set
            # predictions = best_model.predict([test_X_ch1, test_X_ch2, test_X_ch3, test_X_ch4])
            predictions = best_model.predict([test_X_ch0])

            # Save predictions and true labels to a file
            with open(predictions_file, "w") as f:
                f.write("Prediction\tLabel\n")
                for i in range(0, len(predictions)):
                    f.write(f"{predictions[i][0]}\t{test_Y[i]}\n")

    return predictions_file, saved_model_file


def load_and_predict(saved_model_file, test_sec_struct_file, one_pair=False, microRNA_max_len=25, target_max_len=25):
    # Load the saved model
    model = load_model(saved_model_file)
    print(f"Model loaded successfully from {saved_model_file}")

    # Load test data
    test_group_duplexes, test_group_labels = read_dp_outfile_and_extract_features(test_sec_struct_file)
    load_and_create_matrices = load_and_create_matrices1
    test_X, test_Y, _ = load_and_create_matrices(
        test_group_duplexes[0], test_group_labels[0], microRNA_max_len, target_max_len)

    # Separate channels for the test data
    test_X_ch1 = test_X[:, :, :, 0:1]
    test_X_ch2 = test_X[:, :, :, 1:2]
    test_X_ch3 = test_X[:, :, :, 2:3]
    test_X_ch4 = test_X[:, :, :, 3:4]

    # Make predictions on the test set
    predictions = model.predict([test_X_ch1, test_X_ch2, test_X_ch3, test_X_ch4])

    if one_pair is False:
        # check if input has labels
        true_labels = test_Y.flatten() if hasattr(test_Y, "flatten") else test_Y.tolist()
        if all(label == -1 for label in true_labels):
            no_label_provided = True
        else:
            no_label_provided = False

        # save predictions to a file
        predictions_file = f"{test_sec_struct_file.replace('.dp_sec_struct.out', '.predictions.txt')}"
        with open(predictions_file, "w") as f:
            if no_label_provided:
                f.write("Prediction\n")
            else:
                f.write("Prediction\tLabel\n")

            for i in range(len(predictions)):
                if no_label_provided:
                    f.write(f"{predictions[i][0]}\n")
                else:
                    f.write(f"{predictions[i][0]}\t{test_Y[i]}\n")

        print(f"Predictions made successfully on the input file.")
        print(f"Predictions saved to {predictions_file}")
        predicted_scores = [pred[0] for pred in predictions]

        if no_label_provided is False:
            auprc = average_precision_score(true_labels, predicted_scores)
            print(f"AUPRC on the input set: {auprc}")

    return predictions, test_Y

# if __name__ == '__main__':
#     main()
