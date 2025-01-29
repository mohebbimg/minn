# MicroRNA target-site prediction via a Multi-Input Neural Network

This project implements a multi-input neural network (MINN) for predicting microRNA's target sites. The code includes utilities for training, testing, and evaluating the model.

## Project Structure

```
project/
|-- scripts/
|   |-- train_and_test.py          # Main script to train and test the model
|   |-- run_trained_MINN_on_testset_file.py          # Runs the pre-trained model on an input test file.
|   |-- run_trained_MINN_on_pairs_of_microRNA_and_cts.py  # Runs the pre-trained model on a pair of microRNA and cts sequenes.
|-- src/
|   |-- DP_duplex_structure_prediction # Source codes of our DP algorithm for duplex structure prediction, and MFE computations
|   |-- MINN_main_algorithm.py    # Core model creation and training code
|-- model_of_MINN.keras           # Saved pre-trained model file
|-- requirements.yml              # Conda environment setup
```


## Installation Instructions

### Prerequisites
- Python 3.10
- Conda installed on your system

---

### Installation on **Windows**:

1. **Install Conda**:
   - Download and install Miniconda or Anaconda from [Conda's official website](https://docs.conda.io/en/latest/miniconda.html).

2. **Clone the Repository**:
   - Open the Command Prompt and run:
     ```bash
     git clone <repository-url>
     cd <repository-name>
     ```

3. **Create and Activate the Environment**:
   - Run the following commands in the same Command Prompt:
     ```bash
     conda env create -f requirements.yml
     conda activate multi_input_project
     ```

---

### Installation on **macOS**:

1. **Install Conda**:
   - Install Miniconda or Anaconda. Use `brew` to install Miniconda if you prefer:
     ```bash
     brew install --cask miniconda
     ```

2. **Clone the Repository**:
   - Open Terminal and run:
     ```bash
     git clone <repository-url>
     cd <repository-name>
     ```

3. **Create and Activate the Environment**:
   - Run the following commands in the Terminal:
     ```bash
     conda env create -f requirements.yml
     conda activate multi_input_project
     ```

---

### Installation on **Ubuntu/Linux**:

1. **Install Conda**:
   - Download Miniconda or Anaconda. Alternatively, install Miniconda via the command line:
     ```bash
     wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
     bash Miniconda3-latest-Linux-x86_64.sh
     ```

2. **Clone the Repository**:
   - Open Terminal and run:
     ```bash
     git clone <repository-url>
     cd <repository-name>
     ```

3. **Create and Activate the Environment**:
   - Run the following commands in the Terminal:
     ```bash
     conda env create -f requirements.yml
     conda activate multi_input_project
     ```

---

### Final Steps for All Platforms:

Once the environment is set up, verify the installation:
```bash
python --version
conda list
```

## Usage

### Training and Testing
Run the main script to train and evaluate the model:
```bash
python scripts/train_and_test.py training_set.csv validation_set.csv test_set.csv
```

### Using saved model to predict target sites
Run the following script to predict labels for samples in a file, here test_set.csv file:
```bash
python scripts/run_trained_MINN_on_testset_file.py <test_set.csv>
```
Run the following script to predict label for a pair of microRNA and candiate target site (cts) sequences:
```bash
python scripts/run_MINN_on_pairs_of_microRNA_and_cts.py <microRNA seq> <cts seq>
```


## Key Dependencies
The project relies on the following key libraries:
- `tensorflow`
- `keras`
- `numpy`
- `pandas`
- `scipy`
- `scikit-learn`
- `matplotlib`

For a complete list of dependencies, see (environment.yml).

## Documentation

### Scripts
- **train_and_test.py**: Runs training and evaluation. Outputs predictions and saves the best model.

### Source Code
- **src/DP_duplex_structure_prediction.py**: Functions for duplex structure modeling.
- **src/MINN_main_algorithm.py**: Implements the multi-input neural network (MINN) logic and utilities for dataset preprocessing.

### Output
- **Saved Model**: `model_of_MINN.keras`.
- **Predictions File**: Generated predictions for test data in `.txt` format.

## Contributions
Feel free to open issues or submit pull requests for improvements.

