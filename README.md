# microRNA Prediction

This project uses deep learning methods to predict microRNA target sites. It includes model training, evaluation, and data preprocessing steps.

## Installation

### Using `requirements.txt`:

1. Clone the repository:
    ```bash
    git clone https://github.com/your-username/microRNA-prediction.git
    cd microRNA-prediction
    ```

2. Create a Python virtual environment (optional but recommended):
    ```bash
    python3 -m venv venv
    source venv/bin/activate  # On Windows use `venv\Scripts\activate`
    ```

3. Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

### Using Conda (`environment.yml`):

1. Clone the repository:
    ```bash
    git clone https://github.com/your-username/microRNA-prediction.git
    cd microRNA-prediction
    ```

2. Create the conda environment:
    ```bash
    conda env create -f environment.yml
    conda activate microrna-prediction
    ```

## Usage

To preprocess the data:

```bash
python scripts/preprocess_data.py --input raw_data.csv --output processed_data.csv

