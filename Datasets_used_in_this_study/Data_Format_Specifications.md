# Data Format Specifications

## Program Input Files
The program requires three input files:
- **Training Set**
- **Validation Set**
- **Test Set**

All three files share the same structure and format.

## File Structure
Each file is a **comma-separated values (CSV)** file with the following columns:

| Column Name          | Description                                                                 |
|----------------------|-----------------------------------------------------------------------------|
| `microRNA_5p3p`     | The microRNA sequence in 5' to 3' orientation.                              |
| `target_site_5p3p`  | The target site sequence in 5' to 3' orientation.                           |
| `label`             | A numerical label indicating the interaction type:                         |
|                      | - `1`: Positive interaction (validated microRNA-target interaction)        |
|                      | - `0`: Negative interaction (no validated interaction)                    |
| `microRNA_family`   | The family to which the microRNA belongs, based on sequence similarity.     |
| `MFE`               | The Minimum Free Energy (MFE) of the microRNA and target site duplex.       |

## Example Data
Below is an example of the expected input format:

```
microRNA_5p3p,target_site_5p3p,label,microRNA_family,MFE
UGAGGUAGUAGGUUGUAUAGUU,ACUGCCUUGGGAUCCUAUGUUA,1,let-7,-20.5
AGCUUAUCACUGUUUGAUGAAC,UUACGUCGAAUGGUGAAGACAU,0,miR-21,-15.8
```

### Notes:
1. **Sequence Orientation**: Both `microRNA_5p3p` and `target_site_5p3p` must be provided in the 5' to 3' direction.
2. **Label Values**: Ensure labels are strictly binary (`0` or `1`).
3. **MFE Values**: Optional, our sample data files have these values, as we used the mfe values for selecting cts sequences, but these values are not used by our model for predicting labels.
4. **Delimiter**: Ensure the file uses commas as delimiters, with no extra spaces between fields.
