# Data Processing Instructions

## mirTarBase Data Processing

### Step-by-Step Procedure
1. **Input File**:  
   - The file `miTarBase_microRNA_target_sequences.csv` contains the following columns:  
     - `mirtid`, `target_site_seq`, `gene_name_details`, `microRNA_seq`, `experiment_methods`, `microRNA_name`, `functional`, `target_source`, `seq_3p_UTR`.

2. **Extract Positive Samples**:  
   - Each record provides a **microRNA sequence**, a **target site sequence** (positive sample), and the corresponding **3'UTR sequence** containing the target site.

3. **Candidate Target Site (CTS) Generation**:  
   - On the **3'UTR sequence**, slide a window of **25 nucleotides** to generate all possible CTSs.

4. **Thermodynamic Filtering**:  
   - For each CTS, create a duplex with the corresponding microRNA sequence.  
   - Use the `DuplexFold` tool to compute the **Minimum Free Energy (MFE)** of the duplex.  
   - If the MFE is **less than 0**, select the pair of microRNA and CTS as a **negative sample**.

5. **Output Constraints**:  
   - Limit the selection to **up to 10 negative samples** per record in the file.

---

## Helwak Data Processing

### Step-by-Step Procedure
1. **Input File**:  
   - The file `Helwak_nonCanonical_targets.txt` contains records of:  
     - MicroRNAs,  
     - Binding regions of 3'UTRs,  
     - Interaction details.

2. **Extract Positive Samples**:  
   - From each record, locate the **microRNA seed binding site** on the 3'UTR sequence.

3. **Generate Subsequence**:  
   - Extract a **25-nucleotide subsequence** starting from the binding location on the 3'UTR sequence.

4. **Label Positive Pairs**:  
   - Pair each microRNA sequence with its corresponding 25-nucleotide target site sequence.  
   - Label these pairs as **positive samples**.

5. **Result**:  
   - This procedure generates **18,493 pairs** of microRNAs and target sites, all labeled as positive samples.

---

## DIANA-TarBase Data Processing

### Step-by-Step Procedure
1. **Input File**:  
   - The file `DIANA_tarBase_validated_transcripts.txt` is tab-separated with the following columns:  
     - `miRNA`, `gene_name`, `EnsemblId`, `Positive_Negative`, `Mature_mirna_transcript`, `3UTR_transcript`.

2. **Select Negative Cases**:  
   - Filter rows where `Positive_Negative == 0` to identify **validated negative interactions**.

3. **Candidate Target Site (CTS) Generation**:  
   - On the **3'UTR sequence**, slide a window of **25 nucleotides** to generate CTSs.

4. **Thermodynamic Filtering**:  
   - For each CTS, create a duplex with the corresponding **mature microRNA sequence**.  
   - Use the `DuplexFold` tool to compute the **MFE** of the duplex.  
   - If the MFE is **less than 0**, select the pair of microRNA and CTS as a **negative sample**.

5. **Output Constraints**:  
   - Limit the selection to **up to 10 negative samples** per record in the file.

---
