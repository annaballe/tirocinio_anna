# REPORT ( analysis )

---

### Step 1: Parsing VCF Files and Converting to BED

**Input:**

- A compressed VCF file with genomic data on short tandem repeats (STRs).

**Objective:**

- Extract key genomic coordinates and information to produce a BED file for downstream analysis.

**Process:**

1. **Parsing the VCF File:**
    - The first script used was `parsing_test_set.py`.
    - **Usage:** `python parsing_test_set.py --file/--folder [input file or folder path] --output [output path]`
    - This script extracted the following fields from the VCF file: CHROM, POS, END, REP_UNIT, VAR_ID, and genotype information.
    - Allele counts (AC) were calculated from the genotype fields.
    - The parsed data was saved as a tab-separated `.txt` file.
2. **Converting .txt to BED:**
    - The second script used was `txt_to_bed.py`.
    - **Usage:** `python txt_to_bed.py variants.txt output.bed`
    - This script processed the `.txt` file to extract CHROM, zero-based POS (start), END, and VAR_ID fields.
    - Ensured compliance with the BED format, producing a tab-delimited, headerless output.

**Output:**

- Parsed data in `.txt` format.
- Final BED file containing genomic coordinates and variant IDs.

**Additional Commands Used:**

- To extract information directly from the VCF and create a BED file, the following command was executed in the terminal:
    
    ```bash
    bcftools query -f '%CHROM\\t%POS0\\t%END\\t%ID\\n' HG002_GRCh38_TandemRepeats_v1.0.1.vcf.gz > true_set.bed
    
    ```
    

### Step 2: Building the Confusion Matrix and Computing Metrics

**Input:**

- BED files: `call_set_ordered.bed` and `true_set.bed`.

**Objective:**

- Compute overlaps between the call set and the true set to build a confusion matrix for evaluating performance metrics.

**Process:**

1. **Computing Intersections:**
    - Two intersection operations were performed using `intersectBed`:
        - To identify true positives and other overlaps:
            
            ```bash
            intersectBed -a call_set_ordered.bed -b true_set.bed -f 0.80 -r -wao > intersect.bed
            
            ```
            
        - To identify false negatives:
            
            ```bash
            intersectBed -a true_set.bed -b call_set_ordered.bed -f 0.80 -r -wao > FN.bed
            
            ```
            

**Output:**

- `intersect.bed`: Contains overlaps between the call set and the true set.
- `FN.bed`: Contains regions in the true set that were not matched in the call set.

### Step 3: Classifying Variants

**Objective:**

- Annotate and classify variants from the intersections into categories: True Positives (TP), False Positives (FP), and False Negatives (FN).

**Process:**

1. **Classifying Variants in `intersect.bed`:**
    - A script, `classify_variants.py`, was used for this task:
        
        ```python
        import pandas as pd
        import sys
        
        # Read input file (use command-line argument to specify the file)
        input_file = sys.argv[1]
        
        # Load data into a DataFrame
        df = pd.read_csv(input_file, sep='\\t', header=None)  # Assuming tab-separated format
        
        # Add a new column 'state' based on conditions in columns 5, 6, and 8 (adjust as necessary)
        df['state'] = df.apply(
            lambda row: 'fp' if row[8] == "0" else  # If the value in column 8 (OVERLAPP_50) is "0", mark as 'fp'
                       'fn' if row[8] == "." else  # If the value in column 8 (OVERLAPP_50) is ".", mark as 'fn'
                       'tp',  # For other cases, mark as 'tp' (adjust this as needed for your specific conditions)
            axis=1
        )
        
        # Save the output to a new file
        output_file = input_file.replace('.bed', '_classified.bed')
        df.to_csv(output_file, sep='\\t', index=False, header=False)
        
        print(f"Classification complete. Results saved to {output_file}")
        
        ```
        
2. **Classifying Variants in `FN.bed`:**
    - The following script was used:
        
        ```python
        import pandas as pd
        import sys
        
        # Input file (FN.bed) provided as a command-line argument
        fn_bed_file = sys.argv[1]
        
        # Read FN.bed (tab-separated with no header)
        fn_df = pd.read_csv(fn_bed_file, sep="\\t", header=None)
        
        # Verify the number of columns
        print(f"Number of columns in {fn_bed_file}: {len(fn_df.columns)}")
        
        # Assign column names (adjust to 9 columns)
        fn_df.columns = ['chr', 'start', 'end', 'name', 'score', 'strand', 'value1', 'value2', 'overlap']
        
        # Label rows based on the overlap column
        fn_df['label'] = fn_df['overlap'].apply(lambda x: 'fn' if x == 0 else 'tp')
        
        # Save the labeled DataFrame to a file
        output_file = "labeled_FN.bed"
        fn_df.to_csv(output_file, sep="\\t", index=False)
        print(f"Labeled FN.bed saved to {output_file}")
        
        ```
        

**Output:**

- `intersect_classified.bed`: Annotated with FP, FN, and TP labels.
- `labeled_FN.bed`: Contains labeled rows for false negatives.

### Step 4: Computing the Confusion Matrix

**Objective:**

- Aggregate the results from classified intersection files to build the confusion matrix, summarizing TP, FP, FN, and TN counts.

**Process:**

- The following script was used:
    
    ```python
    import pandas as pd
    import sys
    
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python confusion_matrix.py <intersect_classified.bed> <labeled_FN.bed>")
        sys.exit(1)
    
    # Read the files without header
    intersect_df = pd.read_csv(sys.argv[1], sep="\\t", header=None)
    fn_df = pd.read_csv(sys.argv[2], sep="\\t", header=None)
    
    # Print the number of columns for each DataFrame to debug
    print(f"Number of columns in intersect_df: {intersect_df.shape[1]}")
    print(f"Number of columns in fn_df: {fn_df.shape[1]}")
    
    # Assign column names based on the actual number of columns
    # You can modify the column names according to the structure of your .bed file
    if intersect_df.shape[1] == 10:
        intersect_df.columns = ['chr', 'start', 'end', 'name', 'score', 'strand', 'value1', 'value2', 'value3', 'label']
    elif intersect_df.shape[1] == 9:
        intersect_df.columns = ['chr', 'start', 'end', 'name', 'score', 'strand', 'value1', 'value2', 'label']
    else:
        print("Unexpected number of columns in intersect_df.")
        sys.exit(1)
    
    if fn_df.shape[1] == 10:
        fn_df.columns = ['chr', 'start', 'end', 'name', 'score', 'strand', 'value1', 'value2', 'value3', 'label']
    elif fn_df.shape[1] == 9:
        fn_df.columns = ['chr', 'start', 'end', 'name', 'score', 'strand', 'value1', 'value2', 'label']
    else:
        print("Unexpected number of columns in fn_df.")
        sys.exit(1)
    
    # Count the labels for TP, FP, and FN
    tp = len(intersect_df[intersect_df['label'] == 'tp'])
    fp = len(intersect_df[intersect_df['label'] == 'fp'])
    fn = len(fn_df[fn_df['label'] == 'fn'])
    
    # Assume true negatives are based on external logic (example placeholder value)
    # Adjust this based on genome size or experiment-specific details
    tn = 1000000
    
    ```
    

### Step 5: Calculating Evaluation Metrics

**Objective:**

- Compute evaluation metrics (e.g., Precision, Recall, F1 Score) from the confusion matrix.

**Process:**

- A script, `calculate_metrics.py`, was used to compute metrics from the confusion matrix.

**Script for Metric Calculation:**

```python
import sys
import pandas as pd

# Check if the file path is provided as a command-line argument
if len(sys.argv) != 2:
    print("Usage: python

```