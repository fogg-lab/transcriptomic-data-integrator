from typing import Iterable, List
import subprocess
import shutil
import os.path
import re

import pkg_resources
import requests
import pandas as pd
import mygene


R_SCRIPTS_DIR = pkg_resources.resource_filename('transcriptomics_data_query', 'rscripts')
MICROARRAY_NORMALIZATION_SCRIPT = os.path.join(R_SCRIPTS_DIR, 'rma_normalization.R')
RNASEQ_NORMALIZATION_SCRIPT = os.path.join(R_SCRIPTS_DIR, 'rnaseq_normalization.R')
BATCH_CORRECTION_SCRIPT = os.path.join(R_SCRIPTS_DIR, 'batch_correction.R')


def normalize_microarray(input_dir, output_file, remove_cel_dir=False):
    """Normalize microarray expression data given a directory containing CEL.gz files.

    Args:
        input_dir (str): Path to the directory containing CEL.gz files.
        output_file (str): Path to the output file.
        remove_cel_dir (bool, optional): If True, remove the input directory after normalization.
            Defaults to False.
    """
    command = ["Rscript", MICROARRAY_NORMALIZATION_SCRIPT, input_dir, output_file]

    print(f"Executing: {' '.join(command)}")
    subprocess.run(command, check=True)
    print(f"Normalization complete. Output file: {output_file}")

    if remove_cel_dir:
        shutil.rmtree(input_dir)


def normalize_rnaseq(expression_file, clinical_file, output_file):
    """Normalize RNA-seq expression data given a file containing raw counts.

    Args:
        expression_file (str): Path to the input file containing raw counts.
        clinical_file (str): Path to the input file containing clinical data.
        output_file (str): Path to the output file.
    """
    command = ["Rscript", RNASEQ_NORMALIZATION_SCRIPT, expression_file, clinical_file, output_file]

    print(f"Executing: {' '.join(command)}")
    subprocess.run(command, check=True)
    print(f"Normalization complete. Output file: {output_file}")


def normalize(input_path, output_file, clinical_file=None):
    """Normalize microarray or RNASeq expression data.

    Args:
        input_path (str): Path to the input file (for RNASeq) or directory (for microarray).
        output_file (str): Path to the output file.
        clinical_file (str, optional): Path to the input file containing clinical data.
            Not required for microarray.
    """
    if os.path.isdir(input_path):
        normalize_microarray(input_path, output_file, remove_cel_dir=True)
    else:
        normalize_rnaseq(input_path, clinical_file, output_file)


def load_genes_from_file(filename):
    """Read genes from a text file with one gene symbol per line.

    Args:
        filename (str): Path to the text file containing gene symbols.

    Returns:
        list: List of gene symbols.
    """
    with open(filename, 'r', encoding='utf-8') as file:
        genes = [line.strip() for line in file]

    return genes


def get_genes_from_msig_set(gene_set_name, species="human"):
    """Fetches the genes associated with a given gene set name from the MSigDB (e.g, NABA_MATRISOME).

    This function constructs a URL for the specified gene set name and species, then performs a GET request to fetch
    the associated genes in JSON format from the Molecular Signatures Database (MSigDB).

    Args:
        gene_set_name (str): The name of the gene set for which to fetch the associated genes.
        species (str, optional): The species for which to fetch the gene set. Defaults to "human".

    Returns:
        list[str]: A list of gene symbols associated with the specified gene set name.

    Raises:
        HTTPError: If the GET request to the MSigDB results in an error.
    """
    url = f"https://www.gsea-msigdb.org/gsea/msigdb/{species}/download_geneset.jsp?geneSetName={gene_set_name}&fileType=json"
    response = requests.get(url, timeout=10)
    response.raise_for_status()
    return response.json()[gene_set_name]["geneSymbols"]


def convert_genes(genes: Iterable, in_format: str, out_format: str, species: str="human",
                  returnall: bool=False) -> pd.Series:
    """Converts a list of genes between formats 'entrezgene', 'ensembl.gene', and 'symbol'.

    Args:
        genes (Union[List, pd.Series]): A list of genes.
        in_format (str): The format of the input genes.
        out_format (str): The format of the output genes.
        species (str, optional): The species of the genes. Defaults to "human".
        returnall (bool, optional): Whether to return return complete lists of duplicate
            or missing query terms. Defaults to False.

    Returns:
        pd.Series: Query results. Index is the input genes, values are the output genes.
    """
    # Validate in_format and out_format
    valid_formats = ["entrezgene", "ensembl.gene", "symbol"]
    given_invalid = [f for f in [in_format, out_format] if f not in valid_formats]
    if len(given_invalid) > 0:
        raise ValueError(f"Invalid format(s) given: {given_invalid}. Valid formats: {valid_formats}")
    elif in_format == out_format:
        print("Input and output formats are the same. Returning input genes.")
        return pd.Series(genes, index=genes)

    mg = mygene.MyGeneInfo()
    out = mg.querymany(genes, scopes=in_format, fields=out_format, species=species,
                       as_dataframe=True, returnall=returnall)

    return out[out_format]


def select_rows(df, values, column=None):
    """Select rows in DataFrame.

    Args:
        df (pandas.DataFrame): The data frame.
        values (list): The values to select (e.g., ["A2M","A2ML1","ABI3BP"]).
        column (str, optional): The column name (e.g., "gene_symbol"). If None, the index is used.

    Returns:
        pandas.DataFrame: The selected rows.

    Example:
        >>> import transcriptomics_data_query as tdq
        >>> expression_df = pd.DataFrame({"GSM1234": [3.452, 4.123, 5.678, 6.789],
                                          "GSM5678": [1.234, 2.345, 3.456, 4.567]})
        >>> expression_df.index = ["A1BG", "A2M", "CA10", "SEMA6B"]
        >>> expression_df.index.name = "symbol"
        >>> matrisome_genes = tdq.preprocess.get_genes_from_msig_set("NABA_MATRISOME")
        >>> matrisome_expression_df = tdq.preprocess.select_rows(expression_df, matrisome_genes)
        >>> matrisome_expression_df
                GSM1234 GSM5678
        symbol
           A2M    4.123   2.345
        SEMA6B    6.789   4.567
    """
    if column is None:
        return df.loc[df.index.isin(values)]

    return df[df[column].isin(values)]


def drop_nan_row_indices(expr_df: pd.DataFrame):
    """Drop rows where the row index is NaN.

    Args:
        expr_df (pandas.DataFrame): The expression matrix.

    Returns:
        pandas.DataFrame: The expression matrix with NaN rows dropped.
    """
    return expr_df.loc[expr_df.index.dropna()]


def clean_clinical_data(clinical_df: pd.DataFrame, specification: dict, ignore_case: bool=True,
                        drop_no_match_samples: bool=True) -> (pd.DataFrame, pd.DataFrame):
    """Get filtered and cleaned clinical data table based on a filter specification.
        The specification contains column names the patterns to extract values for each column.
        When multiple patterns match, (e.g. "disease" and "no disease"), the longer match is chosen.

    Args:
        clinical_df (pd.DataFrame): The clinical data table.
        specification (dict): The filter specification.
            Example: `specification={'condition': ['tumor', 'normal'], 'age': [r'\d+']}`
        ignore_case (bool, optional): Whether to ignore case when matching patterns. Defaults to True.
        drop_no_match_samples (bool, optional): Whether to drop samples that do not match any patterns
            for any column. Defaults to True.

    Returns:
        pd.DataFrame: The filtered and cleaned clinical data table.

    Raises:
        ValueError: If a column name in the specification is not in the clinical data table.
        ValueError: If a pattern in the specification is not a valid regular expression.
        ValueError: If no patterns match one of the values in a column and 

    Example:
        >>> import transcriptomics_data_query as tdq
        >>> import pandas as pd
        >>> clinical_df = pd.DataFrame({
                "condition": ["cns tumor tissue", "normal tissue", "normal tissue", "tumor tissue"],
                "age": ["45", "50 yo", "eta: 55 anni", "75 years old"],
                "organism": ["homosapiens", "homosapiens", "homosapiens", "homosapiens"]},
                index=["sample1", "sample2", "sample3", "sample4"])
        >>> clinical_df.index.name = "sample_name"
        >>> specification = {'condition': ['tumor', 'normal'], 'patient_age': [r'\d+']}
        >>> cleaned_clinical_df = tdq.preprocess.clean_clinical_data(clinical_df, specification)
        >>> print(cleaned_clinical_df)
                        condition patient_age
            sample_name                      
            sample1         tumor          45
            sample2        normal          50
            sample3        normal          55
            sample4         tumor          75
    """

    # Check for non-existing columns
    for col in specification:
        if col not in clinical_df.columns:
            raise ValueError(f"Column {col} not in clinical data table")

    # Function to find the longest first match among patterns
    def match_longest_pattern(value, patterns):
        longest_match = ""
        for pattern in patterns:
            try:
                match = re.search(pattern, value, flags=re.IGNORECASE if ignore_case else 0)
                if match:
                    match_str = match.group(0)
                    if len(match_str) > len(longest_match):
                        longest_match = match_str
            except re.error:
                raise ValueError(f"Pattern {pattern} is not a valid regular expression")
        if not longest_match:
            if drop_no_match_samples:
                return None
            raise ValueError(f"No patterns match one of the values in {value}")

        return longest_match

    cleaned_df = pd.DataFrame(index=clinical_df.index)

    # Apply filter specification to each column
    for col, patterns in specification.items():
        cleaned_df[col] = clinical_df[col].apply(match_longest_pattern, patterns=patterns)

    # Drop rows with NaN values
    print("Dropping samples with column values that do not match any patterns:")
    print(cleaned_df.index[cleaned_df.isna().any(axis=1)])
    cleaned_df = cleaned_df.dropna()

    return cleaned_df


def join_expression_matrices(expression_dataframes: List[pd.DataFrame]):
    """Concatenate two expression matrices with the same row names.

    Args:
        expression_dataframes (list[pd.DataFrame]): A list of expression matrices.

    Returns:
        pd.DataFrame: The concatenated expression matrix.

    Raises:
        ValueError: If the expression matrices do not all have the same index (row names).
    """

    # Check that all expression matrices have the same index
    indices = [df.index for df in expression_dataframes]
    if not all([i.equals(indices[0]) for i in indices]):
        raise ValueError("Expression matrices do not all have the same index.")

    return pd.concat(expression_dataframes, axis=1)


def join_and_batch(expression_dataframes: List[pd.DataFrame],
                   clinical_dataframes: List[pd.DataFrame]) -> (pd.DataFrame, pd.DataFrame):
    """Join expression matrices, join the corresponding clinical data tables, and assign batches.

    This function adds a batch column to the joined clinical table and assigns batch numbers.
    For instance, if we have `batches([expr_df1, expr_df2], [clin_df1, clin_df2])`,
    samples in `expr_df1`/`clin_df1` are batch 1 and samples in `expr_df2`/`clin_df2` are batch 2.

    Args:
        expression_dataframes (list[pd.DataFrame]): A list of expression matrices.
        clinical_dataframes (list[pd.DataFrame]): A list of clinical data tables.

    Returns:
        pd.DataFrame: The joined expression matrix.
        pd.DataFrame: The joined clinical data table.

    Raises:
        ValueError: If the column names in an expression matrix do not match the index (row names)
            in the clinical data table.
        ValueError: If the clinical data tables do not all have the same column names.

    Example:
        >>> import transcriptomics_data_query as tdq
        >>> import pandas as pd
        >>> # Dummy expression dataframes
        >>> expr_df1 = pd.DataFrame({"sample1": [1, 2], "sample2": [3, 4]},
                                    index=["gene1", "gene2"])
        >>> expr_df2 = pd.DataFrame({"sample3": [5, 6], "sample4": [7, 8]},
                                    index=["gene1", "gene2"])
        >>> expr_df3 = pd.DataFrame({"sample5": [9, 10], "sample6": [11, 12]},
                                    index=["gene1", "gene2"])
        >>> # Dummy clinical dataframes
        >>> clinical_df1 = pd.DataFrame({"condition": ["tumor", "normal"]},
                                        index=["sample1", "sample2"])
        >>> clinical_df2 = pd.DataFrame({"condition": ["normal", "tumor"]},
                                        index=["sample3", "sample4"])
        >>> clinical_df3 = pd.DataFrame({"condition": ["tumor", "tumor"]},
                                        index=["sample5", "sample6"])
        >>> # Lists of expression and clinical dataframes
        >>> expression_dataframes = [expr_df1, expr_df2, expr_df3]
        >>> clinical_dataframes = [clinical_df1, clinical_df2, clinical_df3]
        >>> # Using the batches function
        >>> expr_joined, clinical_joined = tdq.preprocess.batches(expression_dataframes,
                                                                  clinical_dataframes)
        >>> print("Joined Expression Matrix:")
        >>> print(expr_joined)
        >>> print("\nJoined Clinical Data Table:")
        >>> print(clinical_joined)
            Joined Expression Matrix:
                   sample1  sample2  sample3  sample4  sample5  sample6
            gene1        1        3        5        7        9       11
            gene2        2        4        6        8       10       12

            Joined Clinical Data Table:
                    condition  batch
            sample1     tumor      1
            sample2    normal      1
            sample3    normal      2
            sample4     tumor      2
            sample5     tumor      3
            sample6     tumor      3
    """
    # Check if clinical dataframes have same columns
    clinical_columns = clinical_dataframes[0].columns
    if not all([df.columns.equals(clinical_columns) for df in clinical_dataframes]):
        raise ValueError("Clinical data tables do not all have the same column names.")

    # Validate that expression matrices and clinical data table indices match
    for i, (expr_df, clin_df) in enumerate(zip(expression_dataframes, clinical_dataframes)):
        if not expr_df.columns.equals(clin_df.index):
            colname = expr_df.columns[i]
            rowname = clin_df.index[i]
            raise ValueError(f"Column {i} in expression matrix ({colname}) does not match "
                             f"row {i} in the clinical data table ({rowname}).")

    # Concatenate clinical dataframes and assign batch numbers
    for i, df in enumerate(clinical_dataframes):
        df['batch'] = i + 1
    clinical_df_joined = pd.concat(clinical_dataframes, axis=0)

    # Concatenate expression matrices
    expression_df_joined = join_expression_matrices(expression_dataframes)

    return expression_df_joined, clinical_df_joined


def batch_correction(expression_file: str, clinical_file: str, variable: str, data_type: str,
                     factor_levels: List[str], output_file: str):
    """Perform batch correction on expression data.

    Args:
        expression_file (str): Path to the input file containing expression data.
        clinical_file (str): Path to the input file containing clinical data.
        variable (str): The column in the clinical data table to use for batch correction.
        data_type (str): The data type of the expression data (microarray or RNASeq).
        factor_levels (list[str]): The factor levels of the variable to use for batch correction.
            This is used to validate that the clinical data table has the correct factor levels.
        output_file (str): Path to the output file.

    Raises:
        ValueError: If the clinical data table does not have the correct factor levels.
        ValueError: If the data type is not "microarray" or "RNASeq".
        ValueError: If the variable is not in the clinical data table.
        ValueError: If the column names in the expression data do not match the row names
            in the clinical data table.
    """

    if data_type.lower() not in ["microarray", "rnaseq"]:
        raise ValueError("Invalid data type. Must be either 'microarray' or 'RNASeq'.")

    clinical_data = pd.read_csv(clinical_file, sep='\t', index_col=0)
    if variable not in clinical_data.columns:
        raise ValueError(f"Specified condition column '{variable}' not found in clinical data table.")

    if set(clinical_data[variable]) != set(factor_levels):
        raise ValueError(f"The factor levels in the clinical data table do not match the specified factor levels: {factor_levels}.")

    expression_data = pd.read_csv(expression_file, sep='\t', index_col=0)
    if set(expression_data.columns) != set(clinical_data.index):
        raise ValueError("The column names in the expression data do not match the row names in the clinical data table.")

    command = ["Rscript", BATCH_CORRECTION_SCRIPT, data_type, expression_file, clinical_file,
               variable, output_file]

    print(f"Executing: {' '.join(command)}")
    subprocess.run(command, check=True)
    print(f"Batch correction complete. Output file: {output_file}")
