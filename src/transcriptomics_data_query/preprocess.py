import subprocess
import shutil
import os.path

import pkg_resources
import requests
import pandas as pd


CORE_MATRISOME_URL = "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=NABA_CORE_MATRISOME&fileType=json"
ALL_MATRISOME_URL = "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=NABA_MATRISOME&fileType=json"


def normalize_microarray(input_dir, output_file, remove_cel_dir=False):
    """Normalize microarray expression data given a directory containing CEL.gz files.

    Args:
        input_dir (str): Path to the directory containing CEL.gz files.
        output_file (str): Path to the output file.
        remove_cel_dir (bool, optional): If True, remove the input directory after normalization. Defaults to False.
    """
    r_script_path = pkg_resources.resource_filename('transcriptomics_data_query',
                                                    'rscripts/rma_normalization.R')
    print(f"Executing: Rscript {r_script_path} {input_dir} {output_file}")
    subprocess.run(["Rscript", r_script_path, input_dir, output_file], check=True)
    print("Normalization complete.")
    if remove_cel_dir:
        shutil.rmtree(input_dir)


def normalize_rnaseq(expression_file, clinical_file, output_file, norm_method="mrn",
                     design_col=None):
    """Normalize RNA-seq expression data given a file containing raw counts.

    Args:
        expression_file (str): Path to the input file containing raw counts.
        clinical_file (str): Path to the input file containing clinical data.
        output_file (str): Path to the output file.
        norm_method (str, optional): Normalization method (tmm or mrn).
                                     Defaults to "mrn".
    """
    r_script_path = pkg_resources.resource_filename('transcriptomics_data_query',
                                                    'rscripts/rnaseq_normalization.R')

    # Validate normalization method
    if not isinstance(norm_method, str):
        raise TypeError("Normalization method must be a string.")

    norm_method = norm_method.lower()

    if norm_method not in ("tmm", "mrn"):
        raise ValueError("Invalid normalization method. Must be 'tmm' or 'mrn'.")

    if norm_method == "tmm" and design_col is None:
        # Read clinical data
        clinical_data = pd.read_csv(clinical_file, sep='\t')

        # Check the clinical data structure
        if len(clinical_data.columns) <= 1:
            raise ValueError("Clinical data must contain more than 1 column.")

        if clinical_data.columns[0] != "sample_id":
            raise ValueError("The first column of the clinical data must be named 'sample_id'.")

        # List column options, excluding "sample_id"
        design_columns = [col for col in clinical_data.columns if col != "sample_id"]

        print("Available columns for design formula:")
        for idx, col_name in enumerate(design_columns):
            print(f"{idx}: {col_name}")

        selected_idx = int(input("Enter the index of the column to use in the design formula: "))
        design_col = design_columns[selected_idx]

    design_col = str(design_col)    # In case it is a NoneType

    print(f"Executing: Rscript {r_script_path} {norm_method} {expression_file} {clinical_file} {design_col} {output_file}")
    subprocess.run(["Rscript", r_script_path, norm_method, expression_file, clinical_file,
                    design_col, output_file], check=True)
    print("Normalization complete.")


def normalize(input_path, output_file, clinical_file=None, norm_method="mrn", design_col=None):
    """Normalize microarray or RNASeq expression data.

    Args:
        input_path (str): Path to the input file (for RNASeq) or directory (for microarray).
        output_file (str): Path to the output file.
        clinical_file (str, optional): Path to the input file containing clinical data.
            Not required for microarray.
        norm_method (str, optional): Normalization method employed if data is RNASeq (tmm or mrn).
            Defaults to "mrn". For microarray, RMA normalization is used.
    """
    if os.path.isdir(input_path):
        normalize_microarray(input_path, output_file, remove_cel_dir=True)
    else:
        normalize_rnaseq(input_path, clinical_file, output_file, norm_method, design_col)


def get_genes_from_file(filename):
    """Read genes from a text file with one gene symbol per line.

    Args:
        filename (str): Path to the text file containing gene symbols.

    Returns:
        list: List of gene symbols.
    """
    with open(filename, 'r', encoding='utf-8') as file:
        genes = [line.strip() for line in file]

    return genes


def get_matrisome_genes(core_matrisome_only=False):
    """Retrieve the human matrisome genes from MSigDB.

    Args:
        core_matrisome_only (bool, optional): If True, only retrieve the core matrisome genes. Defaults to False.

    Returns:
        list: List of matrisome gene symbols.
    """
    if core_matrisome_only:
        url = CORE_MATRISOME_URL
        field = "NABA_CORE_MATRISOME"
    else:
        url = ALL_MATRISOME_URL
        field = "NABA_MATRISOME"

    response = requests.get(url, timeout=10)
    response.raise_for_status()
    matrisome_genes = response.json()[field]["geneSymbols"]
    return matrisome_genes


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
        >>> matrisome_genes = tdq.preprocess.get_matrisome_genes()
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
