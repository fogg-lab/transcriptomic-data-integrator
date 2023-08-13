import subprocess
import shutil

import pkg_resources
import requests
import pandas as pd


CORE_MATRISOME_URL = "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=NABA_CORE_MATRISOME&fileType=json"
ALL_MATRISOME_URL = "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=NABA_MATRISOME&fileType=json"


def rma_normalization_r(input_dir, output_file, remove_cel_dir=False):
    """Perform RMA normalization on raw data (directory containing CEL files).

    Args:
        input_dir (str): Path to the directory containing CEL files.
        output_file (str): Path to the output file.
        remove_cel_dir (bool, optional): If True, remove the directory containing CEL files after normalization. Defaults to False.
    """
    r_script_path = pkg_resources.resource_filename('transcriptomics_data_query',
                                                    'rscripts/rma_normalization.R')
    subprocess.run(["Rscript", r_script_path, input_dir, output_file], check=True)
    if remove_cel_dir:
        shutil.rmtree(input_dir)


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
