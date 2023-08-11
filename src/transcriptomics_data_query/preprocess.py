import pkg_resources
import subprocess
import requests


CORE_MATRISOME_URL = "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=NABA_CORE_MATRISOME&fileType=json"
ALL_MATRISOME_URL = "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=NABA_MATRISOME&fileType=json"


def rma_normalization_r(input_dir, output_file):
    r_script_path = pkg_resources.resource_filename('transcriptomics_data_query',
                                                    'rscripts/rma_normalization.R')
    subprocess.run(["Rscript", r_script_path, input_dir, output_file], check=True)


def get_genes_from_file(filename):
    """
    Read genes from a text file with one gene symbol per line.

    Parameters
    ----------
    filename : str
        Path to the text file containing gene symbols.

    Returns
    -------
    list
        List of gene symbols.

    """
    with open(filename, 'r', encoding='utf-8') as file:
        genes = [line.strip() for line in file]

    return genes


def get_matrisome_genes(core_matrisome_only=False):
    """
    Retrieve the matrisome genes from the MSigDB.

    Parameters
    ----------
    core_matrisome_only : bool, optional
        If True, only retrieve the core matrisome genes, default is False.

    Returns
    -------
    list
        List of matrisome gene symbols.

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


def select_rows(df, column, values):
    """
    Parameters
    ----------
    df : pandas.DataFrame
        The data frame.
    column : str
        The column name (e.g., "gene_symbol").
    values : list
        The values to select (e.g., ["A2M","A2ML1","ABI3BP"]).
    Returns
    -------
    pandas.DataFrame
        The selected rows.

    Example
    -------
    >>> import transcriptomics_data_query as tdq
    >>> df = pd.DataFrame({"symbol": ["A1BG", "A2M", "CA10", "SEMA6B"]})
    >>> matrisome_genes = tdq.preprocess.get_matrisome_genes()
    >>> df_matrisome = tdq.preprocess.select_rows(df, "symbol", matrisome_genes)
    >>> df_matrisome
        symbol
    0   A2M
    1   SEMA6B

    """
    return df[df[column].isin(values)]
