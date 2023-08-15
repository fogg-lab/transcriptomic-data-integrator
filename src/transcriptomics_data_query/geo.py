import re
from typing import Union
from functools import wraps
import urllib.request
from urllib.error import HTTPError
import os
import shutil

import pkg_resources
import numpy as np
import pandas as pd
from Bio import Entrez
import GEOparse
import requests

from . import util

def get_entrez_email():
    """
    Retrieve the email for NCBI API.

    Returns:
        str: The email address read from the email_for_ncbi_tracking.txt file within the package.
    """
    email_file = pkg_resources.resource_filename('transcriptomics_data_query',
                                                 'email_for_ncbi_tracking.txt')
    with open(email_file, 'r', encoding='utf-8') as file:
        email = file.read().strip()
    return email


def check_entrez_email(func):
    """
    Decorator to check and set the Entrez email if it is None.

    Args:
        func (function): The function to be decorated.

    Returns:
        function: The wrapped function.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        if Entrez.email is None:
            Entrez.email = get_entrez_email()
        return func(*args, **kwargs)

    return wrapper


@check_entrez_email
def accession_from_id(geo_identifier, default_accession=None, exception_on_http_error=False,
                      warn_on_http_error=True):
    """
    Retrieve GEO accession given a GEO identifier.

    Args:
        geo_identifier (str): The GEO identifier for the query.
        exception_on_http_error (bool, optional): If True, raise an exception on HTTP error. Defaults to False.
        warn_on_http_error (bool, optional): If True, print a warning on HTTP error. Defaults to True.

    Returns:
        str or None: The corresponding GEO accession if found, else None.
    """
    try:
        handle = Entrez.efetch(db="gds", id=geo_identifier, retmode="text")
        text = handle.read()
        handle.close()
    except HTTPError as http_err:
        if exception_on_http_error:
            raise http_err
        if warn_on_http_error:
            print(f"HTTP error retrieving accession for GEO ID: {geo_identifier}")
        return default_accession

    # Use regular expression to find accession
    match = re.search(r'Accession: (GSE\d+)', text)
    if match:
        return match.group(1)
    return default_accession


@check_entrez_email
def id_from_accession(geo_accession, exception_on_http_error=False, warn_on_http_error=True):
    """
    Retrieve GEO identifier given a GEO accession.

    Args:
        geo_accession (str): The GEO accession for the query.
        exception_on_http_error (bool, optional): If True, raise an exception on HTTP error. Defaults to False.
        warn_on_http_error (bool, optional): If True, print a warning on HTTP error. Defaults to True.

    Returns:
        str or None: The corresponding GEO identifier if found, else None.
    """
    try:
        handle = Entrez.esearch(db="gds", term=geo_accession)
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"][0] if record["IdList"] else None
    except HTTPError as http_err:
        if exception_on_http_error:
            raise http_err
        if warn_on_http_error:
            print(f"HTTP error retrieving identifier for GEO accession: {geo_accession}")
        return None


def get_accessions_from_ids(geo_ids, default_accession=None, exception_on_http_error=False,
                            warn_on_http_error=True):
    """
    Retrieve a list of GEO accessions given a list of GEO identifiers.

    Args:
        geo_ids (list of str): The GEO identifiers for the query.
        exception_on_http_error (bool, optional): If True, raise an exception on HTTP error. Defaults to False.
        warn_on_http_error (bool, optional): If True, print a warning on HTTP error. Defaults to True.
        default_accession (NoneType or str, optional): Default value to use for study accession if it could not be found (e.g. None or "unknown").

    Returns:
        list of str: The corresponding GEO accessions.
    """
    return [accession_from_id(id, default_accession, exception_on_http_error, warn_on_http_error)
            for id in geo_ids]


@check_entrez_email
def get_study_description(geo_id, exception_on_http_error=False, warn_on_http_error=True):
    """
    Retrieve GEO study description given an identifier.

    Args:
        geo_id (str): The GEO identifier for the query.
        exception_on_http_error (bool, optional): If True, raise an exception on HTTP error. Defaults to False.
        warn_on_http_error (bool, optional): If True, print a warning on HTTP error. Defaults to True.

    Returns:
        str or None: The corresponding study description if found, else None.
    """
    try:
        handle = Entrez.efetch(db="gds", id=geo_id, retmode="text")
        text = handle.read()
        handle.close()
    except HTTPError as http_err:
        if exception_on_http_error:
            raise http_err
        if warn_on_http_error:
            print(f"HTTP error retrieving info for GEO accession: {accession_from_id(geo_id)}")
        return None

    # Extract study description
    description = re.search(r'\d+\. (.+?)\n', text)

    if description:
        return description.group(1)

    return None


def get_descriptions_from_ids(geo_study_ids, convert_to_accessions=True, default_accession=None):
    """
    Retrieve GEO study description given an identifier.

    Args:
        geo_id (str): The GEO identifier for the query.
        exception_on_http_error (bool, optional): If True, raise an exception on HTTP error. Defaults to False.
        warn_on_http_error (bool, optional): If True, print a warning on HTTP error. Defaults to True.
        default_accession (NoneType or str, optional): Default value to use for study accession if it could not be found (e.g. None or "unknown").

    Returns:
        str or None: The corresponding study description if found, else None.
    """
    if convert_to_accessions:
        return {accession_from_id(geo_id, default_accession): get_study_description(geo_id)
                for geo_id in geo_study_ids}

    return {geo_id: get_study_description(geo_id) for geo_id in geo_study_ids}


@check_entrez_email
def search_geo(query, db="gds", max_results=25, exception_on_http_error=False,
               warn_on_http_error=True):
    """
    Retrieve a list of GEO identifiers given a search query.

    Args:
        query (str): The search query string.
        db (str, optional): The database to search. Defaults to "gds."
        max_results (int, optional): The maximum number of results to return. Defaults to 25.
        exception_on_http_error (bool, optional): If True, raise an exception on HTTP error. Defaults to False.
        warn_on_http_error (bool, optional): If True, print a warning on HTTP error. Defaults to True.

    Returns:
        list: List of GEO identifiers corresponding to the query.
    """
    try:
        handle = Entrez.esearch(db=db, term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        print(f"Hits: {record['Count']}")
        return record["IdList"]
    except HTTPError as http_err:
        if exception_on_http_error:
            raise http_err
        if warn_on_http_error:
            print(f"HTTP error retrieving GEO identifiers for query: {query}")
        return []


def download_geo_expression_data(gse: GEOparse.GEOTypes.GSE, output_dir=None, timeout=10):
    """
    Download raw microarray data or RNASeq counts from a GEO accession.

    Args:
        gse (GEOparse.GEOTypes.GSE): The GEO series object.
        output_dir (str, optional): The directory to save the raw data.
                                    Defaults to None (save to current working directory).
        timeout (int, optional): The timeout in seconds for the HTTP request. Defaults to 10.
    """
    accession = gse.name
    if output_dir is None:
        output_dir = os.path.join(os.getcwd(), accession)

    # Create output directory if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    url = gse.metadata["supplementary_file"][0]
    is_microarray = url.endswith(".tar")

    if is_microarray:
        file_path = os.path.join(output_dir, f"{accession}.tar")
    else:
        file_path = os.path.join(output_dir, f"{accession}.txt.gz")

    with urllib.request.urlopen(url) as response:
        # Read the content and write to the file
        with open(file_path, 'wb') as file:
            shutil.copyfileobj(response, file)

    if is_microarray:
        util.extract_tar(file_path, output_dir, delete_tar=True)
    else:
        # extract gz
        import gzip
        with gzip.open(file_path, 'rb') as f_in:
            out_path = os.path.join(output_dir, f"{accession}_expression_matrix.tsv")
            with open(out_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)


def get_geo_clinical_characteristics(gse: GEOparse.GEOTypes.GSE, output_file=None):
    """
    Parse clinical data from a GEO accession.

    Args:
        gse (GEOparse.GEOTypes.GSE): The GEO series object.
        output_file (str, Optional): The file to save the clinical data. Defaults to None.
            If None, file is saved to {accession}_clinical_data.tsv in current working directory.
    """

    characteristics = {sample: gse.gsms[sample].metadata["characteristics_ch1"]
                       for sample in gse.gsms}
    clinical_df = pd.DataFrame.from_dict(characteristics, orient='index')

    # Melt the DataFrame to a long format
    clinical_df = clinical_df.melt(ignore_index=False)

    # Split value into key and value parts
    key_values = clinical_df['value'].str.split(': ', expand=True)

    # Extract the characteristic names and format them
    columns = key_values[0].str.lower().str.replace(' ', '_').unique()

    # Pivot the DataFrame
    clinical_df = key_values.pivot_table(index=clinical_df.index, columns=0, values=1, aggfunc='first')

    # Rename the index to "sample_id" and assign characteristic names to the columns
    clinical_df.index.rename("sample_id", inplace=True)
    clinical_df.columns = columns

    # Add other potentially useful metadata
    for field in ["title", "description", "source_name_ch1"]:
        try:
            field_values = {sample: " // ".join(gse.gsms[sample].metadata[field])
                            for sample in gse.gsms}
            clinical_df[field] = pd.Series(field_values)
        except KeyError:
            pass

    # Save the clinical data to a file
    if output_file is None:
        output_file = os.path.join(os.getcwd(), f"{gse.name}_clinical_data.tsv")

    clinical_df.to_csv(output_file, sep='\t')


def weighted_average_group(df, weights):
    """
    Aggregates groups of rows in a Pandas DataFrame using a weighted average.

    Args:
        df (pd.DataFrame): Input DataFrame containing the data.
        weights (list): List of weights corresponding to the rows of the DataFrame.

    Returns:
        result (pd.DataFrame): Aggregated DataFrame with weighted averages.
    """
    # Convert weights to a NumPy array for efficient multiplication
    weights_np = np.array(weights)

    # Multiply data by weights
    weighted_data = df.multiply(weights_np, axis=0)

    # Calculate the sum of the weights for each unique index
    sum_weights = pd.Series(weights_np, index=df.index).groupby(level=0).sum()

    # Divide the sum of the weighted data by the sum of weights for each unique index
    result = weighted_data.groupby(level=0).sum().divide(sum_weights, axis=0)

    return result


def clean_gpl_annotation_column_values(annotation_column: pd.Series) -> pd.Series:
    """Ensure all values in the annotation column are strings using ' // ' as separator.

    Args:
        annotation_column (pandas.Series): The annotation column.
    Returns:
        pandas.Series: The cleaned annotation column.
    """
    annotation_column = annotation_column.astype(str)
    annotation_column = annotation_column.str.replace(' /// ', ' // ')
    return annotation_column


def get_gene_mapper(gpl: GEOparse.GEOTypes.GPL) -> dict:
    """raise exception if annotation not parsable"""

    colnames = gpl.table.columns
    simplified_colnames = pd.Series(colnames, index=colnames).str.lower().str.replace(' ', '_')

    gene_symbol_indices = simplified_colnames.index[simplified_colnames == "gene_symbol"]
    if len(gene_symbol_indices) == 1:
        # There is a gene symbol column
        key = gene_symbol_indices[0]
        print(f"Using annotation column {key} for gene symbols")
        genes_series = gpl.table.set_index('ID')[key]
        genes_series = clean_gpl_annotation_column_values(genes_series)
        return genes_series.to_dict()

    entrez_indices = simplified_colnames.index[simplified_colnames == "entrez_id"]
    if len(entrez_indices) == 0:
        entrez_indices = simplified_colnames.index[simplified_colnames == "entrez_gene"]
    if len(entrez_indices) == 0:
        entrez_indices = simplified_colnames.index[simplified_colnames == "entrez_gene_id"]
    if len(entrez_indices) == 1:
        # There is an entrez column
        key = entrez_indices[0]
        print(f"Using annotation column {key} for Entrez IDs")
        genes_series = gpl.table.set_index('ID')[key]
        genes_series = clean_gpl_annotation_column_values(genes_series)
        return genes_series.to_dict()

    # Last resort: try to find a column with patterns like "ENS[A-Z]+[0-9]+" (Ensembl IDs)
    pattern = re.compile(r'ENS[A-Z]+[0-9]+')
    for colname in colnames:
        first_five_values = gpl.table[colname].head(5)
        col_is_junk = False
        for value in first_five_values:
            if not pattern.findall(str(value)):
                col_is_junk = True
                break
        if not col_is_junk:
            print(f"Using annotation column {colname} for Ensembl IDs")
            genes_series = gpl.table.set_index('ID')[colname]
            genes_series = genes_series.str.findall(pattern)
            genes_series = genes_series.apply(lambda x: ' // '.join(set(x)))
            return genes_series.to_dict()

    # If above efforts have failed, raise exception.
    raise ValueError("Could not parse the platform annotation table.")


def map_probes_to_genes(expression_df, gse: GEOparse.GEOTypes.GSE):
    """Map probes to genes. The identifiers used for genes will either be symbols,
        Entrez IDs, or Ensembl IDs, depending on what the platform annotation table contains.

    Args:
        expression_df (pandas.DataFrame): Expression data.
        gse (GEOparse.GEOTypes.GSE): The GEO series object.

    Returns:
        pandas.DataFrame: Expression data with probes mapped to genes.

    Notes:
        This function maps probes to genes using the platform annotation, then aggregates the
        expression data for each gene using a weighted average. The weights are calculated
        as 1 / n, where n is the number of genes associated with each probe. This is performed
        to avoid biasing the average towards probes with more genes.
    """

    # Validate that the indices of the expression data are not the default RangeIndex
    if (isinstance(expression_df.index, pd.RangeIndex)
        and expression_df.index.start == 0 and expression_df.index.step == 1):
        raise ValueError("The index of the expression data must be set to the probe IDs. "
                         "If you are reading expression data from a file with pandas.read_csv(), "
                         "try setting the index_col argument to 0.")

    platform_id = gse.gpls[list(gse.gpls.keys())[0]].get_accession()
    gpl = GEOparse.get_GEO(geo=platform_id)

    annotation_mapper = get_gene_mapper(gpl)
    expression_df.index = expression_df.index.map(annotation_mapper)

    expression_df = expression_df[~expression_df.index.isin(["nan", ""])]
    expression_df = expression_df.groupby(expression_df.index).mean()
    expanded_genes = expression_df.index.astype(str).str.split(' // ')

    gene_counts = np.array([len(x) for x in expanded_genes])
    expression_df["weight"] = 1 / gene_counts
    expanded_expression_df = expression_df.reindex(expression_df.index.repeat(gene_counts), method='ffill')
    expanded_genes_flat = [gene for genes in expanded_genes for gene in genes]
    expanded_expression_df.index = expanded_genes_flat

    row_weights = expanded_expression_df["weight"]
    expanded_expression_df.drop(columns=["weight"], inplace=True)

    return weighted_average_group(expanded_expression_df, row_weights)


def extract_gsm(column_name: str):
    """Extract a GSM sample name from a given string, or return the original string if not found."""
    return re.search(r'GSM[0-9]+', column_name).group(0) or column_name


def clean_geo_sample_columns(expr_df: pd.DataFrame):
    """
    Clean the sample columns of a GEO expression matrix.

    Args:
        expr_df (pandas.DataFrame): The expression matrix.

    Returns:
        pandas.DataFrame: The expression matrix with cleaned sample columns.
    """
    return expr_df.rename(columns=extract_gsm)
