import re
from functools import wraps
from urllib.error import HTTPError
import tarfile
import gzip
import os
import shutil

import pkg_resources
from Bio import Entrez
import GEOparse
import requests


def get_entrez_email():
    """
    Retrieve the email for NCBI API.

    Returns
    -------
    str
        The email address read from the email_for_ncbi_tracking.txt file within the package.

    """
    email_file = pkg_resources.resource_filename('transcriptomics_data_query',
                                                 'email_for_ncbi_tracking.txt')
    with open(email_file, 'r', encoding='utf-8') as file:
        email = file.read().strip()
    return email


def check_entrez_email(func):
    """
    Decorator to check and set the Entrez email if it is None.

    Parameters
    ----------
    func : function
        The function to be decorated.

    Returns
    -------
    wrapper : function
        The wrapped function.
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

    Parameters
    ----------
    geo_identifier : str
        The GEO identifier for the query.
    exception_on_http_error : bool, optional
        If True, raise an exception on HTTP error, default is False.
    warn_on_http_error : bool, optional
        If True, print a warning on HTTP error, default is True.

    Returns
    -------
    str or None
        The corresponding GEO accession if found, else None.

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

    Parameters
    ----------
    geo_accession : str
        The GEO accession for the query.
    exception_on_http_error : bool, optional
        If True, raise an exception on HTTP error, default is False.
    warn_on_http_error : bool, optional
        If True, print a warning on HTTP error, default is True.

    Returns
    -------
    str or None
        The corresponding GEO identifier if found, else None.

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

    Parameters
    ----------
    geo_ids : list of str
        The GEO identifiers for the query.
    exception_on_http_error : bool, optional
        If True, raise an exception on HTTP error, default is False.
    warn_on_http_error : bool, optional
        If True, print a warning on HTTP error, default is True.
    default_accession : NoneType or str, optional
        Default value to use for study accession if it could not be found (e.g. None or "unknown").

    Returns
    -------
    list of str
        The corresponding GEO accessions.

    """
    return [accession_from_id(id, default_accession, exception_on_http_error, warn_on_http_error)
            for id in geo_ids]


@check_entrez_email
def get_study_description(geo_id, exception_on_http_error=False, warn_on_http_error=True):
    """
    Retrieve GEO study description given an identifier.

    Parameters
    ----------
    geo_id : str
        The GEO identifier for the query.
    exception_on_http_error : bool, optional
        If True, raise an exception on HTTP error, default is False.
    warn_on_http_error : bool, optional
        If True, print a warning on HTTP error, default is True.

    Returns
    -------
    str or None
        The corresponding study description if found, else None.

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

    Parameters
    ----------
    geo_id : str
        The GEO identifier for the query.
    exception_on_http_error : bool, optional
        If True, raise an exception on HTTP error, default is False.
    warn_on_http_error : bool, optional
        If True, print a warning on HTTP error, default is True.
    default_accession : NoneType or str, optional
        Default value to use for study accession if it could not be found (e.g. None or "unknown").

    Returns
    -------
    str or None
        The corresponding study description if found, else None.

    """
    if convert_to_accessions:
        return {accession_from_id(geo_id, default_accession): get_study_description(geo_id)
                for geo_id in geo_study_ids}

    return {geo_id: get_study_description(geo_id) for geo_id in geo_study_ids}


@check_entrez_email
def search_geo(query, db="gds", max_results=25, exception_on_http_error=False,
               warn_on_http_error=True, print_n_results=True):
    """
    Retrieve a list of GEO identifiers given a search query.

    Parameters
    ----------
    query : str
        The search query string.
    db : str, optional
        The database to search, default is "gds."
    exception_on_http_error : bool, optional
        If True, raise an exception on HTTP error, default is False.
    warn_on_http_error : bool, optional
        If True, print a warning on HTTP error, default is True.

    Returns
    -------
    list
        List of GEO identifiers corresponding to the query.

    """
    try:
        handle = Entrez.esearch(db=db, term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        if print_n_results:
            print(f"Hits: {record['Count']}")
        return record["IdList"]
    except HTTPError as http_err:
        if exception_on_http_error:
            raise http_err
        if warn_on_http_error:
            print(f"HTTP error retrieving GEO identifiers for query: {query}")
        return []


def download_raw_data(accession, output_dir=".", timeout=300):
    """
    Parameters
    ----------
    accession : str
        The GEO accession.
    output_dir : str
        The directory to save the raw data.
    timeout : int
        The timeout in seconds for the HTTP request.
    Returns
    -------
    list
        Paths to decompressed CEL files.

    """
    url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={accession}&format=file"
    file_path = os.path.join(output_dir, f"{accession}.tar")

    with requests.get(url, stream=True, timeout=timeout) as response:
        with open(file_path, 'wb') as file:
            shutil.copyfileobj(response.raw, file)

    with tarfile.open(file_path, 'r') as tar:
        tar.extractall(path=accession)

    cel_files = []
    for root, _, files in os.walk(accession):
        for file in files:
            if file.lower().endswith('.cel.gz'):
                gz_path = os.path.join(root, file)
                cel_path = os.path.splitext(gz_path)[0]
                with gzip.open(gz_path, 'rb') as gz_file:
                    with open(cel_path, 'wb') as cel_file:
                        shutil.copyfileobj(gz_file, cel_file)
                cel_files.append(cel_path)

    return cel_files


def map_probes_to_genes(expression_data, accession):
    """
    Parameters
    ----------
    expression_data : pandas.DataFrame
        Expression data.
    accession : str
        The GEO accession.
    Returns
    -------
    pandas.DataFrame
        Expression data with probes mapped to genes.

    """
    gse = GEOparse.get_GEO(geo=accession)
    platform_id = gse.gpls[list(gse.gpls.keys())[0]].get_accession()
    gpl = GEOparse.get_GEO(geo=platform_id)
    annotation = {row['ID']: row['Gene Symbol'] for _, row in gpl.table.iterrows()}
    expression_data.index = expression_data.index.map(annotation)

    # TODO: Drop NaNs and resolve duplicates

    return expression_data
