from pathlib import Path
import re
from urllib.error import HTTPError
import pkg_resources

from Bio import Entrez
import pandas as pd


def get_entrez_email():
    """
    Retrieve the email for NCBI API.

    Returns
    -------
    str
        The email address read from the email_for_ncbi_tracking.txt file within the package.

    """
    email_file = pkg_resources.resource_filename('mypackage', 'email_for_ncbi_tracking.txt')
    with open(email_file, 'r', encoding='utf-8') as file:
        email = file.read().strip()
    return email
