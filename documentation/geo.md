<!-- markdownlint-disable -->

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/geo.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `transcriptomics_data_query.geo`





---

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/geo.py#L3"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_entrez_email`

```python
get_entrez_email()
```

Retrieve the email for NCBI API. 



**Returns:**
 
 - <b>`str`</b>:  The email address read from the email_for_ncbi_tracking.txt file within the package. 


---

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/geo.py#L12"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `check_entrez_email`

```python
check_entrez_email(func)
```

Decorator to check and set the Entrez email if it is None. 



**Args:**
 
 - <b>`func`</b> (function):  The function to be decorated. 



**Returns:**
 
 - <b>`function`</b>:  The wrapped function. 


---

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/geo.py#L24"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `accession_from_id`

```python
accession_from_id(
    geo_identifier,
    default_accession=None,
    exception_on_http_error=False,
    warn_on_http_error=True
)
```

Retrieve GEO accession given a GEO identifier. 



**Args:**
 
 - <b>`geo_identifier`</b> (str):  The GEO identifier for the query. 
 - <b>`exception_on_http_error`</b> (bool, optional):  If True, raise an exception on HTTP error. Defaults to False. 
 - <b>`warn_on_http_error`</b> (bool, optional):  If True, print a warning on HTTP error. Defaults to True. 



**Returns:**
 
 - <b>`str or None`</b>:  The corresponding GEO accession if found, else None. 


---

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/geo.py#L39"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `id_from_accession`

```python
id_from_accession(
    geo_accession,
    exception_on_http_error=False,
    warn_on_http_error=True
)
```

Retrieve GEO identifier given a GEO accession. 



**Args:**
 
 - <b>`geo_accession`</b> (str):  The GEO accession for the query. 
 - <b>`exception_on_http_error`</b> (bool, optional):  If True, raise an exception on HTTP error. Defaults to False. 
 - <b>`warn_on_http_error`</b> (bool, optional):  If True, print a warning on HTTP error. Defaults to True. 



**Returns:**
 
 - <b>`str or None`</b>:  The corresponding GEO identifier if found, else None. 


---

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/geo.py#L53"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_accessions_from_ids`

```python
get_accessions_from_ids(
    geo_ids,
    default_accession=None,
    exception_on_http_error=False,
    warn_on_http_error=True
)
```

Retrieve a list of GEO accessions given a list of GEO identifiers. 



**Args:**
 
 - <b>`geo_ids`</b> (list of str):  The GEO identifiers for the query. 
 - <b>`exception_on_http_error`</b> (bool, optional):  If True, raise an exception on HTTP error. Defaults to False. 
 - <b>`warn_on_http_error`</b> (bool, optional):  If True, print a warning on HTTP error. Defaults to True. 
 - <b>`default_accession`</b> (NoneType or str, optional):  Default value to use for study accession if it could not be found (e.g. None or "unknown"). 



**Returns:**
 
 - <b>`list of str`</b>:  The corresponding GEO accessions. 


---

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/geo.py#L69"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_study_description`

```python
get_study_description(
    geo_id,
    exception_on_http_error=False,
    warn_on_http_error=True
)
```

Retrieve GEO study description given an identifier. 



**Args:**
 
 - <b>`geo_id`</b> (str):  The GEO identifier for the query. 
 - <b>`exception_on_http_error`</b> (bool, optional):  If True, raise an exception on HTTP error. Defaults to False. 
 - <b>`warn_on_http_error`</b> (bool, optional):  If True, print a warning on HTTP error. Defaults to True. 



**Returns:**
 
 - <b>`str or None`</b>:  The corresponding study description if found, else None. 


---

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/geo.py#L83"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_descriptions_from_ids`

```python
get_descriptions_from_ids(
    geo_study_ids,
    convert_to_accessions=True,
    default_accession=None
)
```

Retrieve GEO study description given an identifier. 



**Args:**
 
 - <b>`geo_id`</b> (str):  The GEO identifier for the query. 
 - <b>`exception_on_http_error`</b> (bool, optional):  If True, raise an exception on HTTP error. Defaults to False. 
 - <b>`warn_on_http_error`</b> (bool, optional):  If True, print a warning on HTTP error. Defaults to True. 
 - <b>`default_accession`</b> (NoneType or str, optional):  Default value to use for study accession if it could not be found (e.g. None or "unknown"). 



**Returns:**
 
 - <b>`str or None`</b>:  The corresponding study description if found, else None. 


---

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/geo.py#L98"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `search_geo`

```python
search_geo(
    query,
    db='gds',
    max_results=25,
    exception_on_http_error=False,
    warn_on_http_error=True,
    print_n_results=True
)
```

Retrieve a list of GEO identifiers given a search query. 



**Args:**
 
 - <b>`query`</b> (str):  The search query string. 
 - <b>`db`</b> (str, optional):  The database to search. Defaults to "gds." 
 - <b>`exception_on_http_error`</b> (bool, optional):  If True, raise an exception on HTTP error. Defaults to False. 
 - <b>`warn_on_http_error`</b> (bool, optional):  If True, print a warning on HTTP error. Defaults to True. 



**Returns:**
 
 - <b>`list`</b>:  List of GEO identifiers corresponding to the query. 


---

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/geo.py#L114"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `download_raw_data`

```python
download_raw_data(accession, output_dir=None, timeout=10)
```



**Args:**
 
 - <b>`accession`</b> (str):  The GEO accession. 
 - <b>`output_dir`</b> (str):  The directory to save the raw data. 
 - <b>`timeout`</b> (int):  The timeout in seconds for the HTTP request. 


---

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/geo.py#L123"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `weighted_average_group`

```python
weighted_average_group(df, weights)
```

Aggregates groups of rows in a Pandas DataFrame using a weighted average. 



**Args:**
 
 - <b>`df`</b> (pd.DataFrame):  Input DataFrame containing the data. 
 - <b>`weights`</b> (list):  List of weights corresponding to the rows of the DataFrame. 



**Returns:**
 
 - <b>`result`</b> (pd.DataFrame):  Aggregated DataFrame with weighted averages. 


---

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/geo.py#L136"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `map_probes_to_genes`

```python
map_probes_to_genes(expression_df, accession)
```



**Args:**
 
 - <b>`expression_df`</b> (pandas.DataFrame):  Expression data. 
 - <b>`accession`</b> (str):  The GEO accession. 



**Returns:**
 
 - <b>`pandas.DataFrame`</b>:  Expression data with probes mapped to genes. 



**Notes:**

> This function maps probes to genes using the GPL annotation, then aggregates the expression data for each gene using a weighted average. The weights are calculated as 1 / n, where n is the number of genes associated with each probe. This is performed to avoid biasing the average towards probes with more genes. 


---

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/geo.py#L153"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `extract_gsm`

```python
extract_gsm(column_name: str)
```

Extract a GSM sample name from a given string, or return the original string if not found. 


---

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/geo.py#L157"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `clean_geo_sample_columns`

```python
clean_geo_sample_columns(expr_df: DataFrame)
```

Clean the sample columns of a GEO expression matrix. 



**Args:**
 
 - <b>`expr_df`</b> (pandas.DataFrame):  The expression matrix. 



**Returns:**
 
 - <b>`pandas.DataFrame`</b>:  The expression matrix with cleaned sample columns. 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
