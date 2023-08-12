<!-- markdownlint-disable -->

<style>
  .collapsible-section {
    border: 1px solid #ccc;
    border-radius: 5px;
    padding: 10px;
  }

  .collapsible-title {
    cursor: pointer;
    color: #007bff;
    font-weight: bold;
    margin: -10px;
    padding: 10px;
  }
</style>

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/geo.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `transcriptomics_data_query.geo`


<details class="collapsible-section">
  <summary class="collapsible-title">Expand section</summary>


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


</details>


<hr style="border:2px solid gray">


<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/preprocess.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>


# <kbd>module</kbd> `transcriptomics_data_query.preprocess`


<details class="collapsible-section">
  <summary class="collapsible-title">Expand section</summary>

---

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/preprocess.py#L4"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `rma_normalization_r`

```python
rma_normalization_r(input_dir, output_file, remove_cel_dir=False)
```

Perform RMA normalization on raw data (directory containing CEL files). 



**Args:**
 
 - <b>`input_dir`</b> (str):  Path to the directory containing CEL files. 
 - <b>`output_file`</b> (str):  Path to the output file. 
 - <b>`remove_cel_dir`</b> (bool, optional):  If True, remove the directory containing CEL files after normalization. Defaults to False. 


---

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/preprocess.py#L14"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_genes_from_file`

```python
get_genes_from_file(filename)
```

Read genes from a text file with one gene symbol per line. 



**Args:**
 
 - <b>`filename`</b> (str):  Path to the text file containing gene symbols. 



**Returns:**
 
 - <b>`list`</b>:  List of gene symbols. 


---

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/preprocess.py#L25"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_matrisome_genes`

```python
get_matrisome_genes(core_matrisome_only=False)
```

Retrieve the matrisome genes from the MSigDB. 



**Args:**
 
 - <b>`core_matrisome_only`</b> (bool, optional):  If True, only retrieve the core matrisome genes. Defaults to False. 



**Returns:**
 
 - <b>`list`</b>:  List of matrisome gene symbols. 


---

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/preprocess.py#L36"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `select_rows`

```python
select_rows(df, values, column=None)
```

Select rows in DataFrame. 



**Args:**
 
 - <b>`df`</b> (pandas.DataFrame):  The data frame. 
 - <b>`values`</b> (list):  The values to select (e.g., ["A2M","A2ML1","ABI3BP"]). 
 - <b>`column`</b> (str, optional):  The column name (e.g., "gene_symbol"). If None, the index is used. 



**Returns:**
 
 - <b>`pandas.DataFrame`</b>:  The selected rows. 



**Example:**

```python
import transcriptomics_data_query as tdq
>>> expression_df = pd.DataFrame({"GSM1234": [3.452, 4.123, 5.678, 6.789],
                                  "GSM5678": [1.234, 2.345, 3.456, 4.567]})
>>> expression_df.index = ["A1BG", "A2M", "CA10", "SEMA6B"]
>>> expression_df.index.name = "symbol"
>>> matrisome_genes = tdq.preprocess.get_matrisome_genes()
>>> matrisome_expression_df = tdq.preprocess.select_rows(expression_df, matrisome_genes)
>>> matrisome_expression_df
         GSM1234 GSM5678
symbol
   A2M     4.123   2.345
SEMA6B     6.789   4.567
```


---

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/preprocess.py#L63"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `drop_nan_row_indices`

```python
drop_nan_row_indices(expr_df: DataFrame)
```

Drop rows where the row index is NaN. 



**Args:**
 
 - <b>`expr_df`</b> (pandas.DataFrame):  The expression matrix. 



**Returns:**
 
 - <b>`pandas.DataFrame`</b>:  The expression matrix with NaN rows dropped. 



</details>


<hr style="border:2px solid gray">


<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/util.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `transcriptomics_data_query.util`


<details class="collapsible-section">
  <summary class="collapsible-title">Expand section</summary>



---

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/util.py#L4"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `is_valid_tar_member`

```python
is_valid_tar_member(member: TarInfo, target_dir: str) → bool
```

Check if a tar member is safe to extract. `target_dir` should be an absolute path. 



**Args:**
 
 - <b>`member`</b> (tarfile.TarInfo):  The tar member to check. 
 - <b>`target_dir`</b> (str):  The absolute path to the target directory. 



**Returns:**
 
 - <b>`bool`</b>:  True if the tar member is safe to extract, False otherwise. 


---

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/util.py#L17"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `extract_tar`

```python
extract_tar(tar_file, target_dir, delete_tar=False)
```

Extract a tar file to a target directory. 



**Args:**
 
 - <b>`tar_file`</b> (str):  Path to the tar file. 
 - <b>`target_dir`</b> (str):  Path to the target directory. 
 - <b>`delete_tar`</b> (bool, optional):  If True, delete the tar file after extraction. Defaults to False. 


</details>

---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._