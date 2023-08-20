<!-- markdownlint-disable -->


<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/geo.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `transcriptomic_data_integrator.geo`


<details class="collapsible-section" style="border: 1px solid #ccc; border-radius: 5px; padding: 10px;">
  <summary class="collapsible-title" style="cursor: pointer; color: #007bff; font-weight: bold; margin: -10px; padding: 10px;">Expand section</summary>





---

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/geo.py#L18"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_entrez_email`

```python
get_entrez_email()
```

Retrieve the email for NCBI API. 



**Returns:**
 
 - <b>`str`</b>:  The email address read from the email_for_ncbi_tracking.txt file within the package. 


---

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/geo.py#L32"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/geo.py#L51"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/geo.py#L83"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/geo.py#L109"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/geo.py#L127"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/geo.py#L160"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/geo.py#L180"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `search_geo`

```python
search_geo(
    query,
    db='gds',
    max_results=25,
    exception_on_http_error=False,
    warn_on_http_error=True
)
```

Retrieve a list of GEO identifiers given a search query. 



**Args:**
 
 - <b>`query`</b> (str):  The search query string. 
 - <b>`db`</b> (str, optional):  The database to search. Defaults to "gds." 
 - <b>`max_results`</b> (int, optional):  The maximum number of results to return. Defaults to 25. 
 - <b>`exception_on_http_error`</b> (bool, optional):  If True, raise an exception on HTTP error. Defaults to False. 
 - <b>`warn_on_http_error`</b> (bool, optional):  If True, print a warning on HTTP error. Defaults to True. 



**Returns:**
 
 - <b>`list`</b>:  List of GEO identifiers corresponding to the query. 


---

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/geo.py#L210"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `download_geo_expression_data`

```python
download_geo_expression_data(gse: GSE, output_dir=None, timeout=10)
```

Download raw microarray data or RNASeq counts from a GEO accession. 



**Args:**
 
 - <b>`gse`</b> (GEOparse.GEOTypes.GSE):  The GEO series object. 
 - <b>`output_dir`</b> (str, optional):  The directory to save the raw data.  Defaults to None (save to current working directory). 
 - <b>`timeout`</b> (int, optional):  The timeout in seconds for the HTTP request. Defaults to 10. 


---

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/geo.py#L252"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_geo_clinical_characteristics`

```python
get_geo_clinical_characteristics(gse: GSE, output_file=None)
```

Parse clinical data from a GEO accession. 



**Args:**
 
 - <b>`gse`</b> (GEOparse.GEOTypes.GSE):  The GEO series object. 
 - <b>`output_file`</b> (str, Optional):  The file to save the clinical data. Defaults to None.  If None, file is saved to {accession}_clinical_data.tsv in current working directory. 


---

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/geo.py#L318"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/geo.py#L344"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `clean_gpl_annotation_column_values`

```python
clean_gpl_annotation_column_values(annotation_column: Series) → Series
```

Ensure all values in the annotation column are strings using ' // ' as separator. 



**Args:**
 
 - <b>`annotation_column`</b> (pandas.Series):  The annotation column. 

**Returns:**
 
 - <b>`pandas.Series`</b>:  The cleaned annotation column. 


---

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/geo.py#L357"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_gene_mapper`

```python
get_gene_mapper(gpl: GPL) → dict
```

raise exception if annotation not parsable 


---

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/geo.py#L405"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `map_probes_to_genes`

```python
map_probes_to_genes(expression_df, gse: GSE)
```

Map probes to genes. The identifiers used for genes will either be symbols,  Entrez IDs, or Ensembl IDs, depending on what the platform annotation table contains. 



**Args:**
 
 - <b>`expression_df`</b> (pandas.DataFrame):  Expression data. 
 - <b>`gse`</b> (GEOparse.GEOTypes.GSE):  The GEO series object. 



**Returns:**
 
 - <b>`pandas.DataFrame`</b>:  Expression data with probes mapped to genes. 



**Notes:**

> This function maps probes to genes using the platform annotation, then aggregates the expression data for each gene using a weighted average. The weights are calculated as 1 / n, where n is the number of genes associated with each probe. This is performed to avoid biasing the average towards probes with more genes. 


---

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/geo.py#L452"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `extract_gsm`

```python
extract_gsm(column_name: str)
```

Extract a GSM sample name from a given string, or return the original string if not found. 


---

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/geo.py#L457"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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


<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/preprocess.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `transcriptomic_data_integrator.preprocess`

<details class="collapsible-section" style="border: 1px solid #ccc; border-radius: 5px; padding: 10px;">
  <summary class="collapsible-title" style="cursor: pointer; color: #007bff; font-weight: bold; margin: -10px; padding: 10px;">Expand section</summary>



Transcriptomic data preprocessing module. 

Functions: 
- normalize_microarray: Normalize microarray expression data in a directory containing CEL.gz files. 
- normalize_rnaseq: Normalize RNA-seq expression data given a file containing raw counts. 
- normalize: Normalize microarray or RNASeq expression data. 
- load_genes_from_file: Read genes from a text file with one gene symbol per line. 
- get_genes_from_msig_set: Fetches genes associated with a given gene set name from MSigDB. 
- convert_genes: Converts a list of genes between different gene identifier formats. 
- select_rows: Select rows in a DataFrame. 
- drop_nan_row_indices: Drop rows where the row index is NaN in an expression matrix. 
- clean_clinical_data: Get filtered and cleaned clinical data table based on a filter specification. 
- join_expression_matrices: Concatenate two or more expression matrices with the same row names. 
- join_and_batch: Join expression matrices, join clinical data tables, and assign batches. 
- batch_correction: Perform batch correction on expression data. 


---

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/preprocess.py#L42"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `normalize_microarray`

```python
normalize_microarray(input_dir, output_file, remove_cel_dir=False)
```

Normalize microarray expression data given a directory containing CEL.gz files. 



**Args:**
 
 - <b>`input_dir`</b> (str):  Path to the directory containing CEL.gz files. 
 - <b>`output_file`</b> (str):  Path to the output file. 
 - <b>`remove_cel_dir`</b> (bool, optional):  If True, remove the input directory after normalization.  Defaults to False. 


---

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/preprocess.py#L61"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `normalize_rnaseq`

```python
normalize_rnaseq(expression_file, clinical_file, output_file)
```

Normalize RNA-seq expression data given a file containing raw counts. 



**Args:**
 
 - <b>`expression_file`</b> (str):  Path to the input file containing raw counts. 
 - <b>`clinical_file`</b> (str):  Path to the input file containing clinical data. 
 - <b>`output_file`</b> (str):  Path to the output file. 


---

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/preprocess.py#L76"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `normalize`

```python
normalize(input_path, output_file, clinical_file=None)
```

Normalize microarray or RNASeq expression data. 



**Args:**
 
 - <b>`input_path`</b> (str):  Path to the input file (for RNASeq) or directory (for microarray). 
 - <b>`output_file`</b> (str):  Path to the output file. 
 - <b>`clinical_file`</b> (str, optional):  Path to the input file containing clinical data.  Not required for microarray. 


---

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/preprocess.py#L91"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `load_genes_from_file`

```python
load_genes_from_file(filename)
```

Read genes from a text file with one gene symbol per line. 



**Args:**
 
 - <b>`filename`</b> (str):  Path to the text file containing gene symbols. 



**Returns:**
 
 - <b>`list`</b>:  List of gene symbols. 


---

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/preprocess.py#L106"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_genes_from_msig_set`

```python
get_genes_from_msig_set(gene_set_name, species='human')
```

Fetches the genes associated with a given gene set name from the MSigDB (e.g, NABA_MATRISOME). 

This function constructs a URL for the specified gene set name and species, then performs a GET request to fetch the associated genes in JSON format from the Molecular Signatures Database (MSigDB). 



**Args:**
 
 - <b>`gene_set_name`</b> (str):  The name of the gene set for which to fetch the associated genes. 
 - <b>`species`</b> (str, optional):  The species for which to fetch the gene set. Defaults to "human". 



**Returns:**
 
 - <b>`list[str]`</b>:  A list of gene symbols associated with the specified gene set name. 



**Raises:**
 
 - <b>`HTTPError`</b>:  If the GET request to the MSigDB results in an error. 


---

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/preprocess.py#L128"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `convert_genes`

```python
convert_genes(
    genes: Iterable,
    in_format: str,
    out_format: str,
    species: str = 'human',
    returnall: bool = False
) → Series
```

Converts a list of genes between formats 'entrezgene', 'ensembl.gene', and 'symbol'. 



**Args:**
 
 - <b>`genes`</b> (Union[List, pd.Series]):  A list of genes. 
 - <b>`in_format`</b> (str):  The format of the input genes. 
 - <b>`out_format`</b> (str):  The format of the output genes. 
 - <b>`species`</b> (str, optional):  The species of the genes. Defaults to "human". 
 - <b>`returnall`</b> (bool, optional):  Whether to return return complete lists of duplicate  or missing query terms. Defaults to False. 



**Returns:**
 
 - <b>`pd.Series`</b>:  Query results. Index is the input genes, values are the output genes. 


---

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/preprocess.py#L159"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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
import transcriptomic_data_integrator as tdi
expression_df = pd.DataFrame({"GSM1234": [3.452, 4.123, 5.678, 6.789],
                              "GSM5678": [1.234, 2.345, 3.456, 4.567]})
expression_df.index = ["A1BG", "A2M", "CA10", "SEMA6B"]
expression_df.index.name = "symbol"
matrisome_genes = tdi.preprocess.get_genes_from_msig_set("NABA_MATRISOME")
matrisome_expression_df = tdi.preprocess.select_rows(expression_df, matrisome_genes)
print(matrisome_expression_df)
```
Output:
```plain
             GSM1234 GSM5678
    symbol
        A2M    4.123   2.345
    SEMA6B    6.789   4.567
```

---

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/preprocess.py#L190"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `drop_nan_row_indices`

```python
drop_nan_row_indices(expr_df: DataFrame)
```

Drop rows where the row index is NaN. 



**Args:**
 
 - <b>`expr_df`</b> (pandas.DataFrame):  The expression matrix. 



**Returns:**
 
 - <b>`pandas.DataFrame`</b>:  The expression matrix with NaN rows dropped. 


---

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/preprocess.py#L202"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `clean_clinical_data`

```python
clean_clinical_data(
    clinical_df: DataFrame,
    specification: dict,
    ignore_case: bool = True,
    drop_no_match_samples: bool = True
) → (<class 'DataFrame'>, <class 'DataFrame'>)
```

Get filtered and cleaned clinical data table based on a filter specification.  The specification contains column names the patterns to extract values for each column.  When multiple patterns match, (e.g. "disease" and "no disease"), the longer match is chosen. 



**Args:**
 
 - <b>`clinical_df`</b> (pd.DataFrame):  The clinical data table. 
 - <b>`specification`</b> (dict):  The filter specification. 
 - <b>`Example`</b>:  `specification={'condition': ['tumor', 'normal'], 'age': [r'\d+']}` 
 - <b>`ignore_case`</b> (bool, optional):  Whether to ignore case when matching patterns. Defaults to True. 
 - <b>`drop_no_match_samples`</b> (bool, optional):  Whether to drop samples that do not match any patterns  for any column. Defaults to True. 



**Returns:**
 
 - <b>`pd.DataFrame`</b>:  The filtered and cleaned clinical data table. 



**Raises:**
 
 - <b>`ValueError`</b>:  If a column name in the specification is not in the clinical data table. 
 - <b>`ValueError`</b>:  If a pattern in the specification is not a valid regular expression. 
 - <b>`ValueError`</b>:  If no patterns match one of the values in a column and  



**Example:**
```python
import transcriptomic_data_integrator as tdi
import pandas as pd
clinical_df = pd.DataFrame({
    "condition": ["cns tumor tissue", "normal tissue", "normal tissue", "tumor tissue"],
    "age": ["45", "50 yo", "eta: 55 anni", "75 years old"],
    "organism": ["homosapiens", "homosapiens", "homosapiens", "homosapiens"]},
    index=["sample1", "sample2", "sample3", "sample4"])
clinical_df.index.name = "sample_name"
specification = {'condition': ['tumor', 'normal'], 'patient_age': [r'\d+']}
cleaned_clinical_df = tdi.preprocess.clean_clinical_data(clinical_df, specification)
print(cleaned_clinical_df)
```
Output:
```plain
                 condition patient_age
     sample_name                      
     sample1         tumor          45
     sample2        normal          50
     sample3        normal          55
     sample4         tumor          75
```

---

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/preprocess.py#L282"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `join_expression_matrices`

```python
join_expression_matrices(expression_dataframes: List[DataFrame])
```

Concatenate two expression matrices with the same row names. 



**Args:**
 
 - <b>`expression_dataframes`</b> (list[pd.DataFrame]):  A list of expression matrices. 



**Returns:**
 
 - <b>`pd.DataFrame`</b>:  The concatenated expression matrix. 



**Raises:**
 
 - <b>`ValueError`</b>:  If the expression matrices do not all have the same index (row names). 


---

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/preprocess.py#L303"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `join_and_batch`

```python
join_and_batch(
    expression_dataframes: List[DataFrame],
    clinical_dataframes: List[DataFrame]
) → (<class 'DataFrame'>, <class 'DataFrame'>)
```

Join expression matrices, join the corresponding clinical data tables, and assign batches. 

 This function adds a batch column to the joined clinical table and assigns batch numbers.  For instance, if we have `batches([expr_df1, expr_df2], [clin_df1, clin_df2])`,  samples in `expr_df1`/`clin_df1` are batch 1 and samples in `expr_df2`/`clin_df2` are batch 2. 



**Args:**
 
     - <b>`expression_dataframes`</b> (list[pd.DataFrame]):  A list of expression matrices. 
     - <b>`clinical_dataframes`</b> (list[pd.DataFrame]):  A list of clinical data tables. 



**Returns:**
 
     - <b>`pd.DataFrame`</b>:  The joined expression matrix. 
     - <b>`pd.DataFrame`</b>:  The joined clinical data table. 



**Raises:**
 
     - <b>`ValueError`</b>:  If the column names in an expression matrix do not match the index (row names)  in the clinical data table. 
     - <b>`ValueError`</b>:  If the clinical data tables do not all have the same column names. 



**Example:**
```python
import transcriptomic_data_integrator as tdi
import pandas as pd
# Dummy expression dataframes
expr_df1 = pd.DataFrame({"sample1": [1, 2], "sample2": [3, 4]},
                        index=["gene1", "gene2"])
expr_df2 = pd.DataFrame({"sample3": [5, 6], "sample4": [7, 8]},
                        index=["gene1", "gene2"])
expr_df3 = pd.DataFrame({"sample5": [9, 10], "sample6": [11, 12]},
                        index=["gene1", "gene2"])
# Dummy clinical dataframes
clinical_df1 = pd.DataFrame({"condition": ["tumor", "normal"]},
                            index=["sample1", "sample2"])
clinical_df2 = pd.DataFrame({"condition": ["normal", "tumor"]},
             index=["sample3", "sample4"])
clinical_df3 = pd.DataFrame({"condition": ["tumor", "tumor"]},
                            index=["sample5", "sample6"])
# Lists of expression and clinical dataframes
expression_dataframes = [expr_df1, expr_df2, expr_df3]
clinical_dataframes = [clinical_df1, clinical_df2, clinical_df3]
# Using the batches function
expr_joined, clinical_joined = tdi.preprocess.batches(expression_dataframes,
                                                      clinical_dataframes)
print("Joined Expression Matrix:")
print(expr_joined)
print("Joined Clinical Data Table:")
print(clinical_joined)
```
Output:
```plain
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
```



---

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/preprocess.py#L389"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `batch_correction`

```python
batch_correction(
    expression_file: str,
    clinical_file: str,
    variable: str,
    data_type: str,
    factor_levels: List[str],
    output_file: str
)
```

Perform batch correction on expression data. 



**Args:**
 
 - <b>`expression_file`</b> (str):  Path to the input file containing expression data. 
 - <b>`clinical_file`</b> (str):  Path to the input file containing clinical data. 
 - <b>`variable`</b> (str):  The column in the clinical data table to use for batch correction. 
 - <b>`data_type`</b> (str):  The data type of the expression data (microarray or RNASeq). 
 - <b>`factor_levels`</b> (list[str]):  The factor levels of the variable to use for batch correction.  This is used to validate that the clinical data table has the correct factor levels. 
 - <b>`output_file`</b> (str):  Path to the output file. 



**Raises:**
 
 - <b>`ValueError`</b>:  If the clinical data table does not have the correct factor levels. 
 - <b>`ValueError`</b>:  If the data type is not "microarray" or "RNASeq". 
 - <b>`ValueError`</b>:  If the variable is not in the clinical data table. 
 - <b>`ValueError`</b>:  If the column names in the expression data do not match the row names  in the clinical data table. 




</details>


<hr style="border:2px solid gray">


<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/util.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `transcriptomic_data_integrator.util`


<details class="collapsible-section" style="border: 1px solid #ccc; border-radius: 5px; padding: 10px;">
  <summary class="collapsible-title" style="cursor: pointer; color: #007bff; font-weight: bold; margin: -10px; padding: 10px;">Expand section</summary>




---

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/util.py#L4"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/src/transcriptomic_data_integrator/util.py#L17"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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
