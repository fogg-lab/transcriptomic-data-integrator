<!-- markdownlint-disable -->

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/preprocess.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `transcriptomics_data_query.preprocess`





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




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
