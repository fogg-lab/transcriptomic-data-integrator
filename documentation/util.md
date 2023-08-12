<!-- markdownlint-disable -->

<a href="https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/src/transcriptomics_data_query/util.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `transcriptomics_data_query.util`





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




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
