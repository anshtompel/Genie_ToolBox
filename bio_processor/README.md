# BIO_FILES_PROCESSOR

> This README describes tool for manipulation with bioformats such as FASTA, GBK and BLAST results.

## Content

1. [Installation](#installation)
2. [Description of operations](#description)
    * [Selection genes from GBK file](#selectgbk)
    * [Multiline FASTA to single line](#multiline)
    * [Change FASTA start](#change_start)
3. [Contacts](#contacts)

## 1.Installation <a name="installation"></a>

1. Clone Genie_ToolBox repo to your local machine using SSH:
``` bash
git clone git@github.com:anshtompel/Genie_ToolBox.git
```
or HTTPS:
``` bash
git clone https://github.com/anshtompel/Genie_ToolBox.git
```
2. Change directory to `bio_processor` and launch `bio_files_processor.py` or launch script as:
```python
from bio_processor import bio_files_processor
```
and call one of four available functions:
   -  `select_genes_from_gbk_to_fasta` 
   -  `convert_multiline_fasta_to_oneline`
   -  `change_fasta_start_pos`
   -  `parse_blast_output`


