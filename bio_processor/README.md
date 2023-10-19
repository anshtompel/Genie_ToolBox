# BIO_FILES_PROCESSOR

> This README describes tool for manipulation with bioformats such as FASTA, GBK and BLAST results.

## 1.Installation 

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
import bio_processor.bio_files_processor
```
and call one of four available functions:
   -  `select_genes_from_gbk_to_fasta` - selects a certain number of genes before and after each of the genes of interest and save in .fasta file;
   -  `convert_multiline_fasta_to_oneline` - converts fasta with multiple sequence lines to .fasta file with a single sequence line; 
   -  `change_fasta_start_pos` - shifts start position of sequence;
   -  `parse_blast_output` - parses blast file and creates file only with the best matches.

## 2. Description of operations
### Selection genes from GBK file 
Function `select_genes_from_gbk_to_fasta` selects a certain number of genes before and after each of the genes of interest and save their protein sequence (translation) from .gbk file and save it in a .fasta file.
Input:
    - `input_gbk` (str): name of input .gbk file. Should be located in working directory.
    - `genes` (str): gene of interest.
    - `output_fasta` (str): name of output file. Default - **'file.fasta'**.
    - `n_before` (int): the number of genes before gene of interest. Default - *1*.
    - `n_after` (int): the number of genes after gene of interest. Default - *1*.
    
Output:
Function creates .fasta file with gene names and their translations except of gene of interest. If the function worked without errors, it will be return *Fasta file was written in `your_working_directory`*.

**NOTES**
1) Function can accept only one gene in the same time.
2) It is requered that input file is located in the tool directory. In the same directory function will be created output file in .fasta format.

### Multiline FASTA to single line 
Function `convert_multiline_fasta_to_oneline` converts .fasta sequence with multiple sequence lines to .fasta file with a single sequence line. The function combines strings with a sequences by removing control character `\n` at the end of the lines of interest.
Input:
- `input_fasta` (str): name of input .fasta file. Should be located in working directory.
- `output_fasta` (Optional[str] = None): name of output file. Default - *None*. If file name is not passed, **'file.fasta'** name will be used.
Output:
The function does not return any value in stadart output. Creates new file in working directory.

### Change FASTA start
Function `change_fasta_start_pos` shifts start position of sequence forward (+n), or back (-n).
Input:
- `input_fasta` (str): name of input .fasta file. Should be located in working directory.
- `shift` (int): an integer (can be negative) is how much sequence should be shifted.
- `output_fast` (Optional[str] = None): output file name. Default - None. If file name is not passed, **'shifted_fasta.fasta'** name will be used.
Output:
The function does not return any value in stadart output. Creates new file in working directory.

### Get the best match from BLAST
Function `parse_blast_output` Parses blast file in .txt format and creates output .txt file only with the best match in *Description* row from every *Query*.
Input:
- `input_file` (str): name of input .txt file. Should be located in working directory.
- `output_file` (Optional[str] = None): output file name. Default - None. If file name is not passed, 'sorted_blast_result.txt' name will be used.
Output:
The function does not return any value in stadart output. Creates new file in working directory.
