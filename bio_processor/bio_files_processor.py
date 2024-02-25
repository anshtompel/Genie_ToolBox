import os
from typing import Optional

AMINOACIDS = {'G', 'A', 'V', 'I', 'L', 'P',
              'S', 'T', 'C', 'M', 'N', 'Q', 
              'F', 'W', 'Y', 'D', 'E', 'K', 
              'R', 'H'}

def dict_to_fasta(fasta_dict: dict, output_fasta: Optional[str] = None) -> None:
    """
    Converts dictionary from fasta data to .fasta file. 
    Input:
    - fasta_dict (dict): dictionary from fasta data, key - name of sequence, value - sequence.
    - output_fasta (Optional[str] = None): name of output file. Default - None. If file name is not passed, 'file.fasta' name will be used. 
    Output:
    The function does not return any value. Creates new file in working directory.
    """
    if output_fasta is None:
        output_fasta = 'file.fasta'
    if not output_fasta.endswith('.fasta'):
        output_fasta += '.fasta'
    with open(output_fasta, mode='w') as output:
        for string in fasta_dict.items():
            output.write(string[0] + '\n')
            output.write(string[1] + '\n')

    
def filter_gbk(file: str) -> list:
    """
    Filters input .gbk file and creates list of genes and their translation sequences.
    Input:
    - file (str): name of input .gbk file. Should be located in working directory.
    Output:
    Returns list of genes and their translation sequences.
    """
    with open(file) as lines:
        seq_list = []
        gene = None
        translation = None
        for line in lines:
            line = line.strip()
            if line.startswith('/gene'):
                gene = line.strip()
                gene = '>'+ gene
                seq_list.append(gene.replace('>/gene=', '>gene='))
            if line.startswith('/translation') or set(line).issubset(AMINOACIDS):
                if gene is None:
                    continue
                else:
                    translation = line.strip('/translation="')
                    seq_list.append(translation)
    return seq_list


def get_gene_names(seq_list: list) -> list:
    """
    Creates list of unique gene names from filter_gbk output list of .gbk file.
    Input:
    - seq_list (list): filter_gbk output list of .gbk file.
    Output:
    Returns list of unique gene names.
    """
    names = []
    for seq in seq_list:
        if '>gene=' in seq:
            names.append(seq.strip('>gene="'))
    return names


def find_range_genes(seq_list: list, names: list, gene: str, n_before: int, n_after: int) -> list:
    """
    Creates list of neighboring genes of gene of interest and their translation sequences. The number of neighboring genes is determined by function parameters n_before (upper border) and n_after (lower border).
    
   Input:
   - seq_list (list): list of genes and their translation sequences from input .gbk  file.
   - names (list): list of unique gene names.
   - n_before (int): the number of genes before gene of interest. Default - 1.
   - n_after (int): the number of genes after gene of interest. Default - 1.
   
   Output:
   Returns list of genes and their translation sequences according to parameters.
    """
    upper_section = []
    lower_section = []
    left = names.index(gene) - n_before   
    right = names.index(gene) + n_after   
    upper = seq_list.index('>gene=\"'+names[left]+'"')
    selected_gene = seq_list.index('>gene="'+gene+'"')
    lower = seq_list.index('>gene=\"'+names[right]+'"')   
    for i in range(upper, selected_gene):
        upper_section.append(seq_list[i])
    for i in range(lower, selected_gene, -1):
        lower_section.append(seq_list[i])       
    sum_list_genes =  upper_section + lower_section
    return sum_list_genes


def select_genes_from_gbk_to_fasta(input_gbk: str, genes: str, output_fasta: str = 'file.fasta' , n_before: int = 1, n_after: int = 1) -> str:
    """
    Selects a certain number of genes before and after each of the genes of interest and save their protein sequence (translation) from .gbk file and save it in a .fasta file.
    
    Input:
    - input_gbk (str): name of input .gbk file. Should be located in working directory.
    - genes (str): gene of interest.
    - output_fasta (str): name of output file. Default -'file.fasta'.
    - n_before (int): the number of genes before gene of interest. Default - 1.
    - n_after (int): the number of genes after gene of interest. Default - 1.
    
    Output:
    Creates .fasta file with gene names and their translations except of gene of interest. If the function worked without errors, it will be return Fasta file was written in your_working_directory.
    """
    output_list = [] 
    filter_gene_list = filter_gbk(input_gbk)
    gene_of_interest_list = get_gene_names(filter_gene_list)
    output_list = find_range_genes(filter_gene_list, gene_of_interest_list, genes, n_before, n_after)
    with open (output_fasta, mode = 'w') as output:
        for line in output_list:
            output.write(line + '\n')
    return f'Fasta file was written in {os.path.join(os. getcwd(), output_fasta)}.'


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: Optional[str] = None) -> None:
    """
    Converts fasta with multiple sequence lines to .fasta file with a single sequence line. 
    Input:
    - input_fasta (str): name of input .fasta file. Should be located in working directory.
    - output_fasta (Optional[str] = None): name of output file. Default - None. If file name is not passed, 'file.fasta' name will be used.
    Output:
    The function does not return any value. Creates new file in working directory.
    """
    with open (input_fasta) as fasta_file:
        seq = []
        fasta_dict = {}
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
                read_id = line
            else:
                seq.append(line)
            join_str = ''.join(seq)
            fasta_dict[read_id] = join_str
        return dict_to_fasta(fasta_dict, output_fasta)


def change_fasta_start_pos(input_fasta: str, shift: int, output_fasta: Optional[str] = None) -> None:
    """
    Shifts start position of sequence according to input paraments. Sequence could be shifted forward (+int), or back (-int).
    Input:
    - input_fasta (str): name of input .fasta file. Should be located in working directory.
    - shift (int): an integer (can be negative) is how much sequence should be shifted.
    - output_fast (Optional[str] = None): output file name. Default - None. If file name is not passed, 'shifted_fasta.fasta' name will be used.
    Output:
    The function does not return any value. Creates new file in working directory.
    """
    if output_fasta is None:
        output_fasta = 'shifted_fasta.fasta' 
    with open(input_fasta) as fasta_sequence:
        sequence_ID = fasta_sequence.readline().strip()
        sequence = fasta_sequence.readline().strip()
        shifted_seq = sequence[shift:]+sequence[:shift]
    with open(output_fasta, mode = 'w') as output:
        output.write(sequence_ID + '\n')
        output.write(shifted_seq + '\n')
        

def parse_blast_output(input_file: str, output_file: Optional[str] = None) -> None:
    """
    Parses blast file in .txt format and creates output .txt file only with the best match in "Description" row from every "Query".
    Input:
    - input_file (str): name of input .txt file. Should be located in working directory.
    - output_file (Optional[str] = None): output file name. Default - None. If file name is not passed, 'sorted_blast_result.txt' name will be used.
    Output:
    The function does not return any value. Creates new file in working directory.
    """
    if output_fasta is None:
        output_fasta = 'sorted_blast_result.txt' 
    with open(input_file) as blast_txt:
        output = []
        blast = blast_txt.readlines()
        for i in range(len(blast)):
            if 'Description' in blast[i]:
                output.append(blast[i+1]) 
    sort_output = sorted(output)
    with open (output_file, mode = 'w') as output:
        for line in sort_output:
            output.write(line + '\n')
