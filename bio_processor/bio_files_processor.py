import os

AMINOACIDS = {'G', 'A', 'V', 'I', 'L', 'P',
              'S', 'T', 'C', 'M', 'N', 'Q', 
              'F', 'W', 'Y', 'D', 'E', 'K', 
              'R', 'H'}



def dict_to_fasta(fasta_dict, output_fasta=None):   
    if output_fasta is None:
        output_fasta = 'file.fasta' 
    if not output_fasta.endswith('.fasta'):
        output_fasta += '.fasta'
    with open(output_fasta, mode='w') as output:
        for string in fasta_dict.items():
            output.write(string[0] + '\n')
            output.write(string[1] + '\n')    
    return output


def convert_multiline_fasta_to_oneline(input_fasta, output_fasta = None):
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

    
def filter_gbk(file):
    with open(file) as lines:
        seq_list = []
        gene = None
        translation = None
        for line in lines:
            line = line.strip()
            if line.startswith('/gene'):
                gene = line.strip('/gene="')
                gene = '>'+ gene
                seq_list.append(gene)
            if line.startswith('/translation') or set(line).issubset(AMINOACIDS):
                if gene is None:
                    continue
                else:
                    translation = line.strip('/translation="')
                    seq_list.append(translation)
    return seq_list


def get_gene_names(seq_list):
    names = []
    for seq in my_list:
        if seq.startswith('>'):
            names.append(seq.strip('>'))
    return names


def find_range_genes(seq_list, names_seq, gene, n_before, n_after):
    upper_section = []
    lower_section = []
    left = names_seq.index(gene) - n_before   
    right = names_seq.index(gene) + n_after   
    upper = seq_list.index('>'+ names_seq[left])
    selected_gene = seq_list.index('>'+ gene)
    lower = seq_list.index('>'+ names_seq[right])   
    for i in range(upper, selected_gene):
        upper_section.append(my_list[i])
    for i in range(lower, selected_gene, -1):
        lower_section.append(my_list[i])       
    sum_list_genes =  upper_section + lower_section
    return sum_list_genes


def write_to_fasta(gene_list, output_fasta):
    with open(output_fasta, mode='w') as output:
        for line in gene_list:
            output.write(line + '\n')


def select_genes_from_gbk_to_fasta(input_gbk, genes, output_fasta = None, n_before = 1, n_after = 1):
    filter_gene_list = filter_gbk(input_gbk)
    gene_of_interest_list = get_gene_names(filter_gene_list)
    output_list = find_range_genes(filter_gene_list, gene_of_interest_list, genes, n_before, n_after)
    if output_fasta is None:
        output_fasta = 'file.fasta' 
    if not output_fasta.endswith('.fasta'):
        output_fasta += '.fasta'
    write_to_fasta(output_list, output_fasta)
    return f'Fasta file was written in {os.path.join(os. getcwd(), output_fasta)}.'
