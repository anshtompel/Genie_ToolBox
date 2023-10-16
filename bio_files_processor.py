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
