import os
from typing import TextIO

def count_gc_content(seq: str) -> float:
    """
    Counts GC content of input sequence reads.
    
    Input:
    - seq (str): input sequence reads from FASTQ file.
    
    Output:
    The GC bases percentage (float) of input sequence (range 0-100).
    """
    gc_content = ((seq.count('G') + seq.count('C')) / len(seq)) * 100
    return gc_content


def filter_gc_bound(fasta_input: dict, gc_params: tuple) -> dict:
    """
    Filters sequences in FASTQ file by GC percentage. 
    
    Input:
    - fasta_input (dict): FASTQ file in dictionary format: key - read ID, value - tuple of sequence and quality.
    - gc_params (tuple): range for filtration, accepts upper or upper/lower border. Both borders are included.
    
    Output:
    Returns input dictionary only with filtered values.
    """
    filtered = {}
    for seq in fasta_input:
        gc_result = count_gc_content(fasta_input[seq][0])
        if gc_params[0] < gc_result < gc_params[1]:
            filtered[seq] = fasta_input[seq][:]
        continue
    return filtered


def filter_length_bounds(result_gc_cound: dict, len_bound_params: tuple) -> dict:
    """
    Filters sequences in FASTQ file by sequence length.
    
    Input:
    - result_gc_cound (dict): FASTA file in dictionary format after filter_gc_bound proccesing: 
    key - read ID, value - tuple of sequence and quality.
    - len_bound_params (tuple): range for filtration, accepts upper or upper/lower border. Both borders are included.
    
    Output:
    Returns input dictionary only with filtered values.
    """
    filtered = {}
    for seq in result_gc_cound:
        len_of_seq = len(result_gc_cound[seq][0])
        if len_bound_params[0] < len_of_seq < len_bound_params[1]:
            filtered[seq] = result_gc_cound[seq][:]
        continue
    return filtered


def count_quality_score(seq: str) -> float:
    """
    Counts mean quality score of input sequence based on the value of each symbol in ASCII table.
    
    Input:
    - seq (str): quality score from FASTQ file. 
    
    Output:
    Mean quality score of input sequence (float).
    
    """
    score = 0
    for symbol in seq:
        score += ord(symbol) - 33
    mean_score = score / len(seq)
    return mean_score     
 

def filter_quality_threshold(result_len_bound: dict, quality_params: int) -> dict:
    """
    Filters sequences in FASTQ file by mean quality score.
    
    Input:
    - result_len_bound (dict): FASTQ file in dictionary format after filter_length_bounds proccesing: 
    key - read ID, value - tuple of sequence and quality.
    - quality_params (int): threshold value of reads quality in phred33 scale.
    
    Output:
    Returns input dictionary only with filtered values.
    """
    filtered = {}
    for seq in result_len_bound:
        ascii_quality = result_len_bound[seq][2]
        if count_quality_score(ascii_quality) >= quality_params:
            filtered[seq] = result_len_bound[seq][:]
        continue
    return filtered


def int_to_tuple(input_parameters) -> tuple:
    """
    Converts input parameters to tuple format.
    If input is already a tuple, it will be return without changes.
    If input parameter is int (only upper threshold is entered), function will return a tuple like (0, 'input').
    
    Input:
    - input_parameters (tuple or int).
    
    Output:
    Lower and upper threshold for functions in tuple format.
    """
    if isinstance(input_parameters, tuple):
        return input_parameters
    return (0, input_parameters)


def fastq_to_dict(file: TextIO) -> dict:
    """
    Converts input FASTQ file to dictionary. Dictionary has four strings: (1) - read ID (str), 
    (2) sequence, commentary and quality (tuple of str). Read ID is identified as started with '@' string and include
    one of this strings in the end: 'BH:ok', 'BH:failed' or 'BH:changed'. 
    
    Input:
    - .fastq file (textIO format): FASTQ file;
    
    Output:
    Dict(key - str; value - tuple(str)) - dictionary from input FASTQ file.
    
    """
    with open (file) as fastq_file:
        seq_com_quality = ()
        fastq_dict = {}
        for line in fastq_file:
            line = line.strip()
            if line.startswith('@') and ('BH:ok' in line or 'BH:failed' in line or 'BH:changed' in line):
                read_id = line
            else:
                seq_com_quality += tuple([line])
                if len(seq_com_quality) == 3:
                    fastq_dict[read_id] = seq_com_quality
                    seq_com_quality = ()
                continue
        return fastq_dict

    
def dict_to_fastq(fastq_dict: dict, output_filename: str) -> TextIO: 
    """
    Converts filtetred dictionary to FASTQ file. Dictionary has four strings: (1) - read ID (str), 
    (2) sequence, commentary and quality (tuple of str). Every string is written sequentially. 
    
    Input:
    - Dict(key - str; value - tuple(str)) - dictionary of filtered reads;
    
    Output:
    .fastq file (textIO format): FASTQ file with filtered reads.
    """
    path = 'fastq_filtrator_resuls'
    if not os.path.exists(path):
        os.makedirs(path)
    with open(os.path.join(path, output_filename), mode='w') as output:
        for read_id in fastq_dict:
            output.write(read_id + '\n')
            for seq_com_quality in fastq_dict[read_id]:
                output.write(seq_com_quality + '\n')
    return output


def run_fastq_filter(input_path = None, output_filename = None, gc_bounds = (0, 100), length_bounds = (0, 2**32), quality_threshold = 0) -> TextIO:
    """
    Performs filter of input FASTQ file according to input parameters. 
    Input will be filtered by: 
        - GC content (gc_bounds);
        - reads length (length_bounds);
        - reads quality score (quality_threshold).
        
    Input:
    - input_path (str): path to .fastq file; include 4 strings: 1 - read ID, 2 - sequence, 3 - comment, 4 - quality. Default - None.
    - output_filename (str): name of output file, by default, it will be saved in the directory 'fastq_filtrator_resuls'. Default name will be name of input file.
    - gc_bounds (tuple or int): GC content filter parameters, it accepts lower and upper (tuple), or only upper threshold value (int). Default value (0, 100).
    - length_bounds (tuple or int): read length filter parameters, it accepts lower and upper (tuple), or only upper threshold value (int). Default value (0, 2**32).
    - quality_threshold (int): upper quality threshold in phred33 scale. Reads with average quality below the threshold are discarded. Default value - 0. 
    
    Output:
    Returns FASTQ only with filtered reads which satisfied all input/default conditions.
    """
    if not input_path.endswith('.fastq'):
        raise ValueError('Incorrect input file extension, should be .fastq')
    
    fastq_dictionary = fastq_to_dict(input_path)
    gc_params = int_to_tuple(gc_bounds)
    len_bound_params = int_to_tuple(length_bounds) 
    
    result_gc_bound = filter_gc_bound(fastq_dictionary, gc_params)
    result_len_bound = filter_length_bounds(result_gc_bound, len_bound_params)
    filtered_result = filter_quality_threshold(result_len_bound, quality_threshold)
    
    if output_filename is None:
        output_filename = os.path.basename(input_path)
        
    return dict_to_fastq(filtered_result, output_filename)