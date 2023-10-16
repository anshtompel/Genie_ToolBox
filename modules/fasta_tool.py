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

def run_fasta_filter(seqs: dict, gc_bounds = (0, 100), length_bounds = (0, 2**32), quality_threshold = 0) -> dict:
    """
    Performs filter of input FASTA file according to input parameters. 
    Input will be filtered by: 
        - GC content (gc_bounds);
        - reads length (length_bounds);
        - reads quality score (quality_threshold).
        
    Input:
    - seqs (dict): FASTA file in dictionary format: key - read ID, value - tuple of sequence and quality.
    - gc_bounds (tuple or int): GC content filter parameters, it accepts lower and upper (tuple), or only upper threshold value (int). Default value (0, 100).
    - length_bounds (tuple or int): read length filter parameters, it accepts lower and upper (tuple), or only upper threshold value (int). Default value (0, 2**32).
    - quality_threshold (int): upper quality threshold in phred33 scale. Reads with average quality below the threshold are discarded. Default value - 0. 
    
    Output:
    Returns FASTA in dictionary format, as input, only with filtered values which satisfied all input/default conditions.
    Dictionary has a format like: key - read ID, value - tuple of sequence and quality.                           
    """
    gc_params = int_to_tuple(gc_bounds)
    len_bound_params = int_to_tuple(length_bounds) 
    result_gc_bound = filter_gc_bound(seqs, gc_params)
    result_len_bound = filter_length_bounds(result_gc_bound, len_bound_params)
    filtered_result = filter_quality_threshold(result_len_bound, quality_threshold)
    return filtered_result
