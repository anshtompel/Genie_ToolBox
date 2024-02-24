import os
from typing import TextIO, Optional, Union
from modules.dna_rna_tools import reverse, complement, reverse_complement, transcribe
from modules.protein_tools import count_protein_mass, count_trypsin_sites, count_seq_length, classify_aminoacids, check_unusual_aminoacids, count_charge,count_aliphatic_index, is_protein
from Bio import SeqIO, SeqUtils


DNA_RNA_OPERATIONS = {
                'reverse': reverse, 
                'complement': complement, 
                'reverse_complement': reverse_complement, 
                'transcribe': transcribe
}

OPERATIONS = {'count_protein_mass':count_protein_mass,
             'count_aliphatic_index': count_aliphatic_index,
             'count_trypsin_sites': count_trypsin_sites,
             'count_seq_length': count_seq_length,
             'classify_aminoacids': classify_aminoacids,
             'check_unusual_aminoacids': check_unusual_aminoacids,
             'count_charge': count_charge}


def run_dna_rna_tools(*strings: str) -> list:
    """
    Procceds DNA and RNA sequences according to the complementary base pairing rule. 
    Function accepts one or arbitary arguments - nucleic acid sequences, 
    and returns string if one argument was entered or list for two and more arguments. 
    The last argument is always operation name.
    
    Input:
    *strings list[str]: list of DNA or RNA sequences.
    The last argument is one of the followed operations:
    
            - reverse: creates reverse sequence(s) of RNA or DNA.
            - complement: create complement sequence(s) of RNA or DNA. 
            - reverse_complement: Create reversed complement sequence(s) of RNA or DNA. 
            - transcribe: create transcribed sequence(s) of DNA.
            
    Function returns sequences in the same letter case as input.
    If invalid sequence is entered, function returns empty list.
    """
    sequences, operation  = strings[:-1], strings[-1]
    return DNA_RNA_OPERATIONS[operation](*sequences)


def run_protein_tools(*peptides: str, operation = None) -> dict:
    """
    'run_protein_tools' function take the protein sequences and the name of the procedure that the user gives and applies this procedure by one of the available functions
    to all the given sequences.
    Calculates protein phisical properties: mass, charge, length, aliphatic index;
    as well as defines biological features: aminoacid composition, trypsin cleavable sites.
    
    Input: a list of protein sequences and one procedure that should be done with these sequences (str type, several values).

    Valid operations: 
    Protein_tools include several operations:
    - count_seq_length: returns length of protein (int);
    - classify_aminoacids: returns collection of classified aminoacids, included in the protein (dict);
    - check_unusual_aminoacids: informs about whether the unusual aminoacis include into the protein (str);
    - count_charge: returns charge value of protein (int);
    - count_protein_mass: calculates mass of all aminoacids of input peptide in g/mol scale (float);
    - count_aliphatic_index: calculates relative proportion of aliphatic aminoacids in input peptide (float);
    - count_trypsin_sites: counts number of valid trypsin cleavable sites.
    
    Output: a dictionary of input sequence as a key, and outputs from the chosen procedure as a value.
    Also this function check the availabilaty of the procedure and raise the ValueError when the procedure is not in the list of available
    functions (see 'OPERATIONS' global variable).
    """
    operation_result = {}
    for peptide in peptides:
        if not is_protein(peptide):
            raise ValueError("One of these sequences is not protein sequence or does not match the rools of input. Please select another sequence.")
        else:
            if operation in OPERATIONS:
                result = OPERATIONS[operation](peptide)
                operation_result.update({peptide: result})
            else:
                raise ValueError("This procedure is not available. Please choose another procedure.")
    return operation_result


#FASTQ Filtraror with BIOPYTHON utilities

def check_gc(fastq_read: str, gc_params: tuple) -> bool:
    """
    Filters sequences in FASTQ file by GC percentage. 
    
    Input:
    - fastq_read (str): FASTQ sequence read.
    - gc_params (tuple): range for filtration, accepts upper or upper/lower border. Both borders are included.
    
    Output:
    Returns boolean value is this read satisfied the criteria.
    """
    gc_result = SeqUtils.GC123(fastq_read)[0]
    if gc_params[0] < gc_result < gc_params[1]:
        return True


def check_length(fastq_read: str, len_bound_params: tuple) -> bool:
    """
    Filters sequences in FASTQ file by sequence length.
    
    Input:
    - fastq_read (str): FASTQ sequence read.
    - len_bound_params (tuple): range for filtration, accepts upper or upper/lower border. Both borders are included.
    
    Output:
    Returns boolean value is this read satisfied the criteria.
    """
    len_of_seq = len(fastq_read)
    if len_bound_params[0] < len_of_seq < len_bound_params[1]:
        return True


def check_quality(fastq_quality: str, quality_params: int) -> bool:
    """
    Filters sequences in FASTQ file by mean quality score.
    
    Input:
    - fastq_quality (str): FASTQ read quality
    - quality_params (int): threshold value of reads quality in phred33 scale.
    
    Output:
    Returns boolean value is this read satisfied the criteria.
    """ 
    mean_quality = sum(fastq_quality.letter_annotations["phred_quality"])/len(fastq_quality.seq)
    if mean_quality >= quality_params:
        return True


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


def run_fastq_filter(input_path: Optional[str] = None, output_filename: Optional[str] = None, gc_bounds: Union[int, tuple] = (0, 100), length_bounds: Union[int, tuple] = (0, 2**32), quality_threshold: int = 0) -> TextIO:
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

    "Specify input path"
    if not input_path.endswith('.fastq'):
        raise ValueError('Incorrect input file extension, should be .fastq')   

    "Specify output path"
    output_path = 'fastq_filtrator_resuls'
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    if output_filename is None:
        output_filename = os.path.basename(input_path)
    "Passed the parameters"
    gc_params = int_to_tuple(gc_bounds)
    len_bound_params = int_to_tuple(length_bounds)    
    "Filter and record results"
    filtererd_fastq = open(os.path.join(output_path, output_filename), mode='w')
    for seq_record in SeqIO.parse(input_path, "fastq"):
        if check_gc(seq_record.seq, gc_params) and check_length(seq_record.seq, len_bound_params) and check_quality(seq_record, quality_threshold):
                SeqIO.write(seq_record, filtererd_fastq, "fastq")  
    filtererd_fastq.close()   
    return filtererd_fastq
