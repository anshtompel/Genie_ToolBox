import os
from typing import TextIO, Optional, Union
from modules.fastq_tool import check_gc, check_lengths, check_quality, int_to_tuple, fastq_to_dict, dict_to_fastq
from modules.dna_rna_tools import reverse, complement, reverse_complement, transcribe
from modules.protein_tools import count_protein_mass, count_trypsin_sites, count_seq_length, classify_aminoacids, check_unusual_aminoacids, count_charge,count_aliphatic_index

COMPLEMENT_DNA = {'a': 't', 'A': 'T', 't': 'a', 'T': 'A',
                  'g': 'c', 'G': 'C', 'c': 'g', 'C': 'G'}

COMPLEMENT_RNA = {'a': 'u', 'A': 'U', 'u': 'a', 'U': 'A',
                  'g': 'c', 'G': 'C', 'c': 'g', 'C': 'G'}

DNA_RNA_OPERATIONS = {
                'reverse': reverse, 
                'complement': complement, 
                'reverse_complement': reverse_complement, 
                'transcribe': transcribe
}

AA_GROUPS = {'Nonpolar': ['G', 'A', 'V', 'I', 'L', 'P'],
             'Polar uncharged': ['S', 'T', 'C', 'M', 'N', 'Q'],
             'Aromatic': ['F', 'W', 'Y'],
             'Polar with negative charge': ['D', 'E'],
             'Polar with positive charge': ['K', 'R', 'H']
                }

AMINOACIDS = {'G', 'A', 'V', 'I', 'L', 'P',
              'S', 'T', 'C', 'M', 'N', 'Q', 
              'F', 'W', 'Y', 'D', 'E', 'K', 
              'R', 'H'
                     }

AMINO_ACIDS_MASSES = {
    'G': 57.05, 'A': 71.08, 'S': 87.08, 'P': 97.12, 'V': 99.13,
    'T': 101.1, 'C': 103.1, 'L': 113.2, 'I': 113.2, 'N': 114.1,
    'D': 115.1, 'Q': 128.1, 'K': 128.2, 'E': 129.1, 'M': 131.2,
    'H': 137.1, 'F': 147.2, 'R': 156.2, 'Y': 163.2, 'W': 186.2    
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
    if not input_path.endswith('.fastq'):
        raise ValueError('Incorrect input file extension, should be .fastq')           
    if output_filename is None:
        output_filename = os.path.basename(input_path)        
    gc_params = int_to_tuple(gc_bounds)
    len_bound_params = int_to_tuple(length_bounds)       
    fastq_dictionary = fastq_to_dict(input_path)
    filtered_fastq = {}   
    for read in fastq_dictionary:
        read_sequence = fastq_dictionary[read][0]
        read_quality = fastq_dictionary[read][2]
        if check_gc(read_sequence, gc_params) and check_length(read_sequence, len_bound_params) and check_quality(read_quality, quality_threshold):
            filtered_fastq[read] = fastq_dictionary[read][:]        
    return dict_to_fastq(filtered_fastq, output_filename)
