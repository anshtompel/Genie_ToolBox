COMPLEMENT_DNA = {'a': 't', 'A': 'T', 't': 'a', 'T': 'A',
                  'g': 'c', 'G': 'C', 'c': 'g', 'C': 'G'}

COMPLEMENT_RNA = {'a': 'u', 'A': 'U', 'u': 'a', 'U': 'A',
                  'g': 'c', 'G': 'C', 'c': 'g', 'C': 'G'}


def check_dna(sequence: str) -> bool:
    """
    Checkes is input string a DNA sequence.
    If input sequence is DNA it will return True.
    
    Input:
    -sequence (str): input nucleic acid sequence, any letter case is suitable.
    
    Output:
    returns True or False (bool), is nucleic acid is DNA.
    """
    flag = True
    valid_nucleotides_dna = {'A', 'T', 'C', 'G'}
    unique_chars = set(sequence.upper())
    if not unique_chars.issubset(valid_nucleotides_dna):
        flag = False
    return flag


def check_rna(sequence: str) -> bool:
    """
    Checkes is input string a RNA sequence.
    
    Input:
    -sequence (str): input nucleic acid sequence, any letter case is suitable.
    
    Output:
    returns True or False (bool), is nucleic acid is RNA.
    """
    flag = True
    valid_nucleotides_rna = {'A', 'U', 'C', 'G'}
    unique_chars = set(sequence.upper())
    if not unique_chars.issubset(valid_nucleotides_rna):
        flag = False
    return flag


def reverse(*sequences: list) -> list:
    """
    Creates reverse sequence(s) of RNA or DNA. Function accepts one or arbitary arguments - nucleic acid sequences,
    and returns string if one argument was entered or list for two and more arguments.
    
    Input:
    - *sequences list[str]: list of DNA or RNA sequences.
    
    Output:
    list[str] or str of reversed sequence(s).
    """
    reverse_sequences = []
    for sequence in sequences:
        if check_dna(sequence) or check_rna(sequence) is True:
            reverse_sequence = sequence[::-1]
            reverse_sequences.append(reverse_sequence)
    if len(reverse_sequences) == 1:
        return reverse_sequences[0]
    return reverse_sequences


def complement_dna(sequence: str) -> str:
    """
    Returns DNA complement sequnence (str).

    Input:
    - sequence (str): input DNA sequnence.
    
    Output:
    DNA complement sequnence (str).
    """
    complement_sequence = ''
    for nucleotide in sequence:
        if nucleotide in COMPLEMENT_DNA:
            complement_sequence += COMPLEMENT_DNA[nucleotide]
    return complement_sequence


def complement_rna(sequence: str) -> str:
    """
    Returns RNA complement sequnence (str).
    
    Input:
    - sequence (str): input RNA sequnence.
    
    Output:
    RNA complement sequnence (str).
    """
    complement_sequence = ''
    for nucleotide in sequence:
        if nucleotide in COMPLEMENT_RNA:
            complement_sequence += COMPLEMENT_RNA[nucleotide]
    return complement_sequence


def complement(*sequences: list) -> list:
    """
    Create complement sequence(s) of RNA or DNA acoording to the complementary base pairing rule. 
    Function accepts one or arbitary arguments - nucleic acid sequences, 
    and returns string if one argument was entered or list for two and more arguments. 
    
    Input:
    - *sequences list[str]: list of DNA or RNA sequences.
    
    Output:
    list[str] or str of complement sequence(s).
    """
    complement_sequences = []
    for sequence in sequences:
        if check_dna(sequence) is True:
            complement_sequences.append(complement_dna(sequence))
        elif check_rna(sequence) is True:
            complement_sequences.append(complement_rna(sequence))
    if len(complement_sequences) == 1:
        return complement_sequences[0]
    return complement_sequences


def reverse_complement(*sequences: list) -> list:
    """
    Create reversed complement sequence(s) of RNA or DNA acoording to the complementary base pairing rule. 
    Function accepts one or arbitary arguments - nucleic acid sequences, 
    and returns string if one argument was entered or list for two and more arguments. 
    
    Input:
    - *sequences list[str]: list of DNA or RNA sequences.
    
    Output:
    list[str] or str of reversed complement sequence(s).
    """
    rev_comp_seqs = []
    for sequence in sequences:
        reverse_seq = reverse(sequence)
        rev_comp_seq = complement(reverse_seq)
        rev_comp_seqs.append(rev_comp_seq)
    if len(rev_comp_seqs) == 1:
        return rev_comp_seqs[0]
    return rev_comp_seqs


def transcribe(*sequences: list) -> list:
    """
    Create transcribed sequence(s) of DNA acoording to the complementary base pairing rule. 
    Function accepts one or arbitary arguments - nucleic acid sequences, 
    and returns string if one argument was entered or list for two and more arguments. 
    Function can procced only DNA seqence(s), if it is received RNA it will return empty list.
    
    Input:
    - *sequences list[str]: list of DNA or RNA sequences.
    
    Output:
    list[str] or str of reversed complement sequence(s).
    """
    transcribed_sequences = []
    for sequence in sequences:
        if check_dna(sequence) is True:
            transcribed_sequence = sequence.replace('t', 'u').replace('T', 'U')
            transcribed_sequences.append(transcribed_sequence)
    if len(transcribed_sequences) == 1:
        return transcribed_sequences[0]
    return transcribed_sequences


DNA_RNA_OPERATIONS = {
                'reverse': reverse,
                'complement': complement,
                'reverse_complement': reverse_complement,
                'transcribe': transcribe
}


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
    return OPERATIONS[operation](*sequences)
