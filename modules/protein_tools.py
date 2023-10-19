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


def is_protein(peptide: str) -> bool:
    """
    'is_protein' function check if the sequence contains only letters in the upper case.
    Input: a protein sequence (a str type).
    Output: boolean value.
    """
    return peptide.isalpha() and peptide.isupper()


def count_seq_length(peptide: str) -> int:
    """
    'count_seq_length' function counts the length of protein sequence.
    Input: a protein sequence (a str type).
    Output: length of protein sequence (an int type).
    """
    return len(peptide)


def classify_aminoacids(peptide: str) -> dict:
    """
    classify_aminoacids' function classify all aminoacids from the input sequence 
    in accordance with the 'AA_ALPHABET' classification. If aminoacid is not included in this list,
    it should be classified as 'Unusual'.
    Input: a protein sequence (a str type).
    Output: a classification of all aminoacids from the sequence (a dict type — 'all_aminoacids_classes' variable).
    """
    all_aminoacids_classes = dict.fromkeys(['Nonpolar', 'Polar uncharged', 'Aromatic', 'Polar with negative charge', 'Polar with positive charge', 'Unusual'], 0)
    for aminoacid in peptide:
        aminoacid = aminoacid.upper()
        if aminoacid not in AMINOACIDS:
            all_aminoacids_classes['Unusual'] += 1
        for aa_key, aa_value in AA_GROUPS.items():
            if aminoacid in aa_value:
                all_aminoacids_classes[aa_key] += 1
    return all_aminoacids_classes


def check_unusual_aminoacids(peptide: str) -> str:
    """
    'check_unusual_aminoacids' function checks the composition of aminoacids and return the list of unusual aminoacids if they present in the sequence. We call the aminoacid
    unusual when it does not belong to the list of proteinogenic aminoacids (see 'ALL_AMINOACIDS' global variable).
    Input: a protein sequence (a str type).
    Output: an answer whether the sequense contains unusual aminoacids (a str type).
    """
    seq_aminoacids = set()
    for aminoacid in peptide:
        seq_aminoacids.add(aminoacid)
    if seq_aminoacids.issubset(AMINOACIDS):
        return 'This sequence contains only proteinogenic aminoacids.'
    else:
        unusual_aminoacids = seq_aminoacids - AMINOACIDS
        unusual_aminoacids_str = ''
        for elem in unusual_aminoacids:
            unusual_aminoacids_str += elem
            unusual_aminoacids_str += ', '
        return f'This protein contains unusual aminoacids: {unusual_aminoacids_str[:-2]}.'


def count_charge(peptide: str) -> int:
    """
    'count_charge' function counts the charge of the protein by the subtraction between the number of positively and negatively charged aminoacids.
    Input: a protein sequence (a str type).
    Output: a charge of the sequence (an int type).
    """
    seq_classes = classify_aminoacids(peptide)
    positive_charge = seq_classes['Polar with positive charge']
    negative_charge = seq_classes['Polar with negative charge']
    sum_charge = positive_charge - negative_charge
    return sum_charge


def count_protein_mass(peptide: str) -> float:
    """
    Calculates mass of all aminoacids of input peptide in g/mol scale.
    Arguments:
    - seq (str): one-letter code peptide sequence, case is not important;
    Output:
    Returns mass of peptide (float).
    """
    aa_mass = 0
    for aminoacid in peptide:
        aa_mass += AMINO_ACIDS_MASSES[aminoacid]
    return aa_mass


def count_aliphatic_index(peptide: str) -> float:
    """
    Calculates aliphatic index - relative proportion of aliphatic aminoacids in input peptide.
    The higher aliphatic index the higher thermostability of peptide.
    Argument:
    - seq (str): one-letter code peptide sequence, letter case is not important.
    Output:
    Returns alipatic index (float).
    """
    aliphatic_aa_fraction = 0
    aliphatic_aa_fraction += peptide.count('A') + 3.9 * peptide.count('L') + 3.9 * peptide.count('I') +  2.9 * peptide.count('V')
    aliph_index = (aliphatic_aa_fraction/count_seq_length(peptide))
    return aliph_index


def count_trypsin_non_cleaved(peptide: str) -> int:
    """
    Counts non-cleavable sites of trypsin: Arginine/Proline (RP) and Lysine/Proline (KP) pairs.
    Argument:
    - seq (str): one-letter code peptide sequence, case is not important.
    Output:
    Returns number of exception sites that cannot be cleaved by trypsin (int).
    """
    not_cleavage_count = 0
    not_cleavage_count += peptide.count('RP') + peptide.count('KP')
    return not_cleavage_count


def count_trypsin_sites(peptide: str) -> int:
    """
    Counts number of valid trypsin cleavable sites:
    Arginine/any aminoacid and Lysine/any aminoacid (except Proline).
    Argument:
    - seq (str): one-letter code peptide sequence, case is not important.
    Output:
    Returns number of valid trypsin cleavable sites (int).
    If peptide has not any trypsin cleavable sites, it will return zero.
    """
    arginine_lysine_count = peptide.count('R') + peptide.count('K')
    cleavage_sites_count = arginine_lysine_count - count_trypsin_non_cleaved(peptide)
    return cleavage_sites_count
  

OPERATIONS = {'count_protein_mass':count_protein_mass,
             'count_aliphatic_index': count_aliphatic_index,
             'count_trypsin_sites': count_trypsin_sites,
             'count_seq_length': count_seq_length,
             'classify_aminoacids': classify_aminoacids,
             'check_unusual_aminoacids': check_unusual_aminoacids,
             'count_charge': count_charge}


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
