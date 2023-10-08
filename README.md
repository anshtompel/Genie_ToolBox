# Genie_ToolBox
> This README describes multitool for analysis DNA, peptide and FASTA sequences.

...you have three wishes and here they are!
Genie ToolBox is an instrument for analysis and manipulation with DNA, protein and FASTA sequences.

## Content

1. [Installation](#installation)
2. [Description of operations](#description)
    * [DNA_RNA_TOOLS](#dnarnatools)
    * [PROTEIN_TOOLS](#proteintools)
    * [FASTA_TOOL](#fastatool)
3. [Usage](#usage)
    * [for `dna_rna_tools`](#dnarnatools_usage)
    * [for `protein_tools`](#proteintools_usage)
    * [for `fasta_tool`](#fastatool_usage)
4. [Contancs](#contacts)

## 1.Installation <a name="installation"></a>

1. Clone repo to your local machine using SSH:
``` bash
git@github.com:anshtompel/Genie_ToolBox.git
```
or HTTPS:
``` bash
https://github.com/anshtompel/Genie_ToolBox.git
```

2. Launch `genie_tool_box.py` script and call one of three available functions:
   - `run_dna_rna_tools` - performs transfornations with DNA or RNA sequences according to the complementary base pairing rule;
   -  `run_protein_tools` - counts different phisical and biological properties of protein/peptide;
   -  `run_fasta_filter` - filters FASTA read sequences by input parameters.

## 2.Description of operations <a name="description"></a>

### DNA_RNA_TOOLS <a name="dnarnatools"></a>
`dna_rna_tools` procceds DNA and RNA sequences according to the complementary base pairing rule.
    Function accepts one or arbitary arguments - nucleic acid sequences,
    and returns string if one argument was entered or list for two and more arguments.
    The last argument is always operation name.

`dna_rna_tools` processes your DNA or RNA sequnece using operations:
- `reverse` - reverse your nucleic acid sequence;
- `complement` - get complementary sequence;
- `reverse_complement` - reverse and make a complement to your sequence
- `transcribe` - return transcibed sequence of DNA coding strand

1. `dna_rna_tools` contains `run_dna_rna_tools` function performed transformations with nucleic acid sequences. `run_dna_rna_tools` function takes an arbitrary number of arguments as input, DNA or RNA sequnces (*str*), as well as one of the valid operations (*str*): **reverse**, **transcribe**, **complement** or **reverse_complement**. The function converts sequences and return `str` if one sequence is received as input, or return `list` if more sequences are received.

2. `dna_rna_tools` preserves **the original case** of the letters of input sequences (e.g. **complement AtGc** = **TaCg**)

3. Functions of `dna_rna_tools` work both with DNA and RNA molecules, except `transcribe` operation, which proceeds only DNA sequences to convert them into RNA.  If `transcribe` function recieved RNA and DNA, `list` of DNA transcibed sequences will be returned. If `transcribe` function recieved only RNA, empty `list` will be returned.

4. `dna_rna_tools` checks the correctness of your input sequences and distinguish RNA and DNA. If input sequence contains both **T** and **U** nucleotides, this sequence will be skipped and won't processed.

### PROTEIN_TOOLS <a name="proteintools"></a>

`Protein_tools` is a tool for basic analysis of protein and polypeptide sequenses. Using this tool you can estimate sequence length, charge, aminoacid compound and mass of the protein, find out the aliphatic index and see if the protein could be cleaved by trypsin. Module includes main function - `run_protein_tools` which take the protein sequences and the name of the procedure that the user gives and applies this procedure by one of the available functions to all the given sequences.

Valid operations:
Protein_tools include several operations:
- `count_seq_length`: returns length of protein (int);
- `classify_aminoacids`: returns collection of classified aminoacids, included in the protein (dict);
- `check_unusual_aminoacids`: informs about whether the unusual aminoacis include into the protein (str);
- `count_charge`: returns charge value of protein (int);
- `count_protein_mass`: calculates mass of all aminoacids of input peptide in g/mol scale (float);
- `count_aliphatic_index`: calculates relative proportion of aliphatic aminoacids in input peptide (float);
- `count_trypsin_sites`: counts number of valid trypsin cleavable sites.

Output is a dictionary: key - input sequence, value - output from the chosen procedure as a value.
Also this function check the availabilaty of the procedure and raise the ValueError when the procedure is not in the list of available functions.
**Protein_tools** has several limitations that can raise the errors in the work of the program. Here are some of them:
1. **Protein_Tools** works only with protein sequences that contains letters of Latin alphabet (the case is not important); also every aminoacid should be coded by one letter. If there are other symbols in the sequence, the tool raise `ValueError` *"One of these sequences is not protein sequence or does not match the rools of input. Please select another sequence."*. In this case you should check if there are punctuation marks, spaces or some other symbols in your sequence.
2. Be careful to work only with the sequences that contain aminoacids that coded with one letter. If your sequense is "SerMetAlaGly", **Protein_tools** reads it as "SERMETALAGLY".
3. The list of available functions is available in section "Options". If you see `ValueError` *"This procedure is not available. Please choose another procedure."*, probably your spelling of the name of function is incorrect. Please check the name of chosen prosedure and make sure that it is available in the **Protein_Tools**.

### FASTA_TOOL <a name="fastatool"></a>
Performs filter of input FASTA file according to input parameters.
Input will be filtered by:
- GC content (`gc_bounds`);
- eads length (`length_bounds`);
- reads quality score (`quality_threshold`).

Input:
- `seqs` (dict): FASTA file in dictionary format: key - read ID, value - tuple of sequence and quality.
- `gc_bounds` (tuple or int): GC content filter parameters, it accepts lower and upper (tuple), or only upper threshold value (int). Default value (0, 100).
- `length_bounds` (tuple or int): read length filter parameters, it accepts lower and upper (tuple), or only upper threshold value (int). Default value (0, 2**32).
- `quality_threshold` (int): upper quality threshold in phred33 scale. Reads with average quality below the threshold are discarded. Default value - 0.

Output:
Returns FASTA in dictionary format, as input, only with filtered values which satisfied all input/default conditions.
Dictionary has a format like: key - read ID, value - tuple of sequence and quality.

## 3. Usage <a name="usage"></a>
#### For `dna_rna_tools`:

Examples of runs with one input sequence:
```python
run_dna_rna_tools('GatCa', 'reverse') # 'aCtaG'
run_dna_rna_tools('GatCa', 'complement') # 'CtaGt'
run_dna_rna_tools('GatCa', 'reverse_complement') # 'tGatC'
run_dna_rna_tools('GatCa', 'transcribe') # 'GauCa'
```
or with a few input sequences:

```python
run_dna_rna_tools('GatCa','gtccT', 'AaTcT', 'reverse') # ['aCtaG', 'Tcctg', 'TcTaA']
```

#### For `protein_tools`:

Examples of runs with one input sequence:
```python
sequence = 'CVWGWAMGEACPNPIKINISAYAKTWYQNG'
run_protein_tools(sequence, operation='count_seq_length') # {'CVWGWAMGEACPNPIKINISAYAKTWYQNG': 30}
run_protein_tools(sequence, operation='count_protein_mass') # {'CVWGWAMGEACPNPIKINISAYAKTWYQNG': 3354.9199999999987}
run_protein_tools(sequence, operation='count_trypsin_sites') # {'CVWGWAMGEACPNPIKINISAYAKTWYQNG': 2}
```

or with arbitary number of arguments:

```python
sequence1 = 'CVWGWAMGEAC'
sequence2 = 'PNPIKINISAYAKTWYQNGPIGRCC'
sequence3 = 'CWVGYTAIRFPHQEMQQNTRFNKP'

run_protein_tools(sequence1, sequence2, sequence3, operation='count_seq_length') # {'CVWGWAMGEAC': 11, 'PNPIKINISAYAKTWYQNGPIGRCC': 25, 'CWVGYTAIRFPHQEMQQNTRFNKP': 24}
```

#### For `fasta_tool`:

```python
run_fasta_filter(EXAMPLE_FASTQ) # with default paraments
run_fasta_filter(EXAMPLE_FASTQ, gc_bounds=50) # with only one paramenter, int as upper border
run_fasta_filter(EXAMPLE_FASTQ, gc_bounds=50, length_bounds=100, quality_threshold=30) # with all input paramenters, int as upper border
run_fasta_filter(EXAMPLE_FASTQ, gc_bounds=(30, 40), length_bounds=(60, 1000), quality_threshold=30) # with only one paramenter, tuple - lower and upper border
```

## 4. Contacts <a name="contacts"></a>

- [Anastasia Shtompel](https://github.com/anshtompel) (Telegram: @Aenye): options `count_protein_mass`, `count_aliphatic_index`, `count_trypsin_sites`, modules: `dna_rna_tools`, `fasta_tools`.
- [Elizaveta Chevokina](https://github.com/e-chevokina) (Telegram: @lzchv): options `count_seq_length`, `classify_aminoacids`, `check_unusual_aminoacids`, `count_charge`.

