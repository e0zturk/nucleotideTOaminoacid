# nucleotideTOaminoacidSBT
Have you got an own nucleotide sequence? don't worry! By means of this repository, your nucleotide/dna sequence will be translated into an amino acid sequence. You can choose whichever of chosen format for output (txt,fasta)...

# requirements
 1. Python 3.8+
 2. codonabbreviations.tsv file must be located in the directory.

# How to work?
1. Download the files to your local files.

2. Note the definitions:
   __input_file__ : A file in text or fasta format containing your DNA/nucleotide sequences
   __manual_DNA_sequence__ : A DNA/nucleotide sequence that can be entered manually without using a file
   __output_format__ : The output file formats you obtain after translation. Two formats are supported: 'fasta', 'txt'
   
3. Navigate to the project folder in Terminal or PowerShell and run the following command: (It can be operated in more different combinations.)

   __A.__
     
   python input_file output_format

   _for example:_

   python example.txt txt

   __B.__
   
   python manuel_DNA_sequance

   _for example:_

   python ATGCGTAGCCGATACGGATCA

   _output_
   MRSRYGS

   __C.__
   
   python manuel_DNA_sequence output_format

   _for example:_

   python ATGCGTAGCCGATACGGATCA fasta


# Possible Errors and Notes
 
1. __Missing input argument__
   
   If you run the script without providing an input file or sequence, Python will raise:
   IndexError: list index out of range
   ->Always provide at least one argument (a sequence string or a file path).

3. __Unsupported file extension__
     
   Only .txt, .fasta, or .fa files are processed as sequence files. Any other extension will be treated as a raw sequence string.
   Invalid characters in sequence
   If the sequence contains characters other than A, T, G, C, U, the program prints:

   Invalid characters found in your sequence
   please check '...'
   Sequence validation failed.
   Please provide a valid DNA/RNA sequence (only A, T, G, C, U allowed).

4. __Missing codon table file__
   
   If codonabbreviations.tsv is not present in the same directory, you will get:

   FileNotFoundError: [Errno 2] No such file or directory: 'codonabbreviations.tsv'

5. __Incomplete codon at the end of sequence__
   
   If the sequence length is not divisible by 3, the last 1–2 nucleotides will be translated as  (unknown amino acid).

6. __Unknown codon__
    
   Any codon not found in the reference table is translated as .

7. __Output file option not specified__
    
   If you don’t provide a second argument ( or ), the translation will only be printed to the console and no file will be saved.

   






   
