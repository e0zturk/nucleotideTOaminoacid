from pathlib import Path
import sys

data = sys.argv[1]

def processed (dt):
    path = Path(dt)
    if path.suffix in [".txt", ".fasta", ".fa"]:
        with open(path,"r") as S:
           sq = S.readlines()[1:]
           seq = ""
           for s in sq:
               seq += s.strip().upper()
        return seq
    else:
        return dt

input_seq = processed(data).upper()

def controlled (preseq):
    if set(preseq).issubset({"A","G","T","U","C"}):
        return preseq
    else:
        print(f"Invalid characters found in your sequence\n"
              f"please check  '{preseq}'")

prechecked = controlled(input_seq)

if prechecked is None:
    print("Sequence validation failed.")

    print("Please provide a valid DNA/RNA sequence (only A, T, G, C, U allowed).")
    sys.exit(1)
else:
    dnasequence = prechecked.replace("U", "T")

referance = Path("Codon-Amino Acid Abbreviations.tsv")

with open(referance, "r") as Abb:
    abbreviations = Abb.readlines()
    ref = []
    mrnacodon = []
    aminoacid = []
    for l in abbreviations[1:]:
        predata = l.strip().split("\t")
        Codon = predata[0]
        Fullname = predata[1]
        Abb3letter = predata[2]
        Abb1letter = predata[3]
        ref.append([Codon,Fullname,Abb3letter,Abb1letter])
        for i in ref:
            mrnacodon.append(i[0])
            aminoacid.append(i[3])

substitutions = dict(zip(mrnacodon, aminoacid))

def translation (subs,sequence):
    final_seq = ""
    for t in range(0,len(sequence),3):
        triplet = sequence[t:t+3]
        final_seq += subs.get(triplet, "X")
    return final_seq


translasyon = translation(substitutions,dnasequence)

if len(sys.argv) < 3 :
    print(translation(substitutions,dnasequence))

if len(sys.argv) > 2 and sys.argv[2] == "fasta":
     with open("uncategorized protein.fasta", "w") as W:
        W.write("> uncategorized protein sequance\n")
        for p in range(0,len(translasyon),80):
            pline = translasyon[p:p+80]
            W.write(f"{pline}\n")
            print("The fasta file has been successfully saved.")

elif len(sys.argv) > 2 and sys.argv[2] == "txt":
    with open("uncategorized protein.txt", "w") as W:
        W.write("> uncategorized protein sequance\n")
        for p in range(0, len(translasyon), 80):
            pline = translasyon[p:p + 80]
            W.write(f"{pline}\n")
            print("The text file has been successfully saved.")