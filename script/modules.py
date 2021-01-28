import getpass
import argparse
import logging
import os
import sys
import os.path as op
import shutil
import subprocess
import gzip
from Bio import SeqIO
from Bio.Seq import Seq


def cast_to_str(obj):
    """Return a string representation of a Seq or SeqRecord.
    Args:
        obj (str, Seq, SeqRecord): Biopython Seq or SeqRecord
    Returns:
        str: String representation of the sequence
    """

    if isinstance(obj, str):
        return obj
    if isinstance(obj, Seq):
        return str(obj)
    if isinstance(obj, SeqRecord):
        return str(obj.seq)
    else:
        raise ValueError('Must provide a string, Seq, or SeqRecord object.')


def read_fasta(sequence_file :str):
    """
    This function is used to read the sequence file.
    Parametre:
    ----------
    sequence_file: gzip fasta file or fasta file// contains sequences
    Return:
    -------
    sequence yield generator.
    """

    #for gziped files:

    if sequence_file.endswith(".gz"):
        with gzip.open(sequence_file, "rt") as file:
            seqDict = SeqIO.to_dict(SeqIO.parse(file, 'fasta'))
            return cast_to_str(seqDict)

    # for no gziped fasta files:
    else:
        seqRecord = SeqIO.read(sequence_file, "fasta")
        sequence = cast_to_str(seqRecord.seq)
        ident = seqRecord.id
        ident = ident.split("|")[1]
        return ident, sequence

def get_path(file):
    path_file = os.path.split(file)[0] +"/"
    if path_file == "/":
        path_file = " "
    return path_file

def align_script(ident, seq_file, output_file):
    """Take in a sequence string and prepares the folder for the HHblits run."""
    path_file = get_path(seq_file)
    print(path_file)
    with open(output_file,"w") as f:

        f.write("#run alignement on the database to found templats\n")
        f.write("hhblits -cpu 30 -i {}. -d ~/Desktop/perso_lynda/databases/uniclust30/UniRef30_2020_06_hhsuite/UniRef30_2020_06 -oa3m ".format(seq_file))
        f.write(path_file)
        f.write("{}.a3m\n\n".format(ident))
        # -M : score POUR LA CONTRAINTE DU FORMAT FASTA 

        f.write("#to add secondary structure information to the MSA\n")
        f.write("addss.pl ")
        f.write(path_file)
        f.write("{}.a3m\n\n".format(ident))
        
        f.write("#generate a hidden Markov model (HMM) from this MSA\n")
        f.write("hhmake -i ")
        f.write(path_file)
        f.write("{}.a3m\n\n".format(ident))
        
        f.write("#to convert on fasta file\n")
        f.write("reformat.pl a3m fas {}{}.a3m ".format(path_file,ident))
        f.write("{}{}_blits.fas\n".format(path_file,ident))
    # Deplacer le fichier dans le repertoire HHblits_script + le renomer 
    os.rename(output_file, ident+"_HHblits.sh")
    os.system("chmod a+xwr *.sh")

if __name__ == '__main__':
    pass  
