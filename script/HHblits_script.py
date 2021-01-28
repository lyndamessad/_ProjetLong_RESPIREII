from modules import *
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


def get_arguments():
    """
    Retrieves the arguments of the program.
    Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-sequence_file', dest='sequence_file', type=str, required=True,
                        help="Sequence is a fasta file (.fasta or .fasta.gz)")


    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="HHblits_script.sh", help="Output file")

    #file_path = sys.argv[]

    return parser.parse_args()

    


def run_HHblits(ident):
    """
        run HHblits
    """
    args = get_arguments()
    align_file = align_script(ident, args.sequence_file, args.output_file)
    #exécuter le script shell
    print("Run HHblits \n")

    os.system("./{}_HHblits.sh".format(ident))
    
    print("HHblits's results in file {}.fas\n".format(ident))


def found_template(hhblits_result, tmp_file, type):
    # output_file de hhblits: ident_blits.fas

    # type 1 : alignement seq+ tmp 
    # type 2 : struct secondaire
    # tmp_file doit contenir la seq target + best tmp.
    """ found the best template and it's alignement"""
    with open(hhblits_result,'r') as f1:
        with open(tmp_file,'w') as f2:
            pass



def found_ss(hhblits_result, ss_file):
    """
    found secondary structure of the target proteine
    """
    with open(hhblits_result,'r') as f1:
        with open(ss_file, 'w') as f2: 
            pass

def main():
    """
    Main program function
    """
    
    # get arguments
    args = get_arguments()
    #ident(sequence name), sequence
    ident, seq_str = read_fasta(args.sequence_file)

    # Preparation du script de hhblits par la fct align_script de modules
    # LAncement de l'alignement sur HHblits
    hhblits_result = run_HHblits(ident)

if __name__ == '__main__':
    main()

    # Verifier l'installation des logiciels: HHblits, Rosetta et I-Tasser

