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

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path

def get_arguments():
    """
    Retrieves the arguments of the program.
    Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))

    parser.add_argument('-i', '-sequence_file', dest='sequence_file', type=isfile, required=True,
                        help=" Required at least one fasta file \
                        use (*) to use all files of the repository : not fasta files will be removed")

    parser.add_argument('-b', '-database', dest='database', type=str, default = "Uniclst30",
                        help="database for template search"  "available databases: Uniclst30: default, pdb70")    
    
    parser.add_argument('-out', '-bash_command', dest='bash_command', type=bool,
                        default=False, help="Set for True to generate file with bash exécuter commands")

    parser.add_argument('-mt', '-tool', dest='tool', type=str,
                        default="both", help="modelling tool. Can choose I-Tasser or Rosetta.\
                        It uses both by default")
    
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="out.pdb".format(), help="modelling result's file")

    parser.add_argument('-r', '-root_dir', dest='root_dir', type=str, default = os.getcwd(),
                        help="Local directory where pdb folder will be created")
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
            ident = ident.split("|")[1]
            return seqDict

    # for no gziped fasta files:
    else:
        seqRecord = SeqIO.read(sequence_file, "fasta")
        sequence = seqRecord.seq
        ident = seqRecord.id
        ident = ident.split("|")[1]
        return ident, sequence

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

