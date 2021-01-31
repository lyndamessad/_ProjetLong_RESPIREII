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

    parser.add_argument('-i', '-sequence_file', dest='sequence_file', type=str, required=True,
                        help=" Required at least one fasta file \
                        use (all) to use all files of the repository : not fasta files will be removed")

    parser.add_argument('-d', '-database', dest='database', type=str, default = "Uniclst30",
                        help="database for template search"  "available databases: Uniclst30: default, pdb70")    
    
    parser.add_argument('-b', '-bash_command', dest='bash_command', type=bool,
                        default=False, help="Set for True to generate file with bash exécuter commands")

    parser.add_argument('-mt', '-tool', dest='tool', type=str,
                        default="I-Tasser + Rosetta", help="modelling tool. Can choose I-Tasser or Rosetta.\
                        It uses both by default")
    
    parser.add_argument('-m', '-ITasser_mode', dest='ITasser_mode', type=bool, default = False,
                        help="Lunch I-Tasser with options(True) ou not. Defaul : Faslse")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="out.pdb".format(), help="modelling result's file")

    parser.add_argument('-r', '-root_dir', dest='root_dir', type=str, default = os.getcwd(),
                        help="Local directory where pdb folder will be created")

    #file_path = sys.argv[]

    return parser.parse_args()

def files_number(seq):
    """
    this function search for number of input files / all directory
    """
    files = []
    if seq.upper()== "ALL":
        files = os.listdir()
        for file in files:
            extend = str(file.split(".")[1])
            if extend.upper() != "FASTA" or ".gz":
                os.system("rm {}".format(file))
            else:
                pass
        files = os.listdir()
        print(files)

    else:
        files.append(seq)
        print(files)

    return files


def read_fasta(sequence_file ):
    """
    This function is used to read the sequence file.
    Parametre:
    ----------
    sequence_file: gzip fasta file or fasta file// contains sequences
    Return:
    -------
    sequence yield generator.
    """

    for file in sequence_file:
        """    
        #for gziped files:

        if file.endswith(".gz"):
            with gzip.open(file, "rt") as f:
                seqDict = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))
                ident = ident.split("|")[1]
                return seqDict
    
        # for no gziped fasta files:
        else:
        """
        seqRecord = SeqIO.read(file, "fasta")
        sequence = seqRecord.seq
        ident = seqRecord.id
        ident = ident.split("|")[1]
    return ident, sequence

args = get_arguments()
#ident(sequence name), sequence
#ident, seq_str = read_fasta(args.sequence_file)
# Itasser mode 
#Itasser_run_mode(args.ITasser_mode)
# Database choice
#database_choice(args.database)
files_input = files_number(args.sequence_file)
ident, seq_str = read_fasta(files_input)
print(ident)
print(seq_str)
os.system("whereis conda")
