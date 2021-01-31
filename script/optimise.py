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
            if extend.upper() != "FASTA":
                os.system("rm {}".format(file))
            else:
                pass
        files = os.listdir()

    else:
        files.append(seq)

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

def Itasser_run_mode(mode:bool):
    # Itasser mode 
    if mode == False:
        run_Itasser(args.sequence_file)
    else:
        run_HHblits(args.sequence_file)

def database_choice(base:str):
    if base.upper() == "UNICLUST30" or "PDB70":
        path = str(os.system("whereis {}".format(base)))
        print(path)
        return base, path
    else:
        print(" choosen database does not exist, please choose Uniclst30 or pdb70")
        print("  database will be Uniclst30")
        base = "Uniclust30"
        path = str(os.system("whereis {}".format(Uniclust30)))
        return path, base

def run_HHblits(ident,seq_file):
    """
        run HHblits
    """

    #os.popen("hhblits -i seq_file  -d DATABASE ....")
    #align_file = align_script(ident, args.sequence_file, args.output_file)
    #exécuter le script shell
    # Database choice
    path, base = database_choice(args.database)
    print("Run HHblits on {} \n".format(seq_file))
    os.system("conda activate {}".format(RESPIRE_ENV))
    os.system("hhblits -cpu 8 -i {} -d ~{}/{} -oa3m {}.a3m ".format(seq_file, path, base, ident)) #not finished yet to check with ferdi 
    #os.system("./{}_HHblits.sh".format(ident))
    
    # for ITASSER   
def itasser_prep_seq (id_seq: str, seq_str:str):
    """ 1.test length of the sequence for modelling
        sequence length have to be between 10 and 1500 bases 
        2. check for I-tasser file for outputs or created

    """
    if len(seq_str) < 10 or len(seq_str) > 1500:
            os.exit('{}: I-TASSER modeling will not run as sequence length ({}) is not in the range [10, 1500]'.format(id, len(seq_str)))

    if not op.exists("I-TASSER"):
            os.makedirs("I-TASSER")
def main():
    """
    Main program function
    """
    
    # get arguments
    args = get_arguments()
    print(args.sequence_file)
    #ident(sequence name), sequence
    #ident, seq_str = read_fasta(args.sequence_file)
    # Itasser mode 
    Itasser_run_mode(args.ITasser_mode)
    

    input_files = files_number(args.sequence_file)

    for file in sequence_file:
        ident, seq_str = read_fasta(files_input)
        run_HHblits(file, ident)

if __name__ == '__main__':
    main()
