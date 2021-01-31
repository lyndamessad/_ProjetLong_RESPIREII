 
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
from ITASSERPrep import *

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
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-file', '-sequence_file', dest='sequence_file', type=isfile, required=True,
                        help="Sequence is a fasta file (.fasta or .fasta.gz)")

    parser.add_argument('-t', '-runtype', dest='runtype', type=str, default = 'local',
                        help="How you will be running I-TASSER - local, slurm, or torque")

    parser.add_argument('-r', '-root_dir', dest='root_dir', type=str, default = os.getcwd(),
                        help="Local directory where I-TASSER folder will be created")

    parser.add_argument('-p', '-itasser_path', dest='itasser_path', type=str,
                        help="Path to I-TASSER folder")

    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="flag.txt", help="Output file")
    return parser.parse_args()

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
            return seqDict

    # for no gziped fasta files:
    else:
        seqRecord = SeqIO.read(sequence_file, "fasta")
        sequence = seqRecord.seq
        ident = seqRecord.id
        return ident, sequence

def main():
    """
    Main program function
    """
    
    # get arguments
    args = get_arguments()

    #ident(sequence name), sequence
    ident, seq_str = read_fasta(args.sequence_file)
    print(seq_str)
    
    # How you will be running I-TASSER - local, slurm, or torque
    runtype = args.runtype
    if runtype.lower() not in ['local', 'torque', 'slurm']:
        raise ValueError('Invalid runtype, must be "local", "torque", "slurm"')
    print(args.root_dir)


if __name__ == '__main__':
    main()

    run = ITASSERPrep()
    run.__init__(self, ident, seq_str, root_dir, itasser_path, itlib_path,
                 execute_dir=None, light=True, runtype='local', print_exec=False, java_home=None,
                 binding_site_pred=False, ec_pred=False, go_pred=False, additional_options=None,
                 job_scheduler_header=None)
    """
    commande python : 
    python script_ITasser.py  sequence.fasta   output_file """

    """
        #Local directory where I-TASSER folder will be created
        root_dir =
        #Path to I-TASSER folder
        itasser_path = 
        # Path to ITLIB folde
        itlib_path = 

        pass
    """ 
    # TODO: make this an executable script to
    # 1) ask for global I-TASSER locations
    # 2) ask for working directory
    # 3) take in multiple inputs and prepare them for I-TASSER runs
    #     a) input types
    #         i) a single FASTA file with single or multiple sequences
    #         ii) multiple FASTA files contained in the working directory
    #         iii) a dataframe with IDs and sequences
    #         iv) a sequence string and an ID (and optional additional identifiers)
    #     b) types of runs
    #         i) NERSC slurm (sbatch) inputs
    #         ii) local torque (qsub) inputs
    #         iii) simple executable background scripts
    # 4) Output executable scripts or submit things to the queue

    # root = '/home/nathan/projects/GEM-PRO/cyano/'
    # files = glob.glob(os.path.join(root,'*.faa'))
    # for f in files:
    #     identifier = os.path.splitext(os.path.basename(f))[0]
    #     ip = ITASSERPrep(id=identifier, root_dir='/home/nathan/projects/GEM-PRO/cyano')
    #
    #     sequence = sl.seq_loader(f, is_file=True)
    #     execute_dir = ip.prep_folder(sequence)
    #     ip.prep_script_local(itasser_loc='/home/nathan/software/I-TASSER4.4',
    #                          itlib_loc='/home/nathan/software/ITLIB',
    #                          datadir=execute_dir)

    # ip = ITASSERPrep(id='W5EP13', root_dir='/home/nathan/projects/GEM-PRO/cyano/')
    #
    # sequence = sl.seq_loader('/home/nathan/Downloads/W5EP13.faa', is_file=True)
    # execute_dir = ip.prep_folder(sequence)
    # ip.prep_script_local(itasser_loc='/home/nathan/software/I-TASSER4.4',
    #                      itlib_loc='/home/nathan/software/ITLIB',
    #                      datadir=execute_dir)


## below is old run_all script in python
# import os
# import shutil
# import subprocess
#
# thedir = '.'
# folders = [name for name in os.listdir(
#     thedir) if os.path.isdir(os.path.join(thedir, name))]
# folders = sorted(folders, reverse=True)
# for_ssb3 = folders[:len(folders) / 2]
#
# for fo in for_ssb3:
#     coach = open('%s_coach.sh' % fo, 'w')
#
#     coach.write('#!/bin/bash\n')
#     coach.write('#PBS -l walltime=05:20:00\n')
#     coach.write('#PBS -q regular\n')
#     coach.write('#PBS -N %s\n' % fo)
#     coach.write('perl ~/software/I-TASSER4.4/I-TASSERmod/runCOACH.pl -pkgdir /home/nathan/software/I-TASSER4.4 -libdir /home/nathan/software/ITLIB -protname %s -model model1.pdb -datadir /home/nathan/projects/GEM-PRO/yome/all_test/%s -GO true\n\n' % (fo, fo))
#
#     coach.close()
#
#     # subprocess.call('qsub %s_coach.sh;' % (fo), shell=True)
#     print('qsub %s_coach.sh;' % (fo)),
