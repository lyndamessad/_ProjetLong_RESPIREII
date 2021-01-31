
from module import *
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
log = logging.getLogger(__name__)

class ITASSERPrep():

    """Prepare a protein sequence for an I-TASSER homology modeling run.
    The main utilities of this class are to:
    * Allow for the input of a protein sequence string and paths to I-TASSER to create execution scripts
    * Automate large-scale homology modeling efforts by creating Slurm or TORQUE job scheduling scripts
    Args:
        ident: Identifier for your sequence. Will be used as the global ID (folder name, sequence name)
        seq_str: Sequence in string format
        root_dir: Local directory where I-TASSER folder will be created
        itasser_path: Path to I-TASSER folder, i.e. '~/software/I-TASSER4.4'
        itlib_path: Path to ITLIB folder, i.e. '~/software/ITLIB'
        execute_dir: Optional path to execution directory - use this if you are copying the homology models to
            another location such as a supercomputer for running
        light: If simulations should be limited to 5 runs
        runtype: How you will be running I-TASSER - local, slurm, or torque
        print_exec: If the execution script should be printed out
        java_home: Path to Java executable
        binding_site_pred: If binding site predictions should be run
        ec_pred: If EC number predictions should be run
        go_pred: If GO term predictions should be run
        additional_options: Any other additional I-TASSER options, appended to the command
        job_scheduler_header: Any job scheduling options, prepended as a header to the file
    """

    def __init__(self, ident, seq_str, root_dir, itasser_path, itlib_path,
                 execute_dir=None, light=True, runtype='local', print_exec=False, java_home=None,
                 binding_site_pred=False, ec_pred=False, go_pred=False, additional_options=None,
                 job_scheduler_header=None):
        
        self.id = ident
        self.seq_str = seq_str

        if len(self.seq_str) < 10 or len(self.seq_str) > 1500:
            log.warning('{}: I-TASSER modeling will not run as sequence length ({}) is not in the range [10, 1500]'.format(self.id, len(self.seq_str)))

        self.root_dir = root_dir
        if not op.exists(root_dir):
            os.makedirs(root_dir)
        if not execute_dir:
            # If no execute_dir is given, use the same dir as the created folder
            self.execute_dir = self.prep_folder(seq_str)
        elif execute_dir:
            orig_data_dir = self.prep_folder(seq_str)
            self.execute_dir = op.join(execute_dir, op.basename(orig_data_dir))

        self.print_exec = print_exec
        self.runtype = runtype
        if light:
            light = 'true'
        else:
            light = 'false'
        self.light = light

        self.model_exists = op.exists(op.join(self.execute_dir, 'model1.pdb'))

        if not additional_options:
            additional_options = ''
        else:
            additional_options += ' '
        if binding_site_pred:
            additional_options += '-LBS true '
        if ec_pred:
            additional_options += '-EC true '
        if go_pred:
            additional_options += '-GO true '
        self.additional_options = additional_options

        if not java_home:
            self.java_home = '${JAVA_HOME}'
        else:
            self.java_home = java_home

        if not job_scheduler_header:
            self.job_scheduler_header = ''
        else:
            self.job_scheduler_header = job_scheduler_header

        if runtype == 'local' or runtype == 'torque':
            self.prep_script_local(itasser_loc=itasser_path,
                                   itlib_loc=itlib_path)

        if runtype == 'slurm':
            self.prep_script_slurm(itasser_loc=itasser_path,
                                   itlib_loc=itlib_path)

    def prep_folder(self, seq):
        """Take in a sequence string and prepares the folder for the I-TASSER run."""
        itasser_dir = op.join(self.root_dir, self.id)

        if not op.exists(itasser_dir):
            os.makedirs(itasser_dir)

        tmp = {self.id: seq}

        fasta.write_fasta_file_from_dict(indict=tmp,
                                         outname='seq',
                                         outext='.fasta',
                                         outdir=itasser_dir)
        return itasser_dir


    def prep_script_local(self, itasser_loc, itlib_loc):
        script_file = '{}.sh'.format(self.id)
        outfile = os.path.join(self.root_dir, script_file)

        itasser = {'executable': op.join(itasser_loc, 'I-TASSERmod/runI-TASSER.pl'),
                   'pkgdir': itasser_loc,
                   'libdir': itlib_loc,
                   'seqname': self.id,
                   'datadir': self.execute_dir,
                   'java_home': self.java_home,
                   'additional_options': self.additional_options,
                   'light': self.light}

        script = open(outfile, 'w')
        script.write('#!/bin/bash -l\n')

        if self.runtype == 'torque':
            script.write('{}'.format(self.job_scheduler_header))

        script.write(("{i[executable]} "
                      "-pkgdir {i[pkgdir]} "
                      "-libdir {i[libdir]} "
                      "-seqname {i[seqname]} "
                      "-datadir {i[datadir]} "
                      "-java_home {i[java_home]} "
                      "{i[additional_options]}"
                      "-light {i[light]}\n\n").format(i=itasser))
        script.close()

        os.chmod(outfile, 0o755)

        if self.print_exec and self.runtype=='local':
            print('nohup ./{} > {}.out &'.format(op.basename(outfile), os.path.join(self.root_dir, self.id)),
                  end='\n\n')

        if self.print_exec and self.runtype == 'torque':
            print('qsub {}'.format(op.basename(outfile), os.path.join(self.root_dir, self.id)),
                  end='; ')

        return outfile

    def prep_script_slurm(self, itasser_loc, itlib_loc):
        script_file = '{}.slm'.format(self.id)
        outfile = os.path.join(self.root_dir, script_file)

        itasser = {'executable': op.join(itasser_loc, 'I-TASSERmod/runI-TASSER.pl'),
                   'pkgdir': itasser_loc,
                   'libdir': itlib_loc,
                   'seqname': self.id,
                   'datadir': self.execute_dir,
                   'java_home': self.java_home,
                   'light': self.light,
                   'additional_options': self.additional_options}

        slurm = open(outfile, 'w')

        slurm.write('#!/bin/bash -l\n')
        slurm.write('{}'.format(self.job_scheduler_header))
        slurm.write(('{i[executable]} '
                     '-pkgdir {i[pkgdir]} '
                     '-libdir {i[libdir]} '
                     '-seqname {i[seqname]} '
                     '-datadir {i[datadir]} '
                     '-java_home {i[java_home]} '
                     '{i[additional_options]}'
                     '-light {i[light]}\n\n').format(i=itasser))

        slurm.close()

        os.chmod(outfile, 0o755)

        if self.print_exec:
            print('sbatch {}'.format(op.basename(outfile)), end='; ')

        return outfile
 


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

    run = ITASSERPrep()
    run.__init__(self, ident, seq_str, args.root_dir, args.itasser_path, args.itlib_path,
                 execute_dir=None, light=True, runtype='local', print_exec=False, java_home=None,
                 binding_site_pred=False, ec_pred=False, go_pred=False, additional_options=None,
                 job_scheduler_header=None)
if __name__ == '__main__':
    main()

    
    """
    commande python : 
    python script_ITasser.py  sequence.fasta   output_file """

    """

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
    """
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
    """
    root = '/home/nathan/projects/GEM-PRO/cyano/'
    files = glob.glob(os.path.join(root,'*.faa'))
    for f in files:
        identifier = os.path.splitext(os.path.basename(f))[0]
        ip = ITASSERPrep(id=ident, root_dir=args.root_dir)
        #sequence = sl.seq_loader(f, is_file=True)
        execute_dir = ip.prep_folder(seq_str)
        ip.prep_script_local(itasser_loc='/home/username/I-TASSER4.4',
                            itlib_loc='/home/username/ITLIB',
                            datadir=execute_dir)
    import os
    import shutil
    import subprocess

    thedir = '../Proteines'
    print(thedir)
    folders = [name for name in os.listdir(thedir)]
    folders = sorted(folders, reverse=True)
    length = int(len(folders)) / int(2)
    print(length)
    for_ssb3 = folders[:int(length)]
    #
    for fo in for_ssb3:
        coach = open('%s_coach.sh' % fo, 'w')
        coach.write('#!/bin/bash\n')
        coach.write('#PBS -N %s\n' % fo)
        coach.write('-pkgdir /I-TASSERmod/runI-TASSER.pl \
                    -libdir /home/yourname/ITLIB \
                    -seqname %s \
                    -model model1.pdb \
                    -datadir /home/yourname/I-TASSER5.0/%s \
                    -GO true \
                    -runstyle parallel \
                    -ligth \n\n '% (fo, fo))
        coach.close()
        subprocess.call('qsub %s_coach.sh;' % (fo), shell=True)
        print('qsub %s_coach.sh;' % (fo)),
        