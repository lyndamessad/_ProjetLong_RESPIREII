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

def test_seq_rep(seq_id: str, seq_str:str):
    """ 1.test length of the sequence for modelling
        sequence length have to be between 10 and 1500 bases 
        2. check for I-tasser file for outputs or created

    """
  if len(seq_str) < 10 or len(seq_str) > 1500:
    os.exit('{}: I-TASSER modeling will not run as sequence length ({}) is not in the range [10, 1500]'.format(id, len(seq_str)))

  if not op.exists("I-TASSER"):
    os.makedirs("I-TASSER")


def prep_script_local(itasser_loc, itlib_loc, seq_id,execute_dir,java_home,runtype =='local'):
  runtype = Itasser_run_mode(args.ITasser_mode)
  script_file = '{}.sh'.format(seq_id)
  outfile = os.path.join(root_dir, script_file)

  itasser = {'executable': op.join(itasser_loc, 'I-TASSERmod/runI-TASSER.pl'),
             'pkgdir': itasser_loc,
             'libdir': itlib_loc,
             'seqname':seq_id,
             'datadir': execute_dir, #????
             'java_home':java_home,  #???? 
             'additional_options': additional_options,
                   #'light': self.light}

  script = open(outfile, 'w')
  script.write('#!/bin/bash -l\n')

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

  if runtype =='local':
    print('nohup ./{} > {}.out & \n\n'.format(op.basename(outfile), os.path.join(root_dir, seq_id)))
    """ 
    if self.print_exec and runtype == 'torque':
      print('qsub {}'.format(op.basename(outfile), os.path.join(root_dir, seq_id)),
      end='; ')
    """

  return outfile
