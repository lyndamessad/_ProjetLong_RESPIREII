""" module contenant les fonctions utiles au script I-tasser"""

from __future__ import print_function
import os
import sys
import os.path as op
import glob
import logging
import gzip
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

log = logging.getLogger(__name__)

def is_non_zero_file(fpath):
    """Check if a file exists, and that it contains some kind of contents
    Args:
        fpath (str): Path to file
    Returns:
        bool: If file exists and contains contents
    """
    return op.isfile(fpath) and os.path.getsize(fpath) > 0


def force_rerun(flag, outfile):
    """Check if we should force rerunning of a command if an output file exists.
    Args:
        flag (bool): Flag to force rerun.
        outfile (str): Path to output file which may already exist.
    Returns:
        bool: If we should force rerunning of a command
    Examples:
        >>> force_rerun(flag=True, outfile='/not/existing/file.txt')
        True
        >>> force_rerun(flag=False, outfile='/not/existing/file.txt')
        True
        >>> force_rerun(flag=True, outfile='./utils.py')
        True
        >>> force_rerun(flag=False, outfile='./utils.py')
        False
    """
    # If flag is True, always run
    if flag:
        return True
    # If flag is False but file doesn't exist, also run
    elif not flag and not op.exists(outfile):
        return True
    # If flag is False but filesize of output is 0, also run
    elif not flag and not is_non_zero_file(outfile):
        return True
    # Otherwise, do not run
    else:
        return False

def cast_to_seq_record(obj, alphabet=IUPAC.extended_protein, id="<unknown id>", name="<unknown name>",
                       description="<unknown description>", dbxrefs=None,
                       features=None, annotations=None,
                       letter_annotations=None):
    """Return a SeqRecord representation of a string or Seq object.
    Args:
        obj (str, Seq, SeqRecord): Sequence string or Biopython Seq object
        alphabet: See Biopython SeqRecord docs
        id: See Biopython SeqRecord docs
        name: See Biopython SeqRecord docs
        description: See Biopython SeqRecord docs
        dbxrefs: See Biopython SeqRecord docs
        features: See Biopython SeqRecord docs
        annotations: See Biopython SeqRecord docs
        letter_annotations: See Biopython SeqRecord docs
    Returns:
        SeqRecord: SeqRecord representation of the sequence
    """

    if isinstance(obj, SeqRecord):
        return obj
    if isinstance(obj, Seq):
        return SeqRecord(obj, id, name, description, dbxrefs, features, annotations, letter_annotations)
    if isinstance(obj, str):
        obj = obj.upper()
        return SeqRecord(Seq(obj, alphabet), id, name, description, dbxrefs, features, annotations, letter_annotations)
    else:
        raise ValueError('Must provide a string, Seq, or SeqRecord object.')


def write_fasta_file_from_dict(indict, outname, outdir=None, outext='.faa', force_rerun=False):
    """Write a FASTA file for a dictionary of IDs and their sequence strings.
    Args:
        indict: Input dictionary with keys as IDs and values as sequence strings
        outname: Name of the output file which will have outext appended to it
        outdir: Path to directory to output sequences to
        outext: Extension of FASTA file, default ".faa"
        force_rerun: If file should be overwritten if it exists
    Returns:
        str: Path to output FASTA file.
    """

    if not outdir:
        outdir = ''
    outfile = outfile_maker(inname='', outname=outname, outdir=outdir, outext=outext)

    if force_rerun(flag=force_rerun, outfile=outfile):
        seqs = []
        for i, s in indict.items():
            seq = cast_to_seq_record(s, id=i)
            seqs.append(seq)
        SeqIO.write(seqs, outfile, "fasta")

    return outfile


def outfile_maker(inname, outext='.out', outname='', outdir='', append_to_name=''):
    """Create a default name for an output file based on the inname name, unless a output name is specified.
    Args:
        inname: Path to input file
        outext: Optional specified extension for output file (with the "."). Default is ".out".
        outfile: Optional specified name of output file.
        outdir: Optional path to output directory.
    Returns:
        str: Path to final output destination.
    Examples:
        >>> outfile_maker(inname='P00001.fasta')
        'P00001.out'
        >>> outfile_maker(inname='/test/other/dir/P00001.fasta', outname='P00001_aligned', outdir='/my/dir/')
        '/my/dir/P00001_aligned.out'
    """

    # TODO: CHECK IF OUTNAME IS A VALID FILE NAME!
    orig_dir, orig_name, orig_ext = split_folder_and_path(inname)

    # If output filename not provided, default is to take name of inname
    if not outname:
        outname = orig_name

    # Create new path in the same directory of old path if a new one isn't specified
    if not outdir:
        outdir = orig_dir

    # Append additional stuff to the filename if specified
    if append_to_name:
        outname += append_to_name

    # Join the output filename and output extension
    final_outfile = op.join(outdir, '{}{}'.format(outname, outext))

    return final_outfile

def split_folder_and_path(filepath):
    """Split a file path into its folder, filename, and extension
    Args:
        path (str): Path to a file
    Returns:
        tuple: of (folder, filename (without extension), extension)
    """
    dirname = op.dirname(filepath)
    filename = op.basename(filepath)
    splitext = op.splitext(filename)
    filename_without_extension = splitext[0]
    extension = splitext[1]

    return dirname, filename_without_extension, extension

def gunzip_file(infile, outfile=None, outdir=None, delete_original=False, force_rerun_flag=False):
    """Decompress a gzip file and optionally set output values.
    Args:
        infile: Path to .gz file
        outfile: Name of output file
        outdir: Path to output directory
        delete_original: If original .gz file should be deleted
        force_rerun_flag: If file should be decompressed if outfile already exists
    Returns:
        str: Path to decompressed file
    """
    if not outfile:
        outfile = infile.replace('.gz', '')

    if not outdir:
        outdir = ''
    else:
        outdir = op.dirname(infile)
    outfile = op.join(outdir, op.basename(outfile))

    if force_rerun(flag=force_rerun_flag, outfile=outfile):
        gz = gzip.open(infile, "rb")
        decoded = gz.read()

        with open(outfile, "wb") as new_file:
            new_file.write(decoded)

        gz.close()
        log.debug('{}: file unzipped'.format(outfile))
    else:
        log.debug('{}: file already unzipped'.format(outfile))

    if delete_original:
        os.remove(infile)

    return outfile


if __name__ == '__main__':
