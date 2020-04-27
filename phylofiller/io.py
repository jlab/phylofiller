import hashlib
import os

import skbio.io
import pandas as pd


def get_assembly_checksum(fp_assembly: str):
    """Get checksum for assemlby sequences.

    The MD5 checksum ignores sequence formatting (linewrap, case), ordering
    and header information alltogether.

    Parameters
    ----------
    fp_assembly : str
        Filepath to fasta or fasta.gz (multiple) sequence file.

    Returns
    -------
    md5 sum of lexicographically sorted, uppercase sequences
    (ignoring formatting, headers and comments)
    """
    lex_sorted_upper_seqs = sorted(
        [str(seq).upper()
         for seq
         in skbio.io.read(fp_assembly, format='fasta')])
    hash_value = hashlib.md5(str(lex_sorted_upper_seqs).encode()).hexdigest()

    return hash_value
