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


def read_metadata(fp_metadata: str, fp_assemblyprefix: str = None):
    """.

    Parameters
    ----------
    fp_metadata : str
    fp_assemblyprefix : str

    Returns
    -------

    Raises
    ------
    Diverse erros if misconfigurations are present.
    """
    EXPECTED_COLUMNS = ['organism', 'fp_assembly']
    RESERVED_COL_NAMES = ['__line_number', '__file_exists']

    meta = pd.read_csv(fp_metadata, sep="\t", index_col=0)

    conflicting_col_names = set(RESERVED_COL_NAMES) & set(meta.columns)
    if len(conflicting_col_names) > 0:
        raise ValueError(
            'Conflicting column name(s) in your metadata: %s please rename!' %
            ', '.join(map(lambda x: '"%s"' % x, conflicting_col_names)))

    if meta.index.name != EXPECTED_COLUMNS[0]:
        raise ValueError(
            'Header of first column must be "%s"' % EXPECTED_COLUMNS[0])

    # temporarily add line numbers to metadata
    meta['__line_number'] = range(2, meta.shape[0]+2)

    # check if organisms have been defined multiple times
    ambig_organisms = meta.loc[meta.index.value_counts() > 1][
        '__line_number'].reset_index().groupby('organism').apply(
            lambda x: x['__line_number'].values)
    if ambig_organisms.shape[0] > 0:
        raise ValueError(
            '%i organisms have been re-defined in your '
            'metadata. Please fix!\n%s' %
            (ambig_organisms.shape[0], str(ambig_organisms)))

    # add assembly prefix if given
    if fp_assemblyprefix is not None:
        meta['fp_assembly'] = meta['fp_assembly'].apply(
            lambda x: os.path.abspath('%s%s' % (fp_assemblyprefix, x)))
    # check if given assembly file paths exists
    meta['__file_exists'] = meta['fp_assembly'].apply(
        lambda x: os.path.exists(x))
    if meta[~meta['__file_exists']].shape[0] > 0:
        raise ValueError(
            "Cannot find files for following assemblies:\n%s" %
            str(meta[~meta['__file_exists']]))

    del meta['__line_number']

    return meta
