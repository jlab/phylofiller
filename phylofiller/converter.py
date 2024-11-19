import pandas as pd
from itertools import groupby
from Bio.Seq import Seq
import collections
import sys
from tqdm import tqdm
import re


def easel_table2pd(lines, verbose=True, add_header_names=dict()) -> pd.DataFrame:
    """
    Parameters
    ----------
    lines : [str]
        list of lines
    verbose : Boolean
        Print warnings
    add_header_names : dict(Column: Label)
        Some tables in easel can have unnamed columns. In this cases, you
        can provide a dictionary with column number as keys and column labels
        as headers, e.g.
        {1: 'significance', 11: 'strand'}
    """
    headers = []
    for m in re.finditer(r"\S+", lines[2]):
        header_start, header_end = m.start(), m.end()
        while lines[1][header_start] == "-":
            header_start -= 1
        while (header_end < len(lines[0])) and (lines[1][header_end] == "-"):
            header_end += 1
        headers.append(lines[0][header_start:header_end].strip())

    for col, label in add_header_names.items():
        if headers[col] != "":
            raise ValueError(
                ("The easel tabel column %i is names '%s', but you try to "
                 "override this with '%s'.") % (col, headers[col], label))
        headers[col] = label

    rows = []
    for row_num, line in tqdm(enumerate(lines), disable=not verbose):
        if row_num <= 1:
            continue
        if row_num > 1 and line.startswith('#'):
            continue  # as this is a trailing comment line
        rows.append(line.rstrip().split())

    table = pd.DataFrame(
        data=rows,
        columns=headers,
        index=range(len(rows))
        )

    # automatic datatype conversion from "object" to int and float to
    # saveguard the user from sorting implicitely lexiographically when
    # sorting numerically is intended
    for col in ['mdl from', 'mdl to', 'seq from', 'seq to', 'pass']:
        if col in table.columns:
            table[col] = table[col].astype(int)
    for col in ['gc', 'bias', 'score', 'E-value']:
        if col in table.columns:
            table[col] = table[col].astype(float)

    if (any(pd.Series(table.columns).value_counts() > 1)) and verbose:
        print("warning: column names are not unique!", file=sys.stderr)

    return table


def parse_easel_output(fp_input: str, verbose=True) -> pd.DataFrame:
    """Parses Infernal (soon also Hmmer) as a pandas DataFrame.

    Parameters
    ----------
    fp_input : str
        Filepath to Infernal output file.

    Returns
    -------
    Pandas.DataFrame holding all information of the infernal run.
    """
    step_lines_alignment = 6
    line_offset_endaln = 1

    program_info = dict()
    with open(fp_input, 'r') as f:
        lines = f.readlines()
        model_name = None
        model_clen = None
        hits = []
        for line_number in range(len(lines)):
            line = lines[line_number]
            if line.startswith('# '):
                if ' INFERNAL ' in line:
                    program_info['software'] = line.split()[1]
                    program_info['software version'] = line.split()[2]
                    # 1.1.2 --> 2.1.1 --> 2*10^0 + 1*10^1 + 1*10^2 --> 112
                    version_fingerprint = sum(
                        map(lambda x: int(x[0])*(10**x[1]),
                            zip(reversed(
                                program_info['software version'].split('.')),
                                range(10))))
                    if (version_fingerprint >= 113) and \
                       (version_fingerprint < 114):
                        step_lines_alignment = 7
                        line_offset_endaln = 2
                elif ' query CM file:' in line:
                    program_info['fp_query'] = line.split(':')[1].strip()
                elif ' target sequence database:' in line:
                    program_info['fp_target'] = line.split(':')[1].strip()
            elif line.startswith('Query:'):
                model_name = line.split()[1]
                model_clen = line.split()[-1].split('=')[-1].replace(']', '')
            elif line.startswith('>> '):
                # next three lines will be an overview table
                hit_info = easel_table2pd(
                    lines[line_number+1:line_number+3+1],
                    verbose=verbose, add_header_names={
                        1: 'significance', 11: 'target strand',
    # from Infernal documentation, page 22
    # It’s not immediately easy to tell from the “to” coordinate whether the
    # alignment ended internally in the query model or target sequence, versus
    # ran all the way (as in a full-length global alignment) to the end(s).
    # To make this more readily apparent, with each pair of query and target
    # endpoint coordinates, there’s also a little symbology. For the normal
    # case of a non-truncated hit: .. means both ends of the alignment ended
    # internally, and [] means both ends of the alignment were full-length
    # flush to the ends of the query or target, and [. and .] mean only the
    # left (5’) or right (3’) end was flush/full length. For truncated hits,
    # the symbols are the same except that either the first and/or the second
    # symbol will be a ˜ for the query and target. If the first symbol is ˜
    # then the left end of the alignment is truncated because the 5’ end of
    # the hit is predicted to be missing (extend beyond the beginning of the
    # target sequence). Similarly, if the second symbol is ˜ then the right
    # end of the alignment is truncated because the 3’ end of the hit is
    # predicted to be missing (extend beyond the end of the target sequence).
    # These two symbols occur just after the “mdl to” column for the query,
    # and after the strand + or - symbol for the target sequence.
                        8: 'mdl truncation', 12: 'seq truncation'})
                # target name is in the >> line
                hit_info['target name'] = line[3:].strip()

                alignment_line_number = line_number + step_lines_alignment
                query_sequence = ""
                target_sequence = ""
                # print("start", step_lines_alignment, alignment_line_number,
                # lines[alignment_line_number])
                while alignment_line_number + 2 < len(lines):
                    query_sequence += lines[alignment_line_number].split()[2]
                    # print("hier", step_lines_alignment,
                    #     line_offset_endaln,
                    #     alignment_line_number+1,
                    #     lines[alignment_line_number+1][:30], file=sys.stderr)
                    target_sequence += lines[
                        alignment_line_number+2].split()[2]
                    alignment_line_number += step_lines_alignment
                    # print(alignment_line_number+1,
                    #       lines[alignment_line_number][:30], file=sys.stderr)
                    if lines[alignment_line_number -
                             line_offset_endaln].startswith('>> '):
                        break
                    elif ((lines[alignment_line_number +
                                 (2 - line_offset_endaln)].startswith(
                            'Internal HMM-only pipeline statistics summary'))
                          or (lines[alignment_line_number +
                                    (2 - line_offset_endaln)].startswith(
                            'Internal CM pipeline statistics summary'))):
                        break
                hit_info['query sequence'] = query_sequence
                hit_info['target sequence'] = target_sequence

                hit_info['model name'] = model_name
                hit_info['model clen'] = model_clen
                hits.append(hit_info)

    # handle situations with zero hits
    if len(hits) > 0:
        hits = pd.concat(hits)
    else:
        hits = pd.DataFrame()

    for field in program_info.keys():
        hits[field] = program_info[field]
    hits.index = range(hits.shape[0])

    return hits


def create_CIGAR(row_reference: str, row_read: str) -> str:
    """Given reference and read alignment row, convert into CIGAR string.

    Note: result is NOT symmetric!

    Parameters
    ----------
    row_reference : str
        The alignment row for the reference, e.g. the genome
    row_read : str
        The alignment row for the read, e.g. a true Illumina read or the
        Infernal model consensus sequence.

    Returns
    -------
    str: the CIGAR string
    """
    GAP = '-'
    # see https://samtools.github.io/hts-specs/SAMv1.pdf for a definition
    # of alignment operators
    assert len(row_reference) == len(row_read), \
        'Reference and read string have different length.'
    assert isinstance(row_reference, str), "Reference is not of type str."
    assert isinstance(row_read, str), "Read is not of type str."

    operations = []
    for i in range(len(row_reference)):
        char_ref = row_reference[i]
        char_read = row_read[i]
        if char_ref in ['.', '-']:
            char_ref = GAP
        if char_read in ['.', '-']:
            char_read = GAP

        if (char_ref != GAP) and (char_read != GAP):
            if char_ref == char_read:
                operations.append('=')
            else:
                operations.append('X')
        elif (char_ref == GAP) and (char_read != GAP):
            operations.append('I')
        elif (char_ref != GAP) and (char_read == GAP):
            operations.append('D')

    runs = [(i, len(list(g))) for i, g in groupby(operations)]
    return ''.join((map(lambda x: '%i%s' % (x[1], x[0]), runs)))


def easle2sam(easel_data: pd.DataFrame) -> str:
    """Converts Infernal output in pd.DataFrame format into SAM string.

    Parameters
    ----------
    easel_data : pd.DataFrame
        The result of parse_easel_output.

    Returns
    -------
    Infernal hits in SAM format as a str.
    """
    # convert datatypes
    hits = easel_data.copy()
    for field in ['mdl from', 'mdl to', 'seq from', 'seq to', 'model clen']:
        hits[field] = hits[field].astype(int)
    for field in ['E-value', 'score', 'bias', 'acc', 'gc']:
        hits[field] = hits[field].astype(float)
    # extract rank number, e.g. (1) becomes 1
    hits['rank'] = hits['rank'].apply(lambda x: int(x[1:-1]))

    # convert T to U
    for field in ['query sequence', 'target sequence']:
        hits[field] = hits[field].apply(lambda x: x.replace('U', 'T').replace(
            'u', 't'))

    # add strand information
    hits['strand'] = hits.iloc[:, 11].apply(
        lambda x: 'rev' if x[0] == '-' else 'fwd')

    # reverse complement if hit is on reverse strand
    for idx in hits.index:
        if (hits.loc[idx, 'seq from'] > hits.loc[idx, 'seq to']) and \
           (hits.loc[idx, 'strand'] == 'rev'):
            for field in ['query sequence', 'target sequence']:
                hits.loc[idx, field] = \
                    str(Seq(hits.loc[idx, field]).reverse_complement())
            helpswap = hits.loc[idx, 'seq from']
            hits.loc[idx, 'seq from'] = hits.loc[idx, 'seq to']
            hits.loc[idx, 'seq to'] = helpswap

    # create CIGAR string
    hits['cigar'] = hits.apply(lambda row: create_CIGAR(
        row['target sequence'], row['query sequence']), axis=1)

    # extend hit to model boundaries, for better visualization in genome
    # browsers
    for idx in hits.index:
        leading = hits.loc[idx, 'mdl from'] - 1
        trailing = hits.loc[idx, 'model clen'] - hits.loc[idx, 'mdl to']
        if hits.loc[idx, 'strand'] == 'rev':
            leading, trailing = trailing, leading
        if leading > 0:
            hits.loc[idx, 'seq from'] -= (leading + 1)
            hits.loc[idx, 'cigar'] = '1M%iD%s' % (
                leading, hits.loc[idx, 'cigar'])
            hits.loc[idx, 'query sequence'] = \
                '?%s' % hits.loc[idx, 'query sequence']
        if trailing > 0:
            hits.loc[idx, 'cigar'] += '%iD' % trailing

    # convert table into SAM formatted strings
    samrows = []
    for i, (idx, hit) in enumerate(hits.sort_values([
            'target name', 'seq from']).iterrows()):
        sam = []
        sequence = hit['query sequence'].replace('-', '').replace('.', '')
        sam.append('%s:%i' % (hit['model name'], i))  # QNAME
        sam.append("0" if hit.iloc[11][0] == '+' else '16')  # FLAG
        sam.append(hit['target name'])  # RNAME
        sam.append(str(hit['seq from']))  # POS
        sam.append(str(60))  # MAPQ
        sam.append(hit['cigar'])  # CIGAR
        sam.append('*')  # RNEXT
        sam.append(str(0))  # PNEXT
        sam.append(str(0))  # TLEN
        sam.append(sequence)  # SEQ
        sam.append('*')  # QUAL
        sam.append('RG:Z:%s' % hit['software'])
        sam.append('CM:f:%f' % hit['score'])
        samrows.append('\t'.join(sam))

    return '\n'.join(samrows)+"\n"


def kreport2feature(fp_input, rank='Genus') -> pd.DataFrame:
    RANKS = collections.OrderedDict({
        'Unclassified': 'U',
        'Root': 'R',
        #  'Domain': 'D',
        'Kingdom': 'D',  # SMJ: count not find examples with K lines in reports
        'Phylum': 'P',
        'Class': 'C',
        'Order': 'O',
        'Family': 'F',
        'Genus': 'G',
        'Species': 'S'})
    """
    https://github.com/DerrickWood/kraken2/wiki/Manual#
    sample-report-output-format
    1. Percentage of fragments covered by the clade rooted at this taxon
    2. Number of fragments covered by the clade rooted at this taxon
    3. Number of fragments assigned directly to this taxon
    4. A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom,
       (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. Taxa that
       are not at any of these 10 ranks have a rank code that is formed by
       using the rank code of the closest ancestor rank with a number
       indicating the distance from that rank. E.g., "G2" is a rank code
       indicating a taxon is between genus and species and the grandparent
       taxon is at the genus rank.
    5. NCBI taxonomic ID number
    6. Indented scientific name

    """
    assert (rank in RANKS.keys()) or rank is None

    report = pd.read_csv(
        fp_input, sep="\t", header=None,
        names=['rel.abund.', '#reads_clade', '#reads_taxon', 'rank',
               'taxid', 'scientific name'])

    # determine lineage for major ranks and taxids
    last_level = -1
    curr_lineage = []
    curr_taxids = []
    res = []
    res_taxids = []
    for i, (_, row) in enumerate(report.iterrows()):
        if row['rank'] not in RANKS.values():
            res.append(None)
            res_taxids.append(None)
            continue
        if (row['rank'] == 'U') and (row['scientific name'] == 'unclassified'):
            res.append(row['scientific name'])
            res_taxids.append(row['taxid'])
            continue

        rank_level = list(RANKS.values()).index(row['rank'])
        rank_name = '%s__%s' % (
            list(RANKS.values())[rank_level].lower().replace('d', 'k'),
            row['scientific name'].strip())
        rank_taxid = row['taxid']
        if rank_level > last_level:
            curr_lineage.append(rank_name)
            curr_taxids.append(rank_taxid)
        elif rank_level == last_level:
            curr_lineage[-1] = rank_name
            curr_taxids[-1] = rank_taxid
        else:
            curr_lineage = curr_lineage[:rank_level-1] + [rank_name]
            curr_taxids = curr_taxids[:rank_level-1] + [rank_taxid]
        res.append(';'.join(curr_lineage))
        res_taxids.append(curr_taxids)
        last_level = rank_level
    report['lineage'] = res
    report['taxid_lineage'] = res_taxids

    if rank is not None:
        report = report[report['rank'] == RANKS[rank]]

    return report
