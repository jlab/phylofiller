import pandas as pd

def easel_table2pd(lines):
    # easel is using a strange tab format. Instead of splitting via a single del. char like \t
    # columns are aligned for visual inspection. Hints of column starts/ends are given by the
    # second line with dashes, but it does not cover ALL columns :-/
    # Thus, I here try to make a consistent annotation of column fields (-) and spearators (#)
    line_fields = lines[1].replace('- ', '-#').replace(' -', '#-').replace(' ', '-')
    # leading whitespaces most likely indicate NO empty column.
    # Thus, I add these positions to the first real column.
    if line_fields.startswith('#'):
        line_fields = '-'+line_fields[1:]

    column_positions = []
    curr_pos = 0
    for i, field in enumerate(line_fields.split('#')):
        column_positions.append((curr_pos, curr_pos + len(field)))
        curr_pos += len(field) + 1

    rows = []
    for row_num, line in enumerate(lines):
        if row_num == 1:
            # this is the ---- line and no content
            continue
        row = []
        for col_num, (start, stop) in enumerate(column_positions):
            if col_num == 0:
                row.append(line[:stop].strip())
            elif col_num+1 == len(column_positions):
                row.append(line[start-1:].strip())
            else:
                row.append(line[start-1:stop].strip())
        rows.append(row)

    table = pd.DataFrame(data=rows[1:], columns=rows[0], index=range(len(rows)-1))

    return table


def easel2sam(fp_input: str):
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
                elif ' query CM file:' in line:
                    program_info['fp_query'] = line.split(':')[1].strip()
                elif ' target sequence database:' in line:
                    program_info['fp_target'] = line.split(':')[1].strip()
            elif line.startswith('Query:'):
                model_name = line.split()[1]
                model_clen = line.split()[-1].split('=')[-1].replace(']', '')
            elif line.startswith('>> '):
                # next three lines will be an overview table
                hit_info =  easel_table2pd(lines[line_number+1:line_number+3+1])

                # target name is in the >> line
                hit_info['target name'] = line[3:].strip()

                alignment_line_number = line_number + 6
                query_sequence = ""
                target_sequence = ""
                while alignment_line_number + 2 < len(lines):
                    query_sequence += lines[alignment_line_number].split()[2]
                    target_sequence += lines[alignment_line_number+2].split()[2]
                    alignment_line_number += 6
                    if lines[alignment_line_number-1].startswith('>> '):
                        break
                    elif lines[alignment_line_number+1].startswith('Internal HMM-only pipeline statistics summary'):
                        break
                hit_info['query sequence'] = query_sequence
                hit_info['target sequence'] = target_sequence

                hit_info['model name'] = model_name
                hit_info['model clen'] = model_clen
                hits.append(hit_info)

    hits = pd.concat(hits)
    for field in program_info.keys():
        hits[field] = program_info[field]

    return hits
