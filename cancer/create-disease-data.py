"""
Convert p53 disease data from nice tsv into an opaque database format.

The input file to this script should be the text file extracted
from the UMD TP53 mutation database,
in particular the curated tp53 database tumours only, US txt version.
http://p53.fr/TP53_database_download/TP53_tumor_database/
TP53_US/UMDTP53_curated_2012_R1_US.txt.zip
The wild codon state is provided for more redundancy.
"""

import argparse
import sqlite3
import csv


def main(args):

    # read the genetic code
    conn = sqlite3.connect(args.code)
    cursor = conn.cursor()
    cursor.execute(
        "select state, residue, codon from code "
        "where residue <> 'Stop' "
        "order by state")
    genetic_code = list(cursor)
    cursor.execute(
        "select codon, state from code "
        "where residue <> 'Stop'")
    codon_to_state = dict(cursor)
    conn.close()

    # Read the tumor info csv file.
    table_tuples = set()
    codon_offset_index = 3
    wild_codon_index = 4
    mutant_codon_index = 5
    wild_aa_index = 6
    mutant_aa_index = 7
    with open(args.infile, 'rU') as fin:
        reader = csv.reader(fin,dialect=csv.excel_tab)
        headers = reader.next()
        for row in reader:

            # filter out mutations that are clearly not amino acid changes
            wild_aa = row[wild_aa_index]
            mutant_aa = row[mutant_aa_index]
            if mutant_aa in ('Stop', 'Fs.'):
                continue
            if wild_aa == mutant_aa:
                continue

            # get the wild codon state
            wild_codon = row[wild_codon_index]
            wild_state = codon_to_state.get(wild_codon, None)
            if wild_state is None:
                continue

            # attempt to add the mutation to the list
            mutant_codon = row[mutant_codon_index]
            mutant_state = codon_to_state.get(mutant_codon, None)
            if mutant_state is None:
                continue

            codon_offset = int(row[codon_offset_index])
            table_tuples.add((codon_offset, wild_state, mutant_state))

    # open the tumor info database for writing
    conn = sqlite3.connect(args.outfile)
    cursor = conn.cursor()

    # Initialize the tumor info table.
    # The mutation id is defined by the UMD p53 mutation database.
    # The offset is the codon offset.
    # The state is the codon defined state defined in the genetic code.
    cursor = conn.cursor()
    s = (
            'create table if not exists tumor ('
            'offset integer, '
            'wild_state integer, '
            'mutant_state integer, '
            'primary key (offset, wild_state, mutant_state))')
    cursor.execute(s)
    conn.commit()

    # populate the primary state history table
    for t in table_tuples:
        s = 'insert into tumor values (?, ?, ?)'
        cursor.execute(s, t)
    conn.commit()

    # close the alignment database table
    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--infile', required=True,
            help='input p53 tumor file')
    parser.add_argument('--code', default='universal.code.db',
            help='genetic code as an sqlite3 database file')
    parser.add_argument('--outfile', default='tumor.db',
            help='output disease data as an sqlite3 database file')
    main(parser.parse_args())

