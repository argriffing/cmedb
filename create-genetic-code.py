"""
This script converts genetic code text files into sqlite3 format.

For example
python create-genetic-code.py
--infile=universal.code.txt
--outfile=universal.code.db
"""

import argparse
import sqlite3

def main(args):

    # open the database file and create the table
    conn = sqlite3.connect(args.outfile)
    cursor = conn.cursor()
    cursor.execute(
            'create table code ('
            'state integer, residue text, codon text, '
            'primary key (state))')
    conn.commit()

    # read the genetic code from the text file and put it into the database
    with open(args.infile) as fin:
        for line in fin.readlines():
            s_state, s_residue, s_codon = line.split()
            state = int(s_state)
            residue = s_residue[0].upper() + s_residue[1:].lower()
            codon = s_codon.upper()
            triple = (state, residue, codon)
            cursor.execute('insert into code values (?, ?, ?)', triple)
    conn.commit()

    # close the database connection
    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', required=True,
            help='genetic code as a text file in an ad hoc format')
    parser.add_argument('-o', '--outfile', required=True,
            help='create this sqlite3 database file')
    main(parser.parse_args())

