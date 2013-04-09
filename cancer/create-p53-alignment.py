"""
Convert from a standard phylip into an opaque database format.

The output alignment database file will have
leaf nodes indexed according to the input tree database file.
The output codon states will be indexed according
to the codon states in the input genetic code database file.
I tried to use dendropy to read this phylip file but it did made errors.
"""

import argparse
import sqlite3


def gen_paragraphs(lines):
    para = []
    for line in lines:
        line = line.strip()
        if not line:
            if para:
                yield para
                para = []
        else:
            para.append(line)
    if para:
        yield para


def read_phylip(fin):
    """
    Yield (taxon name, codons) pairs.
    @param fin: file open for reading
    """

    # Get the paragraphs in the most inefficient way possible.
    # Ignore the first line which is also the first paragraph.
    paras = list(gen_paragraphs(fin))[1:]
    if len(paras) != 25:
        raise Exception('expected p53 alignment of 25 taxa')

    # Each paragraph defines a p53 coding sequence of some taxon.
    # The first line gives the taxon name.
    # The rest of the lines are codons.
    for para in paras:
        taxon_name = para[0]
        codons = ' '.join(para[1:]).split()
        if len(codons) != 393:
            raise Exception('expected 393 codons')
        yield taxon_name, codons


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

    # read the taxon table from the tree database file
    conn = sqlite3.connect(args.tree)
    cursor = conn.cursor()
    name_to_taxon = dict(cursor.execute('select name, node from taxa'))
    conn.close()

    # read the phylip alignment
    alignment_triples = []
    with open(args.infile) as fin:
        for info in read_phylip(fin):
            taxon_name, codons = info
            for i, codon in enumerate(codons):

                # in biology we know that the first thing in a list is thing 1
                offset = i + 1
                taxon = name_to_taxon[taxon_name]
                state = codon_to_state[codon]
                triple = (offset, taxon, state)
                alignment_triples.append(triple)

    # open the alignment database for writing
    conn = sqlite3.connect(args.outfile)
    cursor = conn.cursor()

    # initialize the alignment table
    cursor = conn.cursor()
    s = (
            'create table if not exists alignment ('
            'offset integer, '
            'taxon integer, '
            'state integer, '
            'primary key (offset, taxon))')
    cursor.execute(s)
    conn.commit()

    # populate the primary state history table
    for t in alignment_triples:
        s = 'insert into alignment values (?, ?, ?)'
        cursor.execute(s, t)
    conn.commit()

    # close the alignment database table
    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--infile', required=True,
            help='input phylip file')
    parser.add_argument('--code', default='universal.code.db',
            help='output alignment patterns as an sqlite3 database file')
    parser.add_argument('--tree', default='p53S.tree.db',
            help='input tree for leaf taxon matching')
    parser.add_argument('--outfile', default='alignment.db',
            help='output codon alignment as an sqlite3 database file')
    main(parser.parse_args())

