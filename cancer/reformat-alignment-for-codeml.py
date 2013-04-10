"""
Reformat the data because codeml is a picky eater.
"""

import argparse


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

    # read the phylip alignment
    infos = []
    with open(args.infile) as fin:
        infos = list(read_phylip(fin))

    # write the phylip alignment
    with open(args.outfile, 'w') as fout:
        first_taxon_name, first_codons = infos[0]
        print >> fout, len(infos), 3 * len(first_codons)
        print >> fout
        for taxon_name, codons in infos:
            print >> fout, '%s  %s' % (taxon_name, ''.join(codons))
            print >> fout


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', required=True,
            help='input phylip file')
    parser.add_argument('-o', '--outfile',
            default='alignment.for.codeml.phylip',
            help='output phylip file')
    main(parser.parse_args())

