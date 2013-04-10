"""
Estimate disease vs. non-disease amino acids for each column.

Of all amino acids that can be reached by changing a single
nucleotide of the wild type codon,
count the number of amino acids that are present in the
database as putatively disease-causing
and count the number of amino acids that are not present in the
database and are therefore assumed to be benign.
"""

import argparse
import sqlite3
import math
from collections import defaultdict

def gen_codon_neighbors(codon):
    old_letters = list(codon)
    for i in range(3):
        for nt in 'ACGT':
            new_letters = list(old_letters)
            new_letters[i] = nt
            new_codon = ''.join(new_letters)
            if new_codon != codon:
                yield new_codon


def main(args):

    # read the genetic code
    conn = sqlite3.connect(args.code)
    cursor = conn.cursor()
    cursor.execute(
        "select state, residue, codon from code "
        "where residue <> 'Stop' "
        "order by state")
    genetic_code = list(cursor)
    conn.close()

    # read the disease data
    conn = sqlite3.connect(args.tumor)
    cursor = conn.cursor()
    cursor.execute('select offset, wild_state, mutant_state from tumor')
    disease_data = list(cursor)
    conn.close()

    # map from offset to wild and mutant states
    offsets = set()
    offset_to_wild_states = defaultdict(set)
    offset_to_mutant_states = defaultdict(set)
    for offset, wild_state, mutant_state in disease_data:
        offset_to_wild_states[offset].add(wild_state)
        offset_to_mutant_states[offset].add(mutant_state)
        offsets.add(offset)

    # verify that we only have one wild state per offset
    for item in offset_to_wild_states.items():
        offset, wild_states = item
        if len(wild_states) != 1:
            raise Exception(item)

    # compute maps related to the genetic code
    state_to_codon = {}
    codon_to_residue = {}
    for state, residue, codon in genetic_code:
        state_to_codon[state] = codon
        codon_to_residue[codon] = residue

    # count the amino acids
    informative_disease_aa_count = 0
    informative_benign_aa_count = 0
    for offset in offsets:
        wild_state = list(offset_to_wild_states[offset])[0]
        wild_codon = state_to_codon[wild_state]
        wild_residue = codon_to_residue[wild_codon]
        residue_neighbors = set()
        for neighbor in gen_codon_neighbors(wild_codon):
            residue = codon_to_residue.get(neighbor, None)
            if residue:
                residue_neighbors.add(residue)
        residue_neighbors -= set([wild_residue])
        observed_residues = set()
        for mutant_state in offset_to_mutant_states[offset]:
            mutant_codon = state_to_codon[mutant_state]
            mutant_residue = codon_to_residue[mutant_codon]
            observed_residues.add(mutant_residue)
        informative_disease_aas = residue_neighbors & observed_residues
        informative_benign_aas = residue_neighbors - observed_residues
        informative_disease_aa_count += len(informative_disease_aas)
        informative_benign_aa_count += len(informative_benign_aas)

    # print the counts
    print 'informative disease amino acid count:', informative_disease_aa_count
    print 'informative benign amino acid count:', informative_benign_aa_count

    # print the max log likelihood
    a = informative_disease_aa_count
    b = informative_benign_aa_count
    p = a / float(a + b)
    ll = a*math.log(p) + b*math.log(1-p)
    print 'max log likelihood:', ll


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tumor', default='tumor.db',
            help='input p53 tumor data as an sqlite3 database file')
    parser.add_argument('--code', default='universal.code.db',
            help='genetic code as an sqlite3 database file')
    main(parser.parse_args())

