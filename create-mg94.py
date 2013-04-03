"""
This script creates a time-reversible codon rate matrix in sqlite3 format.

The rate matrix is from Muse-Gaut (1994).
The output rate matrix is scaled to one expected substitution per time unit.
Genetic states corresponding to stop codons are not included
in the output equilibrium distribution or in the output rate matrix.
"""

import argparse
import sqlite3

import numpy as np
import scipy.linalg

import cmedbutil


#XXX unfinished

def get_MG_pre_Q(
        ts, tv, syn, nonsyn, asym_compo,
        nt_distn,
        kappa, omega,
        ):
    """
    This model is nested in FMutSel-F from which this code was copypasted.
    It was found in slowedml/codon1994.py
    """
    if nt_distn.shape != (4,):
        raise Exception(nt_distn.shape)
    A = (omega * nonsyn + syn) * (kappa * ts + tv)
    B = np.dot(asym_compo, nt_distn)
    pre_Q = A * B
    return pre_Q


def main(args):

    # construct and validate the mutational process equilibrium distribution
    ntdistn = np.array([args.A, args.C, args.G, args.T])
    cmedbutil.assert_stochastic_vector(ntdistn)

    # create or open the rate matrix database for writing
    conn = sqlite3.connect(args.outfile)
    cursor = conn.cursor()

    # create the tables
    cursor.execute(
            'create table rates ('
            'source integer, sink integer, rate real, '
            'primary key (source, sink))')
    cursor.execute(
            'create table states ('
            'state integer, name text, '
            'primary key (state))')
    cursor.execute(
            'create table distn ('
            'state integer, prob real, '
            'primary key (state))')
    conn.commit()

    # define the nucleotide states
    nt_states = 'ACGT'

    # populate the rate matrix table
    for i, si in enumerate(states):
        for j, sj in enumerate(states):
            if i == j:
                rate = 0
            elif si+sj in ('AG', 'GA', 'CT', 'TC'):
                rate = args.kappa * distn[j]
            else:
                rate = distn[j]
            if rate > 0:
                cursor.execute(
                        'insert into rates values (?, ?, ?)',
                        (i, j, rate),
                        )
    conn.commit()

    # populate the states table
    for pair in enumerate(states):
        cursor.execute('insert into states values (?, ?)', pair)
    conn.commit()

    # populate the state distribution table
    for pair in enumerate(distn):
        cursor.execute('insert into distn values (?, ?)', pair)
    conn.commit()

    # close the database connection
    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-A', default=0.25, type=cmedbutil.pos_float,
            help='mutational equilibrium probability of nucleotide A')
    parser.add_argument('-C', default=0.25, type=cmedbutil.pos_float,
            help='mutational equilibrium probability of nucleotide C')
    parser.add_argument('-G', default=0.25, type=cmedbutil.pos_float,
            help='mutational equilibrium probability of nucleotide G')
    parser.add_argument('-T', default=0.25, type=cmedbutil.pos_float,
            help='mutational equilibrium probability of nucleotide T')
    parser.add_argument('-k', '--kappa', default=2.0, 
            type=cmedbutil.pos_float,
            help='transition/transversion ratio')
    parser.add_argument('-w', '--omega', default=0.05,
            type=cmedbutil.pos_float,
            help='nonsynonymous/synonymous ratio')
    parser.add_argument('--code', default='universal.code.db',
            help='input genetic code in sqlite3 format')
    parser.add_argument('--outfile', default='rate.matrix.db',
            help='output rate matrix in sqlite3 format')
    main(parser.parse_args())

