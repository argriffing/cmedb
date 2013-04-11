"""
Construct a fast blinking process limit of an existing rate matrix.
"""

#XXX copypasted and under construction

import argparse
import sqlite3

import numpy as np
import scipy.linalg

import cmedbutil



def main(args):

    # construct and validate the mutational process equilibrium distribution
    nt_distn = np.array([args.A, args.C, args.G, args.T])
    cmedbutil.assert_stochastic_vector(nt_distn)
    nt_prob_map = {
            'A' : args.A,
            'C' : args.C,
            'G' : args.G,
            'T' : args.T,
            }

    # read the genetic code
    conn = sqlite3.connect(args.code)
    cursor = conn.cursor()
    cursor.execute(
        "select state, residue, codon from code "
        "where residue <> 'Stop' "
        "order by state")
    genetic_code = list(cursor)
    conn.close()

    # construct the mg94 rate matrix
    nstates = len(genetic_code)
    transitions = ('AG', 'GA', 'CT', 'TC')
    pre_Q = np.zeros((nstates, nstates), dtype=float)
    for a, (state_a, residue_a, codon_a) in enumerate(genetic_code):
        for b, (state_b, residue_b, codon_b) in enumerate(genetic_code):
            if hamming_distance(codon_a, codon_b) != 1:
                continue
            for nta, ntb in zip(codon_a, codon_b):
                if nta != ntb:
                    rate = nt_prob_map[ntb]
                    if nta + ntb in transitions:
                        rate *= args.kappa
            if residue_a != residue_b:
                rate *= args.omega
            pre_Q[a, b] = rate
    Q = pre_Q - np.diag(np.sum(pre_Q, axis=1))
    pre_distn = np.empty(nstates, dtype=float)
    for i, (state, residue, codon) in enumerate(genetic_code):
        pre_distn[i] = np.prod([nt_prob_map[nt] for nt in codon])
    distn = pre_distn / np.sum(pre_distn)

    # compute the expected syn and nonsyn rates for rescaling
    expected_syn_rate = 0.0
    expected_nonsyn_rate = 0.0
    for a, (state_a, residue_a, codon_a) in enumerate(genetic_code):
        for b, (state_b, residue_b, codon_b) in enumerate(genetic_code):
            if hamming_distance(codon_a, codon_b) != 1:
                continue
            rate = distn[a] * Q[a, b]
            if residue_a == residue_b:
                expected_syn_rate += rate
            else:
                expected_nonsyn_rate += rate

    # rescale the rate matrix to taste
    if args.expected_rate is not None:
        expected_rate = expected_syn_rate + expected_nonsyn_rate
        Q *= args.expected_rate / expected_rate
    elif args.expected_syn_rate is not None:
        Q *= args.expected_syn_rate / expected_syn_rate
    else:
        raise Exception

    # check time-reversible rate matrix invariants
    cmedbutil.assert_stochastic_vector(distn)
    cmedbutil.assert_rate_matrix(Q)
    cmedbutil.assert_equilibrium(Q, distn)
    cmedbutil.assert_detailed_balance(Q, distn)

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

    # populate the rate matrix table
    for a, (state_a, residue_a, codon_a) in enumerate(genetic_code):
        for b, (state_b, residue_b, codon_b) in enumerate(genetic_code):
            if a == b:
                rate = 0
            else:
                rate = Q[a, b]
            if rate > 0:
                cursor.execute(
                        'insert into rates values (?, ?, ?)',
                        (state_a, state_b, rate),
                        )
    conn.commit()

    # populate the states table
    for state, residue, codon in genetic_code:
        cursor.execute('insert into states values (?, ?)', (state, codon))
    conn.commit()

    # populate the state distribution table
    for prob, (state, residue, codon) in zip(distn, genetic_code):
        cursor.execute('insert into distn values (?, ?)', (state, prob))
    conn.commit()

    # close the database connection
    conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--proportion-on', type=cmedbutil.nonneg_float,
            help='partition part toleration probability')
    parser.add_argument('--syn-proportion-on', type=cmedbutil.nonneg_float,
            help='synonymous substitution toleration probability')
    parser.add_argument('--outfile', default='fast.blink.rate.matrix.db',
            help='output rate matrix in sqlite3 format')
    main(parser.parse_args())

