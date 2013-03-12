"""
This is a one-off script to create a tree database file.

It uses the brown.trees primate file from paml.
Each node of an unrooted phylogenetic tree
will be associated with an id number.
Each undirected edge of the tree will be associated with an id number.
"""

import sqlite3

def main():
    conn = sqlite3.connect('brown.tree.db')
    cursor = conn.cursor()

    # create the tables
    cursor.execute(
            'create table topo ('
            'edge integer, va integer, vb integer, '
            'primary key (edge))')
    cursor.execute(
            'create table blen ('
            'edge integer, blen real, '
            'primary key (edge))')
    cursor.execute(
            'create table taxa ('
            'node integer, name text, '
            'primary key (node))')
    conn.commit()

    # populate the tree topology table
    edge_va_vb_list = (
            (0, 0, 5),
            (1, 1, 5),
            (2, 2, 6),
            (3, 3, 6),
            (4, 4, 7),
            (5, 6, 7),
            (6, 5, 7),
            )
    for edge_va_vb in edge_va_vb_list:
        cursor.execute('insert into topo values (?, ?, ?)', edge_va_vb)
    conn.commit()

    # populate the branch length table
    edge_blen_list = (
            (0, 0.1),
            (1, 0.2),
            (2, 0.5),
            (3, 0.4),
            (4, 0.3),
            (5, 0.7),
            (6, 0.8),
            )
    for edge_blen in edge_blen_list:
        cursor.execute('insert into blen values (?, ?)', edge_blen)
    conn.commit()

    # populate the taxon name table
    node_name_list = (
            (0, 'human'),
            (1, 'chimpanzee'),
            (2, 'gibbon'),
            (3, 'orangutan'),
            (4, 'gorilla'),
            )
    for node_name in node_name_list:
        cursor.execute('insert into taxa values (?, ?)', node_name)
    conn.commit()

    # close the database connection
    conn.close()

if __name__ == '__main__':
    main()

