"""
This is a one-off script to create a tree database file.
"""

from StringIO import StringIO
import sqlite3

import dendropy


g_newick_string = (
        '((((((Has,Ptr),Ppy),(((Mmu,Mfu),Mfa),Cae)),(Mim,Tgl)),'
        '((((((Mum,Rno),Mun),(Cgr,Mau)),Sju),(Cpo,Mmo)),(Ocu,Opr))),'
        '(Sar,((Fca,Cfa),((Bta,Oar),Dle))));')


def main():

    # use dendropy to read this newick file
    t = dendropy.Tree(stream=StringIO(g_newick_string), schema='newick')
    leaves = t.leaf_nodes()
    nodes = list(t.postorder_node_iter())
    non_leaves = [n for n in nodes if n not in leaves]
    ordered_nodes = leaves + non_leaves

    # node index lookup
    node_id_to_index = dict((id(n), i) for i, n in enumerate(ordered_nodes))

    # get edges
    edge_va_vb_list = []
    edges = list(t.postorder_edge_iter())
    for i, edge in enumerate(edges):
        if edge.head_node and edge.tail_node:
            va = node_id_to_index[id(edge.head_node)]
            vb = node_id_to_index[id(edge.tail_node)]
            edge_va_vb_list.append((i, va, vb))

    # get a list of (node, name) pairs for the table
    node_name_list = [(i, str(n.taxon)) for i, n in enumerate(leaves)]

    # open the database file
    conn = sqlite3.connect('p53S.tree.db')
    cursor = conn.cursor()

    # create the tables
    cursor.execute(
            'create table topo ('
            'edge integer, va integer, vb integer, '
            'primary key (edge))')
    cursor.execute(
            'create table taxa ('
            'node integer, name text, '
            'primary key (node))')
    conn.commit()

    # populate the tree topology table
    for edge_va_vb in edge_va_vb_list:
        cursor.execute('insert into topo values (?, ?, ?)', edge_va_vb)
    conn.commit()

    # populate the taxon name table
    for node_name in node_name_list:
        cursor.execute('insert into taxa values (?, ?)', node_name)
    conn.commit()

    # close the database connection
    conn.close()


if __name__ == '__main__':
    main()

