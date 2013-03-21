cmedb
=====

vaporware computational molecular evolution command line suite
using small sqlite databases for input and output


formats of various data types for computational molecular evolution
-------------------------------------------------------------------

tree:

    $ sqlite3 brown.tree.db 
    SQLite version 3.7.13 2012-06-11 02:05:22
    Enter ".help" for instructions
    Enter SQL statements terminated with a ";"
    sqlite> .schema
    CREATE TABLE blen (edge integer, blen real, primary key (edge));
    CREATE TABLE taxa (node integer, name text, primary key (node));
    CREATE TABLE topo (edge integer, va integer, vb integer, primary key (edge));
    sqlite> 

rate matrix:

    $ sqlite3 rate.matrix.db
    SQLite version 3.7.13 2012-06-11 02:05:22
    Enter ".help" for instructions
    Enter SQL statements terminated with a ";"
    sqlite> .schema
    CREATE TABLE distn (state integer, prob real, primary key (state));
    CREATE TABLE rates (source integer, sink integer, rate real, primary key (source, sink));
    CREATE TABLE states (state integer, name text, primary key (state));
    sqlite> 

single sampled substitution history:

    $ sqlite3 history.db
    SQLite version 3.7.13 2012-06-11 02:05:22
    Enter ".help" for instructions
    Enter SQL statements terminated with a ";"
    sqlite> .schema
    CREATE TABLE history (segment integer, va integer, vb integer, blen real, state integer, primary key (segment));
    sqlite> 

multiple sampled substitution histories:

    $ sqlite3 histories.db 
    SQLite version 3.7.13 2012-06-11 02:05:22
    Enter ".help" for instructions
    Enter SQL statements terminated with a ";"
    sqlite> .schema
    CREATE TABLE histories (history integer, segment integer, va integer, vb integer, blen real, state integer, primary key (history, segment));
    sqlite> 

2d vertex layout:

    $ sqlite3 layout.db
    SQLite version 3.7.13 2012-06-11 02:05:22
    Enter ".help" for instructions
    Enter SQL statements terminated with a ";"
    sqlite> .schema
    CREATE TABLE layout (node integer, x real, y real, primary key (node));
    sqlite> 

segmentation of a sampled history into contiguous isostate regions:

    $ sqlite3 segmentation.db 
    SQLite version 3.7.13 2012-06-11 02:05:22
    Enter ".help" for instructions
    Enter SQL statements terminated with a ";"
    sqlite> .schema
    CREATE TABLE bipartite (isostate integer, parity integer, primary key (isostate));
    CREATE TABLE isostate (segment integer, isostate integer, primary key (segment));
    sqlite> 

