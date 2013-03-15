cmedb
=====

vaporware computational molecular evolution command line suite
using small sqlite databases for input and output


formats of various data types for computational molecular evolution
-------------------------------------------------------------------

tree format:

    $ sqlite3 brown.tree.db 
    SQLite version 3.7.13 2012-06-11 02:05:22
    Enter ".help" for instructions
    Enter SQL statements terminated with a ";"
    sqlite> .schema
    CREATE TABLE blen (edge integer, blen real, primary key (edge));
    CREATE TABLE taxa (node integer, name text, primary key (node));
    CREATE TABLE topo (edge integer, va integer, vb integer, primary key (edge));
    sqlite> 

rate matrix format:

    $ sqlite3 rate.matrix.db
    SQLite version 3.7.13 2012-06-11 02:05:22
    Enter ".help" for instructions
    Enter SQL statements terminated with a ";"
    sqlite> .schema
    CREATE TABLE distn (state integer, prob real, primary key (state));
    CREATE TABLE rates (source integer, sink integer, rate real, primary key (source, sink));
    CREATE TABLE states (state integer, name text, primary key (state));
    sqlite> 

sampled substitution history format:

    $ sqlite3 histories.db 
    SQLite version 3.7.13 2012-06-11 02:05:22
    Enter ".help" for instructions
    Enter SQL statements terminated with a ";"
    sqlite> .schema
    CREATE TABLE histories (history integer, segment integer, va integer, vb integer, blen real, state integer, primary key (history, segment));
    sqlite> 

