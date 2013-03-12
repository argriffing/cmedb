cmedb
=====

vaporware computational molecular evolution command line suite
using small sqlite databases for input and output


schemas
-------

tree format:

    $ sqlite3 brown.tree.db 
    SQLite version 3.7.13 2012-06-11 02:05:22
    Enter ".help" for instructions
    Enter SQL statements terminated with a ";"
    sqlite> .schema
    CREATE TABLE blen (edge integer, blen real, primary key (edge));
    CREATE TABLE taxa (node integer, name text, primary key (node));
    CREATE TABLE topo (edge integer, va integer, vb integer, primary key (edge));


