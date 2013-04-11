cmedb
=====

~~vaporware~~ computational molecular evolution command line suite
using small sqlite databases for input and output


jargon used within the code and for script filenames
----------------------------------------------------

* mclr : Monte Carlo likelihood ratio
* bp : 'blinking process'
* ec : endpoint-conditioned (or leaf-conditioned for trees)
* ctmc : continuous-time Markov chain (often required to be time-reversible)


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

single sampled substitution history on a tree, not endpoint-conditioned:

    $ sqlite3 history.db
    SQLite version 3.7.13 2012-06-11 02:05:22
    Enter ".help" for instructions
    Enter SQL statements terminated with a ";"
    sqlite> .schema
    CREATE TABLE history (segment integer, va integer, vb integer, blen real, state integer, primary key (segment));
    sqlite> 

multiple sampled substitution histories on a tree, not endpoint-conditioned:

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

genetic code:

    $ sqlite3 universal.code.db
    SQLite version 3.7.13 2012-06-11 02:05:22
    Enter ".help" for instructions
    Enter SQL statements terminated with a ";"
    sqlite> .schema
    CREATE TABLE code (state integer, residue text, codon text, primary key (state));
    sqlite> 
    
a single path (not tree) endpoint-conditioned history sample:

    $ sqlite3 pathsample.db
    SQLite version 3.7.13 2012-06-11 02:05:22
    Enter ".help" for instructions
    Enter SQL statements terminated with a ";"
    sqlite> .schema
    CREATE TABLE history (segment integer, state integer, blen real, primary key (segment));
    sqlite> 

multiple path (not tree) endpoint-conditioned history samples:

    $ sqlite3 pathsamples.db 
    SQLite version 3.7.13 2012-06-11 02:05:22
    Enter ".help" for instructions
    Enter SQL statements terminated with a ";"
    sqlite> .schema
    CREATE TABLE histories (history integer, segment integer, state integer, blen real, primary key (history, segment));
    sqlite> 

endpoint-conditioned expectations from a time-reversible rate matrix:

    $ sqlite3 expectations.db                   
    SQLite version 3.7.13 2012-06-11 02:05:22
    Enter ".help" for instructions
    Enter SQL statements terminated with a ";"
    sqlite> .schema
    CREATE TABLE usage (initial integer, final integer, source integer, sink integer, usage real, primary key (initial, final, source, sink));
    CREATE TABLE wait (initial integer, final integer, state integer, wait real, primary key (initial, final, state));
    sqlite> 

an alignment sampled at nodes related by an unrooted tree,
using a time-reversible continuous Markov process:

    $ sqlite3 sampled.alignment.db 
    SQLite version 3.7.13 2012-06-11 02:05:22
    Enter ".help" for instructions
    Enter SQL statements terminated with a ";"
    sqlite> .schema
    CREATE TABLE alignment (offset integer, taxon integer, state integer, primary key (offset, taxon));
    sqlite> 

site patterns, with multiplicities, extracted from an alignment:

    $ sqlite3 patterns.db 
    SQLite version 3.7.13 2012-06-11 02:05:22
    Enter ".help" for instructions
    Enter SQL statements terminated with a ";"
    sqlite> .schema
    CREATE TABLE multiplicities (pattern integer, multiplicity integer, primary key (pattern));
    CREATE TABLE patterns (pattern integer, taxon integer, state integer, primary key (pattern, taxon));
    sqlite> 

full history sample of a complicated reversible Markov process on a path:

    $ sqlite3 blink.path.db                     
    SQLite version 3.7.13 2012-06-11 02:05:22
    Enter ".help" for instructions
    Enter SQL statements terminated with a ";"
    sqlite> .schema
    CREATE TABLE blink_history (part integer, segment integer, state integer, duration real, primary key (part, segment));
    CREATE TABLE primary_history (segment integer, state integer, duration real, primary key (segment));
    sqlite> 

log likelihoods under a model, for multiple path histories:

    $ sqlite3 log.likelihoods.db 
    SQLite version 3.7.13 2012-06-11 02:05:22
    Enter ".help" for instructions
    Enter SQL statements terminated with a ";"
    sqlite> .schema
    CREATE TABLE log_likelihoods (history integer, log_likelihood real, primary key (history));
    sqlite> 

p53 tumor data:

    $ echo ".schema" | sqlite3 cancer/tumor.db 
    CREATE TABLE tumor (offset integer, wild_state integer, mutant_state integer, primary key (offset, wild_state, mutant_state));

