
* the following codeml control file parameters are described
* in Khan et al. Maximum likelihood analysis of mammalian p53... (2011).

clock = 0
model = 0
cleandata = 0
fix_kappa = 0
fix_omega = 0


* other stuff added later

* user tree
runmode = 0

* get the standard errors of estimates
getSE = 1

* use codon sequences
seqtype = 1

* use a muse-gaut model
CodonFreq = 4
estFreq = 1

* input files
seqfile = alignment.for.codeml.phylip
outfile = mlb
treefile = p53S.tree
