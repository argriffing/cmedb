
# relative paths to scripts and data
myscripts = ../..
mycancer = ..

# parameter values estimated by codeml maximum likelihood
kappa_mle = 3.17632
omega_mle = 0.21925
T_mle = 0.18883
C_mle = 0.30126
A_mle = 0.25039
G_mle = 0.25952


#TODO under construction

all: likelihood.txt

tree.db:
	python $(myscripts)/newick-to-db.py \
		--newick=$(mycancer)/codeml.estimated.tree \
		--outfile=tree.db

codon.partition.db:
	python $(myscripts)/bp-codon-state-partition.py \
		--code=$(myscripts)/universal.code.db \
		--outfile=codon.partition.db

rate.matrix.db:
	python $(myscripts)/create-mg94.py \
		--code=$(myscripts)/universal.code.db \
		--omega=$(omega_mle) \
		--kappa=$(kappa_mle) \
		-A $(A_mle) \
		-C $(C_mle) \
		-G $(G_mle) \
		-T $(T_mle) \
		--expected-rate=1.0 \
		--outfile=rate.matrix.db

alignment.db: tree.db
	python $(mycancer)/create-p53-alignment.py \
		--code=$(myscripts)/universal.code.db \
		--infile=$(mycancer)/testseq \
		--tree=tree.db \
		--outfile=alignment.db

patterns.db: alignment.db
	python $(myscripts)/alignment-to-patterns.py \
		--alignment=alignment.db \
		--outfile=patterns.db

likelihood.txt: patterns.db tree.db rate.matrix.db
	python $(myscripts)/felsenstein-likelihood.py \
		--tree=tree.db \
		--rates=rate.matrix.db \
		--leaf-patterns=patterns.db \
		--method=felsenstein \
		> likelihood.txt


# wait I don't need this here...
tumor.db:
	$(mycancer)/create-disease-data.py \
		--code=$(myscripts)/universal.code.db \
		--infile=$(mycancer)/UMDTP53_curated_2012_R1_US.txt \
		--outfile=tumor.db


clean:
	rm *.db
	rm *.txt

