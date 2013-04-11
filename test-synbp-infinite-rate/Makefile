
myscripts = ..
nsamples = 10

all: synbp.log.likelihoods.db limit.log.likelihoods.db

codon.partition.db:
	python $(myscripts)/bp-codon-state-partition.py \
		--code=$(myscripts)/universal.code.db \
		--outfile=codon.partition.db

limit.rate.matrix.db:
	python $(myscripts)/create-mg94.py \
		--code=$(myscripts)/universal.code.db \
		--omega=1.6 \
		--expected-syn-rate=0.25 \
		--outfile=limit.rate.matrix.db

plain.rate.matrix.db:
	python $(myscripts)/create-mg94.py \
		--code=$(myscripts)/universal.code.db \
		--omega=1 \
		--expected-syn-rate=1 \
		--outfile=plain.rate.matrix.db

histories.db: limit.rate.matrix.db
	python $(myscripts)/ctmc-segment-bridge-sampling.py \
		--initial=0 --final=1 --elapsed=5 \
		--table=histories --method=modified-rejection \
		--rates=limit.rate.matrix.db --nsamples=$(nsamples) \
		--outfile=histories.db

synbp.log.likelihoods.db: plain.rate.matrix.db histories.db codon.partition.db
	python $(myscripts)/path-histories-synbp-likelihoods.py \
		--rates=plain.rate.matrix.db --histories=histories.db \
		--partition=codon.partition.db \
		--rate-on=20000 --rate-off=30000 \
		--syn-rate-on=10000 --syn-rate-off=30000 \
		--outfile=synbp.log.likelihoods.db

limit.log.likelihoods.db: limit.rate.matrix.db histories.db codon.partition.db
	python $(myscripts)/path-histories-likelihoods.py \
		--rates=limit.rate.matrix.db --histories=histories.db \
		--outfile=limit.log.likelihoods.db

#
# Construct a toy mg94 rate matrix from a toy genetic code.
#

toy.partition.db:
	python $(myscripts)/create-random-rate-matrix.py \
		--outfile=codon.partition.db

limit.rate.matrix.db:
	python $(myscripts)/create-mg94.py \
		--code=$(myscripts)/universal.code.db \
		--omega=1.6 \
		--expected-syn-rate=0.25 \
		--outfile=limit.rate.matrix.db

limit.txt: limit.rate.matrix.db
	python $(myscripts)/ctmc-endpoint-probability.py \
		--initial=0 --final=1 --elapsed=5 \
		--rates=limit.rate.matrix.db \
		> limit.txt

synbp.txt: plain.rate.matrix.db toy.partition.db
	python $(myscripts)/synbp-endpoint-probability.py \
		--initial=0 --final=1 --elapsed=5 \
		--rates=plain.rate.matrix.db \
		--partition=partition.db \
		--rate-on=20000 --rate-off=30000 \
		--syn-rate-on=10000 --syn-rate-off=30000 \
		> synbp.txt

clean:
	rm *.db
	rm *.txt

