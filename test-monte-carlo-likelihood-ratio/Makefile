
myscripts = ..
nstates = 5
nparts = 3
nsamples = 100000

all: mclr.txt ctmc.txt synbp.txt

toy.rate.matrix.db:
	python $(myscripts)/create-random-rate-matrix.py \
		--nstates=$(nstates) --outfile=toy.rate.matrix.db


toy.partition.db: toy.rate.matrix.db
	python $(myscripts)/bp-random-state-partition.py \
		--outfile=toy.partition.db --nparts=$(nparts) \
		--rates=toy.rate.matrix.db

toy.histories.db: toy.rate.matrix.db
	python $(myscripts)/ctmc-segment-bridge-sampling.py \
		--initial=0 --final=1 --method=modified-rejection \
		--elapsed=5 --rates=toy.rate.matrix.db --nsamples=$(nsamples) \
		--table=histories --outfile=toy.histories.db

toy.blink.log.likelihoods.db: toy.rate.matrix.db toy.histories.db \
	toy.partition.db
	python $(myscripts)/path-histories-synbp-likelihoods.py \
		--rates=toy.rate.matrix.db --histories=toy.histories.db \
		--rate-on=1 --rate-off=3 --partition=toy.partition.db \
		--syn-rate-on=2 --syn-rate-off=3 \
		--outfile=toy.blink.log.likelihoods.db

toy.reference.log.likelihoods.db: toy.rate.matrix.db toy.histories.db
	python $(myscripts)/path-histories-likelihoods.py \
		--rates=toy.rate.matrix.db --histories=toy.histories.db \
		--outfile=toy.reference.log.likelihoods.db

mclr.txt: toy.blink.log.likelihoods.db toy.reference.log.likelihoods.db
	python $(myscripts)/monte-carlo-likelihood-ratio.py \
		--numerator-log-likelihood=toy.blink.log.likelihoods.db \
		--denominator-log-likelihood=toy.reference.log.likelihoods.db \
		> mclr.txt

ctmc.txt: toy.rate.matrix.db
	python $(myscripts)/ctmc-endpoint-probability.py \
		--elapsed=5 --initial=0 --final=1 --rates=toy.rate.matrix.db \
		> ctmc.txt

synbp.txt: toy.rate.matrix.db toy.partition.db
	python $(myscripts)/synbp-endpoint-probability.py \
		--elapsed=5 --initial=0 --final=1 --rates=toy.rate.matrix.db \
		--partition=toy.partition.db --rate-on=1 --rate-off=3 \
		--syn-rate-on=2 --syn-rate-off=3 \
		> synbp.txt

clean:
	rm *.db
	rm *.txt

