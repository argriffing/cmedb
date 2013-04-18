
myscripts = ..

nstates = 8
nparts = 4
nsegments = 6
rate_on = 1
rate_off = 3


# Construct an arbitrary reversible continuous-time Markov process.
#
# Construct an arbitrary partition of the states.
#
# Define an interval of known length
# with arbitrary initial and final observed states.
#
# Define a rate-on and a rate-off for the blinking process.
#
# Sample an endpoint-conditioned path history.
#
# Compute the blinking process likelihood for the sampled path.
#
# Pick a root arbitrarily somewhere on the path,
# but not exactly at the time of a substitution.
#
# Construct a tree history from the sampled path history,
# with the arbitrary root.
#
# Compute the blinking process likeihood for the tree history.
#
# Check that the blinking process likelihoods are the same
# for the path history and for the tree history.


all: toy.tree.ll.db toy.tree.bp.ll.db path.ll.txt path.bp.ll.txt

toy.rate.matrix.db:
	python $(myscripts)/create-random-rate-matrix.py \
		--nstates=$(nstates) \
		--outfile=toy.rate.matrix.db

toy.partition.db: toy.rate.matrix.db
	python $(myscripts)/bp-random-state-partition.py \
		--rates=toy.rate.matrix.db \
		--nparts=$(nparts) \
		--outfile=toy.partition.db

toy.path.history.db toy.tree.history.db:
	python $(myscripts)/sample-toy-path-history.py \
		--initial=0 --final=1 --blen=3.2 \
		--nstates=$(nstates) --nsegments=$(nsegments) \
		--path-history=toy.path.history.db \
		--tree-history=toy.tree.history.db


# compute the non-blinking process log likelihoods

path.ll.txt: toy.rate.matrix.db toy.path.history.db
	python $(myscripts)/path-history-likelihood.py \
		--rates=toy.rate.matrix.db \
		--history=toy.path.history.db \
		> path.ll.txt

toy.tree.ll.db: toy.rate.matrix.db toy.tree.history.db
	python $(myscripts)/tree-histories-likelihoods.py \
		--rates=toy.rate.matrix.db \
		--histories=toy.tree.history.db \
		--outfile=toy.tree.ll.db


# compute the blinking process log likelihoods

path.bp.ll.txt: toy.rate.matrix.db toy.partition.db toy.path.history.db
	python $(myscripts)/path-history-blink-likelihood.py \
		--rates=toy.rate.matrix.db \
		--partition=toy.partition.db \
		--path-history=toy.path.history.db \
		--method=dynamic \
		--rate-on=$(rate_on) --rate-off=$(rate_off) \
		> path.bp.ll.txt

toy.tree.bp.ll.db: toy.rate.matrix.db toy.partition.db toy.tree.history.db
	python $(myscripts)/tree-histories-blink-likelihoods.py \
		--rates=toy.rate.matrix.db \
		--partition=toy.partition.db \
		--histories=toy.tree.history.db \
		--rate-on=$(rate_on) --rate-off=$(rate_off) \
		--outfile=toy.tree.bp.ll.db


clean:
	rm *.db
	rm *.txt

