"""
Compute a blinking-process probability of primary state endpoints.

This is on a path, not a tree.
Only the initial and final primary states are observed.
The primary state is not observed along the path,
and the blinking process is not observed anywhere.
This script is slow, and it computes an exact probability
by summing over the uninknown blinking process initial state distribution
and integrating over all of the unknown blinking process transitions
and over all of the primary path transitions.
This is not feasible for actual genetic codes.
.
The practical purpose of this script is to work in tandem with
artificial genetic codes for the purposes of testing "Monte Carlo
likelihood ratio" (MCLR) implementations.
If Monte Carlo likelihood ratio estimation seems to converge to the
exact likelihood ratios for small artificial genetic codes,
then the MCLR implementation would presumably not be buggy and
could be used for inference with non-toy genetic codes.
"""

#XXX under construction

