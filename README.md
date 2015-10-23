
# Founder Inference Bake-off #

To infer founder sequence(s) and time of infection:

    bin/infer.sh <sequences.fa>

## Methodology ##

Uses Prank[^1] to create a multiple sequence alignment and to estimate founder sequence.

Uses Poisson goodness-of-fit test to evaluate whether the sample is
consistent with a single founder hypothesis.  Within a sample, only
the earliest sequences are used to detect multiplicty.

If multiple founders seem likely, the phylogeneic tree inferred by
Prank is split at the root and the leaves on each half are used to
infer seperate founders.

Time of infection is inferred with Beast.  A strict clock and constant
population size are used a priors.  Tree height is constrained to be
consistent with samples dates on the sequences.

[^1]: Löytynoja, A., & Goldman, N. (2008). A model of evolution and structure for multiple sequence alignment. Philosophical Transactions of the Royal Society of London. Series B, Biological Sciences, 363(1512), 3913–3919. http://doi.org/10.1098/rstb.2008.0170
