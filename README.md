
# Founder Inference Bake-off #

To infer founder sequence(s) and time of infection:

	bin/infer.sh <sequences.fa>

Example usage:

	$ ~/src/matsen/founder-inference-bakeoff/bin/infer.sh  ~/src/matsen/hiv-sim/sample_data/caprisa002/gp120/caprisa002_pda_gp120_1m6m_aln.fasta 
	> caprisa002_pda_gp120_1m6m_aln_founder_1
	GCTCCAGCTGGTTATGCGATTCTAAGGTGTAATAATAAGACATTCAATGGGACAGGACCATGCAACAATGTCAGCACAGTACAATGTACACATGGAATTAAGCCAGTGGTATCAACTCAACTACTGTTAAATGGTAGCCTAGCAGAAGAGGAGATAATAATTAGATCTGAAAATCTGACAAGCAATCACAAAACAATAATAGTACAGCTTAATAGATCCATAGAAATTGTGTGCATAAGACCCGGCAATAACACAAGAAAAAGTGTAAGGATAGGACCAGGACAAACATTCTATGCAACAGGTGACATAATAGGAGACATAAGAAAAGCATATTGTAACATTAGTGCAGAAAGATGGAATGAAACTTTAGAATGGGTAAAGGAAAAATTAGCAGAACACTTTCCTAATAAGACAATAAGATATAAGCCATCTTCAGGAGGGGACCCAGAAGTTACAATGCATAGCTTT
	> caprisa002_pda_gp120_1m6m_aln_founder_2
	GCTCCAGCTGGTTATGCGATTCTAAGGTGTAATAATAAGACATTCAATGGGACAGGACCATGCAACAATGTCAGCACAGTACAATGTACACATGGAATTAAGCCAGTGGTATCAACTCAACTACTGTTAAATGGTAGCCTAGCAGAAGAGGAGATAATAATTAGATCTGAAAATCTGACAAACAATCACAAAACAATAATAGTACAGCTTAATAGATCCATAGAAATTGTGTGCATAAGACCCGGCAATAACACAAGAAAAAGTGTAAGGATAGGACCAGGACAAACATTCTATGCAACAGGTGACATAATAGGAGACATAAGAAAAGCATATTGTAACATTAGTGCAGAAAGATGGAATGAAACTTTAGAATGGGTAAAGAAAAAATTGGCAGAACACTTTCCTAATAAGACAATAAGATATCAACCATCTTCAGGAGGGGACCCAGAAGTTACAATGCATAGCTTT
	Estimated date of infection: 2005/06/19  (95% credible interval [2005/05/23 - 2005/07/15])


## Methodology ##

Uses Prank<sup>1</sup> to create a multiple sequence alignment and to estimate founder sequence.

Uses Poisson goodness-of-fit test<sup>2,3</sup> to evaluate whether the sample is
consistent with a single founder hypothesis.  Within a sample, only
the earliest sequences are used to detect multiplicity.

If multiple founders seem likely, the phylogenetic tree inferred by
Prank is split at the root and the leaves on each half are used to
infer separate founders.

Time of infection is inferred with Beast<sup>4</sup>.  A strict clock and constant
population size are used a priors.  Tree height is constrained to be
consistent with samples dates on the sequences.

--------
<sup>1</sup>: Löytynoja, A., & Goldman, N. (2008). A model of evolution and structure for multiple sequence alignment. Philosophical Transactions of the Royal Society of London. Series B, Biological Sciences, 363(1512), 3913–3919. http://doi.org/10.1098/rstb.2008.0170
<br/>
<sup>2</sup>: Giorgi EE, Funkhouser B, Athreya G, Perelson AS, Korber BT, Bhattacharya T. Estimating time since infection in early homogeneous HIV-1 samples using a Poisson model. BMC Bioinformatics 2010 Oct 25;11:532. [PMID: 20973976](http://www.ncbi.nlm.nih.gov/pubmed/20973976)
<br/>
<sup>3</sup>: Giorgi EE and Bhattacharya T. A note on two-sample tests for comparing intra-individual genetic sequence diversity between populations. Biometrics December 2012; 68:4. [PMID: 23004569](http://www.ncbi.nlm.nih.gov/pubmed/23004569)
<br/>
<sup>4</sup>: Drummond, Alexei J., and Andrew Rambaut. "BEAST: Bayesian evolutionary analysis by sampling trees." BMC evolutionary biology 7.1 (2007): 214.