#!/usr/bin/env Rscript

# plot distances form founder infered by prank to actual founder.
#
# Works of files created by a simulation process to mimic some of
# those that we obtain from the RV217 project.  The RV217 sample ahave
# around 10-12 sequences per sample, all taken from the c2v3c3 region
# gp-120 from the HIV-1 genome.

# install any missing packages
# http://stackoverflow.com/a/19873732/1135316
if (!suppressMessages(require("pacman"))) install.packages("pacman")
pacman::p_load(ape, dplyr, ggplot2, stringr, Biostrings, readr, tidyr, seqinr, DECIPHER, functional)

# turn off annoying progress bar
options(dplyr.show_progress=FALSE)

# cribbed from http://manuals.bioinformatics.ucr.edu/home/ht-seq#TOC-Multiple-Sequence-Alignments-MSAs-
removeGaps <- function(align, gap="-") {
    myseq <- maskMotif(align, gap)
    myseq <- paste(as(myseq, "XStringViews"), collapse="")
    return(myseq)
}

# Extract the identifier of the root node from a tree.
# Typically used to extract the name of the root node from a phylogeny inferred by PRANK.
root.id <- function(file) {
    tree <- ape::read.tree(file, comment.char = "")
    return(tree$node.label[1])
}

# Extract the inferred sequence at the root node of PRANK output
#
# 'dir' parameter is a character string naming a directory where the PRANK
# output files can be found.  The directory is expected to have a file
# 'prank.best.anc.dnd' containing a tree generated by PRANK, and a
# file 'prank.best.anc.fas' containing the sequences at each node of

# Read the founder sequences used to seed all the simulated RV!27 samples.
# If our simulation environment gets more complex we will need to
# parameterize this routine, but for now there is only one founder.
founder.sequence <- function(path=file.path('..', '..', 'templates', 'HIV1C2C3.fasta')) {
    fasta <- read.dna(path, format = "fasta", comment.char = "")
    return(paste(fasta[1,],collapse=''))
}



# Calculate the majority-rule consensus of sequences from the BEAST
# posterior traits after a burning period.
#
# NOTE This consensus is unrelated to the score used to assess the
# quality of the founder inferred by BEAST.
# That score is the mean of all pairwise alignments
# between sequences in the posterior and the actual founder, 
# NOT of the alignment between a single consensus and the
# founder.
beast.sequence <- function(dir, burnin=0.10) {
    # burning in the percentage of samples to skip at the beginning of the posterior.
    stopifnot(burnin >= 0 && burnin <= 1.0)

    consensus <- function(seqs) {
        seqinr::consensus(as.matrix(seqs), method = "majority", type="DNA")
    }
    tryCatch({
        read_tsv(file.path(dir, 'ancestralSequences.log'), col_names=c("state", "trait"), skip=3) %>%
        add_rownames() %>%
        filter(as.numeric(rowname) > (n() * burnin)) %>%
        summarize(trait=consensus(trait)) %>% as.character()
    }, error = function(e) {
        as.character(NA)
    })
}



# NOTE: MULTIPLE RETURN TYPES!
# if scoreOnly is FALSE, create a DNAbin object representing a pairwise alignment between two sequences
# otherwise just return the N-W alignment score.
pwalign <- function(s1, s2, scoreOnly=FALSE) {
    ## Nucleotide global, local, and overlap alignments
    pattern <- DNAString(gsub('-', '', s1))
    subject <- DNAString(gsub('-', '', s2))

    # First use a fixed substitution matrix
    mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = FALSE)
    globalAlign <- pairwiseAlignment(pattern, subject, type="global", scoreOnly=scoreOnly, substitutionMatrix = mat,
                                     gapOpening = 10, gapExtension = 4)

        ## # extracting the aligned sequences from a pairwiseAlignment object is unncessarily complex.
        ## # see https://support.bioconductor.org/p/70913/#70966

        ## seqset <- Biostrings:::.makePostalignedSeqs(globalAlign)[1]
        ## return(do.call(rbind, lapply(seqset, as.DNAbin)))
    return(globalAlign)
}

pwscore <- function(s1,s2) {
    tryCatch({pwalign(s1,s2,scoreOnly=TRUE)}, error = function(e) {NA})
}

# NOTE This distance calculation does not use the beast consensus sequence.
# Instead it is based on the mean score of all pairwise alignments
# between sequences in the posterior and the actual founder
beast.distance <- function(dir) {
    burnin_pct = 0.99
    n <- try({
        read_tsv(file.path(dir, 'ancestralSequences.log'), col_names=c("state", "trait"), skip=3) %>%
            add_rownames() %>%
            filter(as.numeric(rowname) > (n() * burnin_pct)) %>%
            select(c(trait)) %>%
            rowwise() %>% mutate(dist=pwscore(trait, founder)) %>%
            ungroup() %>% summarize(dist=mean(dist))
    }, error = function(e) {
        NA
    })
    return(as.numeric(n))
}

# calculate a majority-rule consensus from the input sequences.
# This is the simplest, most obvious founder inference method that works surprisingly well early in infection.
consensus.sequence <- function(dir, allowIUPAC=FALSE) {
    tryCatch({
        if (allowIUPAC) {
            readDNAStringSet(file.path(dir, 'sample_aln.fa')) %>%
                (function(s)  {chartr("-", "N", s)}) %>%
                consensusMatrix() %>%
                consensusString(ambiguityMap=IUPAC_CODE_MAP, threshold=0.25)
        } else {
            readDNAStringSet(file.path(dir, 'sample_aln.fa')) %>%
                as.matrix() %>% seqinr::consensus(method = "majority", type="DNA") %>%
                as.character() %>% paste0(collapse='')
        }
    }, error = function(e) {
        as.character(NA)
    })
}

# proportion of ambiguous positions in the consensus of all input sequences
nAmbig <- function(dir) {
    s <- consensus.sequence(dir, allowIUPAC=TRUE)
    t <- table(strsplit(toupper(s), '')[[1]])
    sum(t[!(names(t) %in% c('A','C','G','T'))])/nchar(s)
}

# proportion of sites in multiple alignment of all input sequences that are gaps
nGaps <- function(dir) {
    file.path(dir, 'sample_aln.fa') %>%
        read.dna(format = 'fasta', comment.char = '') %>%
        sapply(as.character) %>%
        (function(s) grepl('-',s)) %>%
        (function(b) sum(b)/length(b))
}

# total number of input sequences
nSeqs <- function(dir) {
    file.path(dir, 'sample_aln.fa') %>%
        read.dna(format = 'fasta', comment.char = '') %>%
        nrow()
}


# pairwiseGaps - total number of gaps in a pairwise alignment between a sequence and the an inferred prank-codon founder sequence
pairwiseGaps <- function(s, dir) {
    f <- prank.sequence(dir, prank.model='codon')
    n <- sum(grepl('-',sapply(pwalign(s, f), as.character)))
    return(n)
}

    

# Calculate the distance between the consensus of the INPUT sequences and the actual founder.
consensus.distance <- function(dir, model='K80', founder=founder) {
    # calculate a consensus of the input sequences
    tryCatch({
       consensus.sequence(dir) %>%
            sapply(function(s) removeGaps(DNAString(s))) %>%
            pwscore(founder)
    }, error = function(e) {
        as.numeric(NA)
    })
}


sample_seq <- function(dir, n=2) {
    readDNAStringSet(file.path(dir, 'sample_aln.fa')) %>%
        sample(n, replace=FALSE) %>%
        as.character() %>%
        (function(s) list(sname=names(s),seq=s)) %>% data.frame(stringsAsFactors=FALSE)
}


# Extract the inferred ancestral sequence at the root of the prank tree.
#
# This routine is parameterized by the name of the prank model, which
# is used to construct the name of the subdirectory when the prank
# results are expected.
prank.sequence <- function(dir, prank.model='codon') {
    tryCatch({
	dir <- file.path(dir, paste0('prank_', prank.model))
	treefile <- file.path(dir, 'prank.best.anc.dnd')
    	rootid <- root.id(treefile)
    	fastafile <- file.path(dir, 'prank.best.anc.fas')
    	fasta <- read.dna(fastafile, format = 'fasta', comment.char = '')

    	root <- fasta[which(labels(fasta)==rootid),]
	paste(root,collapse='')
    }, error = function(e) {
	NA
    })

}


prank.distance <- function(dir, prank.model='codon', model='K80', founder=founder) {
    #print(file)
    tryCatch({
        prank.sequence(dir, prank.model) %>% pwscore(founder)
    }, error = function(e) {prank.seq
        as.numeric(NA)
    })
}


get_sequences <- function(dir) {
    c(beast=beast.sequence(dir),
      pcodon=prank.sequence(dir,'codon'),
      pdna=prank.sequence(dir, 'dna'),
      unguided=prank.sequence(dir, 'unguided'),
      consensus=consensus.sequence(dir))
}

# save a multiple alignment to a file, and plot the pw alignment scores to the founder.
alignment <- function(dir) {
    DNAStringSet(get_sequences(dir)) %>%
        sapply(function(x) removeGaps(align=x)) %>%
        DNAStringSet() %>%
        AlignSeqs(verbose=FALSE) 
}

# pairwise alignment distances between the actual founder and each of the inferred founders in a directory
distances.dist.dna <- function(dir) {
    d <- alignment(dir) %>%
        as.matrix() %>% ape::as.alignment() %>% as.DNAbin() %>%
        dist.dna(as.matrix = TRUE, pairwise.deletion=TRUE)
    t(d['founder',2:ncol(d)])
}

# pairwise alignment scores between the actual founder and each of the inferred founders in a directory
distances.pwscore <- function(dir) {
    get_sequences(dir) %>% 
        sapply(function(s) pwscore(founder,s)) %>%
        sapply(function(s) -(s-nchar(founder))/nchar(founder)) %>%
        as.list() %>% as_data_frame() %>%
        do({names(.) <- paste0('dist.',names(.)); .})
}

# Measure the diversity of the aligned input sample.
# taken as the median of all the pairwise distance between individual sequences
diversity <- function(dir) {
    tryCatch({
        file.path(dir, 'sample_aln.fa') %>%
            readDNAStringSet() %>%
            as.matrix() %>% ape::as.alignment() %>% as.DNAbin() %>%
            dist.dna(pairwise.deletion=TRUE) %>% median() 
    }, error = function(e) {
        as.numeric(NA)
    })
}

## : the number of ambiguous bases in some consensus allowing ambiguity

# save a multiple alignment to a file, and plot the pw alignment scores to the founder.
write.alignment <- function(df) {
    print(class(df))
    print(dim(df))
          
    with(df, {
             print(dir)
             ## plotfile <- sprintf('plot.%s.%s.%s.replicate_%d.generation_%s.png',rate,fitness,ifelse(indel,"indel","noindel"),replicate,wpi)
             ## print(plotfile)
             ## png(file=plotfile, width=18, height=12, units = "in", res=300)
             ## print(barplot(data.frame(df)))
             ## dev.off()

             alnfile <- sprintf('alignment.%s.%s.%s.replicate_%d.generation_%s.fasta',rate,fitness,ifelse(indel,"indel","noindel"),replicate,wpi)
             aln <- alignment(dir)
             print(aln)
             writeXStringSet(aln, alnfile)
         })
}

# plot lines across the meanvalue of each inference methid, grouped by rate, fitness, indel, wpi, and method
lines_by_method <- function(tmp) {
    tmp %>%
        gather(method, dist, starts_with("dist.")) %>%
        group_by(rate, fitness, indel, wpi, method) %>%
        summarize(dist=mean(dist)) %>%
        mutate(method=sub('dist.','',method)) %>%
        rowwise() %>%
        mutate(method=sprintf("%s_%s", method,rate)) %>%
        ungroup() %>%
        mutate(method=factor(method)) %>%
        # filter(wpi==1200) %>%
        #filter(indel & fitness=='none') %>%
        ggplot(aes(x=wpi, y=dist, group=method, color=method)) +
        facet_grid(indel ~ fitness, labeller=facet_labeller) +
        geom_line() +
        #geom_point() +
        scale_color_discrete(name="Method",
                             breaks=c("beast", "consensus", 'pcodon', 'pdna', 'unguided'),
                             labels=c("Beast", "Consensus", 'Prank codon', 'Prank dna', 'Prank unguided'))
}




facet_labeller <- function(variable,value) {
    switch(variable,
           indel=ifelse(value, "With indels", "Without indels"),
           as.character(value))
}

scatterplot <- function(tmp) {
    tmp %>% gather(method, dist, starts_with("dist.")) %>%
        mutate(method=sub('dist.','',method)) %>%
        mutate(method=factor(method)) %>%
        filter(as.integer(as.character(wpi)) < 5000) %>%
        ggplot(aes(x=wpi, y=dist, group=method, color=method)) +
        facet_grid(indel ~ fitness, labeller=facet_labeller) +
        geom_point(position = position_jitter(w = 0.1)) +
        scale_color_discrete(name="Method",
                            breaks=c("beast", "consensus", 'pcodon', 'pdna', 'unguided'),
                            labels=c("Beast", "Consensus", 'Prank codon', 'Prank dna', 'Prank unguided')) 
        # geom_smooth(method="lm", fill=NA) 
}



boxplot <- function(tmp) {
    tmp %>% gather(method, dist, starts_with("dist.")) %>%
        mutate(method=sub('dist.','',method)) %>%
        mutate(method=factor(method)) %>%
        filter(as.integer(as.character(wpi)) < 20000) %>%
        mutate(wpi=factor(wpi)) %>%
        #filter(indel) %>%
        ggplot(aes(x=wpi, y=dist, fill=method)) +
        facet_grid(indel ~ fitness, labeller=facet_labeller) +
        theme(axis.text.x = element_text(size = 8, colour = "red", angle = 45)) +
        geom_boxplot() +
        xlab("Generations") +
        ylab("N-W Pairwise Distance") +

        scale_fill_discrete(name="Method",
                            breaks=c("beast", "consensus", 'pcodon', 'pdna', 'unguided', 'relaxed', 'strict'),
                            labels=c("Beast", "Consensus", 'Prank codon', 'Prank dna', 'Prank unguided', 'Beast relaxed', 'Beast strict')) 
}



barplot <- function(tmp) {
    tmp %>% gather(method, dist, starts_with("dist.")) %>%
        mutate(method=sub('dist.','',method)) %>%
        mutate(method=factor(method)) %>%
        filter(as.integer(as.character(wpi)) < 20000) %>%
        mutate(wpi=factor(wpi)) %>%
        #filter(indel) %>%
        ggplot(aes(x=factor(wpi), y=dist, fill=method)) +
        facet_grid(indel ~ fitness, labeller=facet_labeller) +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_text(aes(y=dist, ymax=dist, label=method), position= position_dodge(width=0.9), vjust=-.5, color="black") +
        scale_y_continuous("Distance") +
        scale_x_discrete("Generations") +
        scale_fill_discrete(name="Method",
                            breaks=c("beast", "consensus", 'pcodon', 'pdna', 'unguided'),
                            labels=c("Beast", "Consensus", 'Prank codon', 'Prank dna', 'Prank unguided')) 
}