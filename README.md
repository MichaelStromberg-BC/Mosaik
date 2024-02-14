![Mosaik](https://github.com/MichaelStromberg-BC/Mosaik/raw/main/mosaik.png)

During my Ph.D., this was my baby. This version represents the last version that was published before I graduated.

This aligner was historical for many reasons. During the 1000 Genomes Project, this aligner was heavily used for Illumina, SOLID, and 454 reads. It was also the first aligner in that project to produce gapped alignments (these days that's an accepted standard). In addition, it was the first aligner that used a neural network to produce well-behaved alignment/mapping qualities.

N.B. This is **NOT** the latest version. Wan-Ping Lee took over the Mosaik project after my departure and added a performant banded Smith-Waterman algorithm that made the aligner sing.

## Introduction

MOSAIK is a reference-guided assembler comprising of two main modular programs:

* MosaikBuild
* MosaikAligner

MosaikBuild converts various sequence formats into Mosaik’s native read format. MosaikAligner pairwise aligns each read to a specified series of reference sequences and produces BAMs as outputs.

At this time, the workflow consists of supplying sequences in FASTA, FASTQ, Illumina Bustard & Gerald, or SRF file formats and producing results in the BAM format. More information can be found in the [MOSAIK 1.0 Documentation](Mosaik-1.0-Documentation.pdf)

## What's new?

* A new neural-net for mapping quality (MQ) calibration is introduced. Initial testing using simulated reads shows that this method improve the accuracy compared to the previous MQ scheme.

* A local alignment search option has been added to help rescue mates in paired-end/mate-pair reads that may be missing due to highly repetitive regions in the genome.

* SOLiD support has finally come of age. MOSAIK imports and aligns SOLiD reads in colorspace, but now seamlessly converts the alignments back into basespace. No more downstream bioinformatics headaches.

* Robust support for the BAM alignment file formats.

* The command line parameters have been cleaned up and sensible default parameters have been chosen. This cuts down the ridiculously long command-lines to simply specifying an input file and an output file in most cases.

## What makes MOSAIK different?

Unlike many current read aligners, **MOSAIK** produces gapped alignments using the Smith-Waterman algorithm.

**MOSAIK** is written in highly portable C++ and currently targetted for the following platforms: Microsoft Windows, Apple Mac OS X, FreeBSD, and Linux operating systems. Other platforms can easily be supported upon request.

**MOSAIK** is multithreaded. If you have a machine with 8 processors, you can use all 8 processors to align reads faster while using the same memory footprint as when using one processor.

**MOSAIK** supports multiple sequencing technologies. In addition to legacy technologies such as Sanger capillary sequencing, our program supports next generation technologies such as Roche 454, Illumina, AB SOLiD, and experimental support for the Helicos Heliscope.

## Publication

Lee W-P, Stromberg MP, Ward A, Stewart C, Garrison EP, Marth GT (2014) MOSAIK: A Hash-Based Algorithm for Accurate Next-Generation Sequencing Short-Read Mapping. PLoS ONE 9(3): e90581. https://doi.org/10.1371/journal.pone.0090581
