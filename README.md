# RRBS SCP

## Introduction

This guide provides steps for performing alignment and differential methylation analysis of MspI digested RRBS reads.

  * RRBS alignment notes
  * Sample methylation violin plot
  * Sample methylation dendrogram
  * Sample coverage
  * Fragment methylation calculation

### Requirements

  * Unix environment
  * Python 2.6+
  * R

### Software/Packages

  * Bowtie2 (<a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml"target="_blank">http://bowtie-bio.sourceforge.net/bowtie2/index.shtml</a>)
  * Pysam (<a href="https://code.google.com/archive/p/pysam/" target="_blank">https://code.google.com/archive/p/pysam/</a>)
  * BSSeeker2 (<a href="https://github.com/BSSeeker/BSseeker2" target="_blank">https://github.com/BSSeeker/BSseeker2</a>)
  * ggplot2 (<a href="http://ggplot2.org/" target="_blank">http://ggplot2.org/</a>)

## RRBS Alignment