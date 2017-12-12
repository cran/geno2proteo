\name{geno2proteo-package}
\alias{geno2proteo-package}
\alias{geno2proteo}
\docType{package}
\title{
\packageTitle{geno2proteo}
}
\description{
\packageDescription{geno2proteo}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{geno2proteo}
\packageIndices{geno2proteo}


~~ An overview of how to use the package and the most important functions ~~

The package needs three data files. One contains the genetic table, another 
one contains the DNA sequences of a certain genome, and the third one contains
the gene annotations in the same format as the GTF file used in ENSEMBL.
The standard genetic table was provided in this package. The other two data 
files can be downloaded from ENSEMBL web site. For details about those data 
files and how to obtain them, see the introduction document of this package. 
Of course you can create your own data files, as long as you observe the format 
of the files which are specified in the introduction document as well. 

The package also needs Perl installed and being available for use in order to 
run some of the functions. 

Once you have the three data files, you need to run the function 
\code{generatingCDSaaFile} to generate another data file, which 
will be used by some of the functions in this package.

Four main functions were implemented in this package so far.
The function \code{genomicLocToProteinSequence} will find the protein 
sequences and DNA sequences of the coding regions within a list of genomic 
loci given as input. \code{genomicLocToWholeDNASequence} will obtain the 
whole DNA sequences of any genomic loci given in the input data. 
\code{proteinLocsToGenomic} will find the genomic coordinates for a list of
sections in proteins as input. \code{proteinLocsToProteinSeq} will find the
the protein sequences of a list of the protein sections.

}
\author{
\packageAuthor{geno2proteo}

Maintainer: \packageMaintainer{geno2proteo}
}

