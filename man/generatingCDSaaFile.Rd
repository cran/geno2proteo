\name{generatingCDSaaFile}
\alias{generatingCDSaaFile}
\title{
Generating a file containing the DNA and AA sequences of all the protein 
coding regions (CDSs) in a genome
}
\description{
This function will find the DNA and protein sequences for each CDS region 
listed in the ENSEMBL gene annotation file (gtf file) provided, and store the
CDS regions and the corresponding DNA and protein sequences in an output data 
file. The output data file will be needed by some of the functions in this 
package.   
}
\usage{
generatingCDSaaFile(geneticCodeFile_line, gtfFile, DNAfastaFile, 
outputFolder = "./", perlExec='perl')
}
\arguments{
    \item{geneticCodeFile_line}{
A text file containing a genetic coding table of the codons and the amino acids
coded by those codons. A file containing the standard genetic coding table is 
provided in this package, which is the file
"geneticCode\_standardTable\_lines.txt" 
in the folder "geno2proteo/extdata/" which will be located in the folder that
you install the package. This will be 
available after you have installed the package. 

Alternatively you can create your own genetic table if needed. The format of 
the genetic table used in this package 
is one line for one codon. The first column is the codon (namely 3 DNA letters)
and the second column 
is the name (in a single letter) of the amino acid coded by the codon. You may
have more columns 
if you want but only the first two columns are used by the package, and the
columns are separated by a tab.
}
    \item{gtfFile}{
A text file in GTF format containing the gene annotations of the species
that you are interested in. You may obtain this file from the ENSEMBL 
web site. The text file can be compressed by GNU Zip (e.g. gzip). For the 
details about how to get the 
data file from ENSEMBL web site, see the documentation of this package.
}
    \item{DNAfastaFile}{
A text file in fasta format containing the DNA sequence for the genome that 
you want to use 
for analysing the data. The text file can be compressed by GNU Zip (e.g. gzip).
You may download the file directly 
from the ENSEMBL web site. For the details about how to get the 
data file from ENSEMBL web site, see the documentation of this package.
}
    \item{outputFolder}{
    A folder where the result file will be stored. The default value is the 
    current folder "./".
}
    \item{perlExec}{
    Its value should be the full path of the executable file which can be used 
    to run Perl scripts (e.g. "/usr/bin/perl" in a linux computer
    or "C:/Strawberry/perl/bin/perl" in a Windows computer). 
    The default value is "perl".
}
}
\details{
This function generates a data file containing the genomic locations, 
DNA sequences and protein sequences of all of the coding regions (CDSs), which
will be used by some of other functions in this packages.
}
\value{
This function does not return any specific value but generates a file in the 
output folder containing the DNA and AA sequences of the CDS regions. The first
part of the name of the result file is same as the gene annotation file 
(namely the gtf file), and the last part is `\_AAseq.txt.gz'. 
}

\author{
Yaoyong Li
}

\examples{
    # the data folder in this package
    dataFolder = system.file("extdata", package="geno2proteo") 

    geneticCodeFile_line = file.path(dataFolder, 
                                    "geneticCode_standardTable_lines.txt")
    gtfFile = file.path(dataFolder, 
                "Homo_sapiens.GRCh37.74_chromosome16_35Mlong.gtf.gz")
    DNAfastaFile =  file.path(dataFolder, 
            "Homo_sapiens.GRCh37.74.dna.chromosome.16.fa_theFirst3p5M.txt.gz")
    
    outputFolder = "./"; # using the current folder as output folder
    # calling the function.
    generatingCDSaaFile(geneticCodeFile_line=geneticCodeFile_line, 
        gtfFile=gtfFile, DNAfastaFile=DNAfastaFile, outputFolder=outputFolder)
    
    filename00 = sub(".*/", "", gtfFile)
    # get the output file's name.
    outputFile = paste(outputFolder, "/", filename00, "_AAseq.txt.gz", sep="")
    # read the content of the output file into a data frame.
    aaSeq = read.table(outputFile, sep="\t", stringsAsFactors=FALSE)

}

