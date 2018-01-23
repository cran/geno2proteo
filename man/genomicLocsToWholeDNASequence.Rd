\name{genomicLocsToWholeDNASequence}
\alias{genomicLocsToWholeDNASequence}
\title{
Obtaining the DNA sequences of a list of genomic loci
}
\description{
The function takes a list of genomic loci and tries to find the whole DNA 
sequences within each of the loci.
}
\usage{
genomicLocsToWholeDNASequence(inputLoci, DNAfastaFile, 
                                tempFolder = "./", perlExec='perl')
}
\arguments{
    \item{inputLoci}{
A data frame containing the genomic loci as the input. Each row is for one 
genomic locus. The first column is the chromosome name, the 2nd and 3rd columns
are the start and end coordinates of the locus in the chromosome, and the 4th 
column specifies the strand of chromosome ("+" and "-" for forward and reverse
strand, respectively). Other columns are optional and will not be used by the 
function. Note that the chromosome name can be either in the ENSEMBL style, 
e.g. 1, 2, 3, \dots, and X, Y and MT, or in another popular style, namely 
chr1, chr2, chr3, \dots, and chrX, chrY and chrM. But they cannot be mixed in 
the input of one function call.
}
    \item{DNAfastaFile}{
The name of a fasta file containing the whole DNA sequence of the genome used.
For details about this data file see the documentation of this package.
}
    \item{tempFolder}{
A temporary folder into which the program can write some temporary files 
which will be deleted when the function running is finished. The default value
is the current folder.
}
    \item{perlExec}{
Its value should be the full path of the executable file which can be used 
to run Perl scripts (e.g. "/usr/bin/perl" in a linux computer or 
"C:/Strawberry/perl/bin/perl" in a Windows computer). The default value is 
"perl".
}
}
\details{
This function obtains the whole DNA sequences of a list of genomic loci. 
Note that, in contrast, another function \code{genomicLocToProteinSequence} in
this package can return the DNA sequences of the coding regions within the 
given genomic loci. 
}
\value{
The function returns a data frame containing the original genomic loci as in 
the input and after them, one additional column for the DNA sequence of the 
corresponding genomic locus.
}

\author{
Yaoyong Li
}

\examples{

    dataFolder = system.file("extdata", package="geno2proteo")
    inputFile_loci=file.path(dataFolder, 
        "transId_pfamDomainStartEnd_chr16_Zdomains_22examples_genomicPos.txt")
    DNAfastaFile =  file.path(dataFolder, 
        "Homo_sapiens.GRCh37.74.dna.chromosome.16.fa_theFirst3p5M.txt.gz")

    inputLoci = read.table(inputFile_loci, sep="\t", stringsAsFactors=FALSE)
    
    tmpFolder = tempdir()

    DNASeqNow = genomicLocsToWholeDNASequence(inputLoci=inputLoci, 
                            DNAfastaFile=DNAfastaFile, tempFolder=tmpFolder)

}

