test_generatingCDSaaFile <- function() {

    dataFolder = system.file("extdata", package="geno2proteo")

    geneticCodeFile_line = file.path(dataFolder, 
        "geneticCode_standardTable_lines.txt")
    gtfFile = file.path(dataFolder, 
        "Homo_sapiens.GRCh37.74_chromosome16_35Mlong.gtf.gz")
    DNAfastaFile =  file.path(dataFolder, 
        "Homo_sapiens.GRCh37.74.dna.chromosome.16.fa_theFirst3p5M.txt.gz");

    outputFolder = dataFolder;
    generatingCDSaaFile(geneticCodeFile_line=geneticCodeFile_line, 
        gtfFile=gtfFile, DNAfastaFile=DNAfastaFile, outputFolder=outputFolder)
    
    kk = sub('.*/', '', gtfFile)
    # get the output file name
    outputFile = paste(outputFolder, '/', kk, '_AAseq.txt.gz', sep='') 
    aaSeq = read.table(outputFile, sep='\t', stringsAsFactors=F)
    
    # the first exon of ENST00000215798
    #checkEquals(aaSeq[164,12], 'XVRLPRQDALVLEGVRIGSEADPAPLLGGRLLL') 
    #checkEquals(aaSeq[171, 12], 
    #'WLRVLPCKHEFHRDCVDPWLMLQQTCPLCKFNVLGEHQGWGPSAYSACSSPDASLPVLLPLPCREPLLR*')
    
    #checkEquals(aaSeq[2852,12], 'XWMDRGTRDEHLPSCPGCPAVASNTRCCPPAACLPGISLS')
    #checkEquals(aaSeq[2854,12], 'WVAGLLDPSMDK*') 
    
    checkEquals(aaSeq[1,12], 'MGARGALLLALLLARAGLRKP') # the first one
    checkEquals(substr(aaSeq[3,12],1, 20), 'GPCGRRVITSRIVGGEDAEL')
    
    # the last one.
    checkEquals(aaSeq[nrow(aaSeq),12], 'VLYFPDRWWHATLNLDTSVFISTFLG')
    
    checkEquals(aaSeq[nrow(aaSeq)-1,12], 
    'RWFLYPPEKTPEFHPNKTTLAWLRDTYPALPPSARPLECTIRAGE')
    
    #checkEqualsNumeric(divideBy(4, 1.2345), 3.24, tolerance=1.0e-4)
}

test_genomicLocToProteinSequence <- function() {

    dataFolder = system.file("extdata", package="geno2proteo")
    inputFile_loci=file.path(dataFolder, 
        "transId_pfamDomainStartEnd_chr16_Zdomains_22examples_genomicPos.txt")
    CDSaaFile=file.path(dataFolder, 
        "Homo_sapiens.GRCh37.74_chromosome16_35Mlong.gtf.gz_AAseq.txt.gz")

    inputLoci = read.table(inputFile_loci, sep='\t', 
        stringsAsFactors=F, header=T)
    proteinSeq = genomicLocToProteinSequence(inputLoci=inputLoci, 
        CDSaaFile=CDSaaFile)
  
    kkt = sapply(1:nrow(proteinSeq), function(n) {a = proteinSeq[n,]; 
        grepl(as.character(a[1,7]), as.character(a[1,14]) )}  )
    kkt1 = proteinSeq[kkt,]

    kk = as.character(kkt1[,'pepSeq'])

     # the first one
    checkEquals(kk[1], 'PVTFEDVALYLSREEWGRLDHTQQNFYRDVLQK')

    # the last one
    checkEquals(substr(kk[length(kk)], 1, 20), 'AGPVALGDIPFYFSREEWGT')

}

test_genomicLocToWholeDNASequence <- function() {

    dataFolder = system.file("extdata", package="geno2proteo")
    inputFile_loci=file.path(dataFolder, 
    "transId_pfamDomainStartEnd_chr16_Zdomains_22examples_genomicPos.txt")
    DNAfastaFile =  file.path(dataFolder, 
    "Homo_sapiens.GRCh37.74.dna.chromosome.16.fa_theFirst3p5M.txt.gz");

    inputLoci = read.table(inputFile_loci, sep='\t', stringsAsFactors=F, 
        header=T)

    DNASeqNow = genomicLocToWholeDNASequence(inputLoci=inputLoci, 
        DNAfastaFile=DNAfastaFile)

    kk = as.character(DNASeqNow[,ncol(DNASeqNow)])

    checkEquals(substr(kk[1], 50, 70), 'GACGGCTGGACCACACGCAGC')

    checkEquals(kk[10], 
    'TACCACTGCCTCGACTGCGGCAAGAGCTTCAGCCACAGCTCGCACCTCACCGCGCACCAGCGCACCCAC')

}

test_proteinLocsToGenomic <- function() {

    dataFolder = system.file("extdata", package="geno2proteo")
    inputFile_loci=file.path(dataFolder, 
        'transId_pfamDomainStartEnd_chr16_Zdomains_22examples.txt')
    CDSaaFile=file.path(dataFolder, 
        'Homo_sapiens.GRCh37.74_chromosome16_35Mlong.gtf.gz_AAseq.txt.gz')

    inputLoci = read.table(inputFile_loci, sep='\t', 
        stringsAsFactors=F, header=T)
    genomicLoci = proteinLocsToGenomic(inputLoci=inputLoci, 
        CDSaaFile=CDSaaFile)

    kk = paste(genomicLoci[1,1:11], collapse='_')
    checkEquals(kk, 
'chr16_3166431_3166529_+_Exon_4_Exon_4_ENST00000382192_123_155_PF01352_ZNF205')
    kk = paste(genomicLoci[22,1:11], collapse='_')
    checkEquals(kk, 
'chr16_3188540_3190823_+_Exon_1_Exon_3_ENST00000416391_41_127_PF02023_ZNF213')
    kk = paste(genomicLoci[13,1:11], collapse='_')
    # one in reverse strand
    checkEquals(kk, 
'chr16_3274039_3274074_-_Exon_4_Exon_4_ENST00000396868_335_346_PF13912_ZNF200')

}

test_proteinLocsToGenomic_proteinID <- function() {

    dataFolder = system.file("extdata", package="geno2proteo")
    inputFile_loci=file.path(dataFolder, 
    'transId_pfamDomainStartEnd_chr16_Zdomains_22examples_proteinID.txt')
    CDSaaFile=file.path(dataFolder, 
    'Homo_sapiens.GRCh37.74_chromosome16_35Mlong.gtf.gz_AAseq.txt.gz')

    inputLoci = read.table(inputFile_loci, sep='\t', 
        stringsAsFactors=F, header=T)
    genomicLoci = proteinLocsToGenomic(inputLoci=inputLoci, 
        CDSaaFile=CDSaaFile)

    kk = paste(genomicLoci[1,1:11], collapse='_')

    checkEquals(kk, 
'chr16_3166431_3166529_+_Exon_4_Exon_4_ENSP00000371627_123_155_PF01352_ZNF205')

    kk = paste(genomicLoci[22,1:11], collapse='_')
    checkEquals(kk, 
'chr16_3188540_3190823_+_Exon_1_Exon_3_ENSP00000403892_41_127_PF02023_ZNF213')

}


test_proteinLocsToProteinSeq <- function() {

    dataFolder = system.file("extdata", package="geno2proteo")
    inputFile_loci=file.path(dataFolder, 
        "transId_pfamDomainStartEnd_chr16_Zdomains_22examples.txt")
    CDSaaFile=file.path(dataFolder, 
        "Homo_sapiens.GRCh37.74_chromosome16_35Mlong.gtf.gz_AAseq.txt.gz")

    inputLoci = read.table(inputFile_loci, sep='\t', 
        stringsAsFactors=F, header=T)
    ProtSeqNow = proteinLocsToProteinSeq(inputLoci=inputLoci, 
        CDSaaFile=CDSaaFile)

    kk = as.character(ProtSeqNow[,ncol(ProtSeqNow)])
    checkEquals(kk[1], 'PVTFEDVALYLSREEWGRLDHTQQNFYRDVLQK')
    checkEquals(substr(kk[22], 1, 43), 
    'AGPVALGDIPFYFSREEWGTLDPAQRDLFWDIKRENSRNTTLG')
    checkEquals(kk[13], 'FKCPECGKTFPK') # one in reverse strand

}



