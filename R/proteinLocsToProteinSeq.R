proteinLocsToProteinSeq <-
    function(inputLoci, CDSaaFile) {

    #first get the genomic location of those protein regions
    if(!is.data.frame(inputLoci)) {
        inputLoci = data.frame(inputLoci, stringsAsFactors=FALSE)
    }   
    genomicLoci=proteinLocsToGenomic(inputLoci=inputLoci, CDSaaFile=CDSaaFile)
    
    proteinSeq = genomicLocsToProteinSequence(inputLoci=genomicLoci, 
                        CDSaaFile=CDSaaFile)

    idCol =  as.character(proteinSeq[,7])
    idColNow=as.character(proteinSeq[,7+ncol(inputLoci)])
    # if using the protein ID, then have to convert them to transcript ID
    if(grepl('^ENSP',idCol)[1]) { 
        exonsPep_0 = read.table(CDSaaFile, stringsAsFactors=FALSE, 
                sep='\t', na.strings=NULL)
        kk = exonsPep_0[,5]
        names(kk) = exonsPep_0[,13]
        idCol=kk[idCol]
    }
    kk = sapply(1:length(idCol), function(n) grepl(idCol[n], idColNow[n])[1])
    kkt = proteinSeq[kk,]

    kk = apply(inputLoci, 1, function(a) paste(gsub(' ', '', a), collapse='_'))
    kkt1 = apply(kkt[,7:(7+ncol(inputLoci)-1)], 1, function(a) 
                    paste(gsub(' ', '', a), collapse='_'))
    kkt2 = rep('', length(kk))
    names(kkt2) = kk
    kkt2[kkt1] = as.character(kkt[,'pepSeq'])

    lociExt_0 = cbind(inputLoci, pepSeq = kkt2, stringsAsFactors=FALSE)
    rownames(lociExt_0) = NULL
    
    lociExt_0    
}

