genomicLocsToProteinSequence <-
    function(inputLoci, CDSaaFile) {

    loci_0 = inputLoci
    if(!is.data.frame(loci_0)) {
        loci_0 = data.frame(loci_0, stringsAsFactors=FALSE)
    }
    loci_0[,1] = as.character(loci_0[,1]);
    loci_0[,2] = as.numeric(as.character(loci_0[,2])); 
    loci_0[,3] = as.numeric(as.character(loci_0[,3]))

    # make sure that the original loci were named by their row index.
    rownames(loci_0) = 1:nrow(loci_0) 
    # removing those loci which don't have a proper location.
    kk = is.na(loci_0[,1:4]); kk = rowSums(kk)>0;  loci = loci_0[!kk,] 

    if(nrow(loci)==0) {
        return("there is no valid loci!!")
    }

    exonsPep = read.table(CDSaaFile, stringsAsFactors=FALSE, sep='\t')

    kk = exonsPep[,1]
    if(sum(grepl('^chr', loci[1,1])) > 0 ) {
        kk[kk=='MT'] = 'M'
        kk = paste('chr', kk, sep='')
        exonsPep[,1] = kk
    }
    
    exonsPep_gr = GenomicRanges::GRanges(exonsPep[,1], 
            IRanges::IRanges(exonsPep[,2], exonsPep[,3]), strand=exonsPep[,4])  
    loci_GR =  GenomicRanges::GRanges(loci[,1], 
    IRanges::IRanges(as.numeric(loci[,2]), as.numeric(loci[,3])), 
        strand=loci[,4])
        
    # the overlap of the loci with CDS, kk000433[[n]]
    kk = GenomicRanges::findOverlaps(loci_GR, exonsPep_gr, minoverlap=1) 
    kkt1 = S4Vectors::queryHits(kk)
    kkt2 = S4Vectors::subjectHits(kk)
    kk01 = cbind(kkt1,kkt2)
    # for each locus, the index of overlapping CDS,  kk0004330[[n]]
    kk02 = tapply(kk01[,2], kk01[,1], sort) 
    # the index of the loci in the overlapping  
    kk021 = as.numeric(names(kk02)) 

    # now checking the exons overlapping with one locus and get the DNA and 
    # AA sequences from those exons
    DNASeqL = list()
    ProteinSeqL = list()
    DNASeqBeforeL = list()
    DNASeqAfterL = list()
    for(locId in 1:length(kk021)) {
        lociNow = loci[kk021[locId], ]
        exonsNow = exonsPep[kk02[[locId]],]
        
        DNAseqNow = vector(length= nrow(exonsNow))
        AAseqNow = vector(length= nrow(exonsNow))
        DNAseqBeforeNow = vector(length= nrow(exonsNow))
        DNAseqEndNow = vector(length= nrow(exonsNow))
        # initialising the start and end with the exon's one
        startNow = exonsNow[,2]; endNow=exonsNow[,3]; strandNow = exonsNow[,4] 
        # get the DNA and AA seq from each exon for the current locus
        for(ind in 1:nrow(exonsNow) ) { 
            # first get the start and end positions of the locus relative 
            # to the 5' side start of the current exon. 
            if(strandNow[ind]=='+') { #for the forward strand
                stCur = startNow[ind] # get the start of exon
                if (lociNow[1,2] > startNow[ind]) stCur = lociNow[1,2]
                endCur = endNow[ind] # get the end of exon
                if(lociNow[1,3] < endNow[ind] ) endCur = lociNow[1,3]

                stCur = stCur - startNow[ind]+1
                endCur = endCur - startNow[ind]+1 
            } else { # for the reverse strand, everything is in reverse now.
                stCur = endNow[ind] # get the start from the exon
                if (lociNow[1,3] < endNow[ind]) stCur = lociNow[1,3]
                endCur = startNow[ind] # get the end
                if(lociNow[1,2] > startNow[ind] ) endCur = lociNow[1,2]
                # get the relative distance to the start of the exon in 
                # the reverse strand.
                stCur = -stCur + endNow[ind]+1 
                endCur = -endCur + endNow[ind]+1 
            }
            # DNA sequences is straightforward
            DNAseqNow[ind] = substr(exonsNow[ind,8], stCur, endCur) 
            
            #now for the protein sequence, which is a bit tricky.
            # the DNA from the beginning of exon to the start of locus
            dna2St = substr(exonsNow[ind,8], 1, stCur)
            # the DNA from the beginning of exon to the end of locus
            dna2End = substr(exonsNow[ind,8], 1, endCur)
            # get the DNA which are the remaining of the preceding exon
            prevDNA = gsub(' ', '', exonsNow[ind,9]) 
            # adding the preceding DNA
            dna2St = paste(prevDNA, dna2St, sep='')  
            dna2End = paste(prevDNA, dna2End, sep='') # similarly
            # how many of the bases in front of the current base are needed 
            # for the start base to form an AA.
            kkt01 = (nchar(dna2St)-1)%% 3
            # how many of the bases after the current one are needed for 
            # the end base to form an AA 
            kkt02 = nchar(dna2End)%% 3; 
            kkt02 = 3 - kkt02; 
            if(kkt02==3) kkt02=0 
            
            if(kkt01==0) { 
                DNAseqBeforeNow[ind]=''
            } else {
                DNAseqBeforeNow[ind] = substr(dna2St, 
                        stCur +nchar(prevDNA)-kkt01, stCur+nchar(prevDNA)-1)
            }
            if(kkt02==0) { 
                DNAseqEndNow[ind]=''
            } else {
                DNAseqEndNow[ind] = substr(exonsNow[ind,8], 
                    endCur+1, endCur+kkt02)
            }
            # get the start of the AA in the current AA sequence
            kkt03 = ceiling(nchar(dna2St)/ 3) 
            # get the end of the AA in the current AA sequence. Note not 
            # using the end base if it have to borrow bases after it.
            kkt04 = floor(nchar(dna2End)/ 3) 
            
            AAseqNow[ind] = substr(exonsNow[ind,12], kkt03, kkt04)
        }
        # now need to put together the DNA and protein sequences belonging 
        # to one transcript.
        exons2transNow= tapply(1:nrow(exonsNow), exonsNow[,5], function(a) a)
        
        DNASeqL[[locId]] = sapply(exons2transNow, function(a) 
                                        paste(DNAseqNow[a], collapse='') )
        ProteinSeqL[[locId]] = sapply(exons2transNow, function(a) 
                                        paste(AAseqNow[a], collapse='') )
        DNASeqBeforeL[[locId]] = sapply(exons2transNow, function(a)  
                                                DNAseqBeforeNow[a[1]])
        DNASeqAfterL[[locId]] = sapply(exons2transNow, function(a) 
                                                DNAseqEndNow[a[length(a)]])
    } 
    seqLst = list()
    for(ind in 1:length(kk02) ) {
        seqLst[[ind]] = cbind(dnaSeq=DNASeqL[[ind]], 
            dnaBefore=DNASeqBeforeL[[ind]], dnaAfter= DNASeqAfterL[[ind]], 
            pepSeq= ProteinSeqL[[ind]])
    }  
    seqLst_unique = lapply(seqLst, function(a) {
        kkt00 = tapply(rownames(a), a[,1], 
                        function(a1) paste(a1, collapse=';'))
        kkt = names(kkt00);
        names(kkt00)=''; kkt01 = match(kkt, a[,1]); kkt02=a[kkt01,]; 
        if(is.matrix(kkt02)) {
            cbind(transId=kkt00,kkt02) 
        } else {
            kkt02 = matrix( c(kkt00,kkt02), nrow=1); 
            colnames(kkt02) = c('transId', colnames(a)); kkt02
        } 
    } )
    #get the loci that have a matching protein sequence.
    kk03 = loci[as.numeric(names(kk02)),]    
    #then get the augmented loci for the protein sequences
    kk = sapply(seqLst_unique, nrow)
    kkt = rep(1:length(seqLst_unique), times=kk)
    kk031 = kk03[kkt,]
    kk = lapply(seqLst_unique, function(a) t(a) )
    kkt = unlist(kk)
    kk032 = t(matrix(kkt, nrow=5))
    colnames(kk032) = rownames(kk[[1]])
    kk032 = cbind(kk031, kk032)
    # get the other rows in the original list
    kk = rownames(kk032)
    kkt = as.numeric(kk); kkt = floor(kkt)
    kkt01 = setdiff(rownames(loci_0), as.character(kkt) )

    kkt02 = loci_0[kkt01,]

    kk= cbind(kkt02, matrix(nrow=nrow(kkt02), ncol=5) ) 
    colnames(kk) = colnames(kk032)
    kkt = rbind(kk032, kk)

    kk = rownames(kkt)
    kk = as.numeric(kk)
    kkt01 = kkt[order(kk),]
    lociExt = kkt01
    lociExt
}
