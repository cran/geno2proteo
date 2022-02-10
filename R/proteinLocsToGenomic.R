proteinLocsToGenomic <- function(inputLoci, CDSaaFile) {

    exonsPep_0 = read.table(CDSaaFile, stringsAsFactors=FALSE, sep='\t', 
                            na.strings=NULL)

    kk = exonsPep_0[,1]
    if(sum(grepl('^chr', kk)) ==0 ) {
        kk[kk=='MT'] = 'M'
        kk = paste('chr', kk, sep='')
        exonsPep_0[,1] = kk
    }

    loci_0 = inputLoci
    if(!is.data.frame(loci_0)) {
        loci_0 = data.frame(loci_0, stringsAsFactors=FALSE)
    }

    loci_0[,1] = as.character(loci_0[,1]);
    loci_0[,2] = as.numeric(as.character(loci_0[,2])); 
    loci_0[,3] = as.numeric(as.character(loci_0[,3]))

    rownames(loci_0) = 1:nrow(loci_0)
    
    if(grepl('^ENSP', loci_0[1,1])[1] ) {
        kk = exonsPep_0[,5]
        exonsPep_0[,5] = exonsPep_0[,13]
        exonsPep_0[,13] = kk
    }
    # get the loci in the list that are in the data with protein sequence.
    loci = loci_0[loci_0[,1]%in% exonsPep_0[,5],] 
    if(nrow(loci)==0) {
        return("there is no valid loci!!")
    }
    # get the protein sequence of the transcripts which are in the loci.
    exonsPep = exonsPep_0[exonsPep_0[,5] %in% loci[,1],] 

    kk = paste(exonsPep[,5], exonsPep[,6], sep='_')
    rownames(exonsPep) = kk

    ## now need get the genomic coordinates of the protein loci.
    kk = nchar(exonsPep[,12]) # get the length of the peptides for each exon
    kkt = tapply(kk, exonsPep[,5], function(a) a)
    kkt1 = lapply(kkt, function(a) cumsum(a) )
    #get the start and end of the peptides for each exon in a transcript
    pepExonSE = lapply(kkt1, function(a) {
        kk= length(a); if(kk==1) cbind(1, a) else cbind(c(1, a[1:(kk-1)]+1), a)
    } ) 
    # assume that the order of exons have been done already in the data file, 
    # from the 1st to the last in the fashion of from 5' to 3' end. 
    kkt13 = names(pepExonSE)
    for(ind in 1:length(kkt13) )  { 
        kkt14=kkt13[ind]; kk=nrow(pepExonSE[[ind]]); 
        kkt15=paste(kkt14, '_', 1:kk, sep=''); 
        rownames(pepExonSE[[ind]])=kkt15 
    }
    # get the exon and the position in the exon for the start and 
    # end of peptide.
    loci2Exon = list() 
    for(ind in 1:nrow(loci) ) {
        a = loci[ind,]
        kk = pepExonSE[[a[1,1]]]; 
        kkt = which(a[1,2]>= kk[,1] & a[1,2] <= kk[,2])[1]; 
        kkt01 = a[1,2] - kk[kkt,1] +1
        kkt1 = which(a[1,3]>= kk[,1] & a[1,3] <= kk[,2])[1]; 
        kkt02 = a[1,3]- kk[kkt1,1] +1
        kkt03 = rbind(c(kkt, kkt01), c(kkt1, kkt02))
        rownames(kkt03) = rownames(kk)[kkt03[,1]]
        loci2Exon[[ind]]= kkt03
    }
    names(loci2Exon) = loci[,1]
    # now get the genomic position of the start and end of the peptide
    pepGenomicPos = matrix(nrow=length(loci2Exon), ncol=6)
    for(ind00 in 1:2) {
        if(ind00 ==1) { #for the start
            kk = sapply(loci2Exon, function(a) rownames(a)[1]); # for the start
            # start with the shift in aa
            kk9991 = sapply(loci2Exon, function(a) a[1,2]); names(kk9991) = kk
        } else { # for the end
            kk = sapply(loci2Exon, function(a) rownames(a)[2]); # for the end
            # end with the shift in aa
            kk9991 = sapply(loci2Exon, function(a) a[2,2]); names(kk9991) = kk 
        }
        kk9992 = exonsPep[names(kk9991),]
        # clean up the previous Ns for the current CDS.
        kk = kk9992[,9]; kk[kk==' '] = ''; kk9992[,9] = kk 
        for(ind in seq_along(kk9991)) { # for each locus
        kkt = kk9991[ind]; # the index of the aa in the peptide 
        if(is.na(kkt)) next;

        kkt = 3*kkt; # the length of condon 
        kkt01 = kk9992[ind,]; # the current exon
        rownames(kkt01)  = names(kk9991)[ind] # get the original row names.

        prevNs = kkt01[1,9] # the Ns as the remaining from the preceding exon
        # if it's the first AA for the current exon and there are some 
        # remaining Ns for the current exon, then the start position should be 
        # in the previous exon if there is one.
        if(kkt==3 & prevNs!='') { 
            # get the current exon's name
            kkt02 = strsplit(rownames(kkt01), split='_')[[1]] 
            if(kkt01[1,4] == '+'){
                kkt02[2] = as.numeric(kkt02[2]) - 1 
            } else {
                # get the preceding exon's name
                kkt02[2] = as.numeric(kkt02[2]) + 1 
            }
            kkt02 = paste(kkt02, collapse='_')
            # if there is the preceding exon.
            if(kkt02 %in% rownames(exonsPep)) { 
                # get the preceding exon,as the current exon.
                kkt01 = exonsPep[kkt02,] 
                if(kkt01[1,4]=='+') {
                    kkt02 = kkt01[1,3] -nchar(prevNs)+1
                } else {
                    kkt02 = kkt01[1,2] + nchar(prevNs) - 1
                }
            } else { # if there is no preceding exon, the first AA would be X.
                # Here is one way to deal with, i.e. using the start of exon 
                # and ignoring the added Ns
                if(kkt01[1,4]=='+') {
                    kkt02 = kkt01[1,2]
                } else {
                    kkt02 = kkt01[1,3]
                }
            }
        } else {
            # considering the Ns added to the start of the DNA seq
            kkt = kkt - nchar(prevNs)  
            if(ind00==1) { # for the start position
                kkt = kkt - 3; #the distance to the starting genomic position
            } else { # for the end position
                kkt = kkt - 1
            }
            if(kkt<0) kkt=0;

            if(kkt01[1,4]=='+') {
                # the genomic position starting the peptide
                kkt02 = kkt01[1,2] + kkt
            } else { # minus strand
                kkt02 = kkt01[1,3] -kkt
            } 
        } # for the else

        pepGenomicPos[ind,ind00+1] = kkt02 # get the genomic position now
        kk = rownames(kkt01)
        pepGenomicPos[ind,ind00+4] = sub('^.*_', 'Exon_', kk)

        if(ind00==1) { # for the start position
            pepGenomicPos[ind,1] = kkt01[1,1] # the chromosome
            pepGenomicPos[ind,4] = kkt01[1,4] # the strand
        }
    } # end of the loop for ind
    } # end of the loop for ind00

    colnames(pepGenomicPos) = c('chr', 'start', 'end', 'strand', 
                                'start_exon', 'end_exon')   

    # swap the coordinates in the negative strand
    kk = is.na(pepGenomicPos); kk = rowSums(kk)>0
    # get the number
    kkt = matrix(pepGenomicPos[!kk, ], ncol=ncol(pepGenomicPos) )
    # not get the number
    kkt00 = matrix(pepGenomicPos[kk, ], ncol=ncol(pepGenomicPos) ) 

    kkt01 = kkt[,4] == '-'
    kkt02 = kkt[,2]
    kkt[kkt01,2] = kkt[kkt01,3] 
    kkt[kkt01,3] = kkt02[kkt01] 

    kkt = rbind(kkt, kkt00)
    kkt1 = rownames(loci)
    kkt1= c(kkt1[!kk], kkt1[kk])
    rownames(kkt) = kkt1 # get the rowname from loci.
    kk = matrix(nrow=nrow(loci_0), ncol=ncol(kkt))
    rownames(kk) = rownames(loci_0)
    for(ind in 1:ncol(kkt00)) kk[rownames(kkt),ind] = kkt[,ind]

    pepGenomicPos = cbind(kk, loci_0, stringsAsFactors=FALSE)

    colnames(pepGenomicPos)[1:6] = c('chr', 'start', 'end', 
                                'strand', 'start_exon', 'end_exon')
    
    pepGenomicPos
}
