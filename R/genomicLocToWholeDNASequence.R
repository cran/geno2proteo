genomicLocToWholeDNASequence <-
    function(inputLoci, DNAfastaFile, tempFolder='./', perlExec='perl') {

    scriptFolder= system.file("exec", package="geno2proteo")
    # the perl script used to get the DNA for the CDS from the fasta file.
    perlScript_getDNAseq=file.path(scriptFolder, 
                    "getDNAseqsFromGenome_loci.perl")
    # get the genomic loci
    loci_0 = inputLoci
    if(!is.data.frame(loci_0)) loci_0 = data.frame(loci_0,
        stringsAsFactors=FALSE)
    loci_0[,1] = as.character(loci_0[,1]);
    loci_0[,2] = as.numeric(as.character(loci_0[,2])); 
    loci_0[,3] = as.numeric(as.character(loci_0[,3]))
    
    # make sure that the original loci were named by their row index.
    rownames(loci_0) = 1:nrow(loci_0)
    # removing those loci which don't have a proper location.
    kk = is.na(loci_0[,1:4]); kk = rowSums(kk)>0;  
    loci = loci_0[!kk,] 

    if(nrow(loci)==0) {
        return("there is no valid loci!!")
    }
    filenamePrefix = 'inputFile_loci'
    filenameGenomicLoc= paste(tempFolder, '/', 
            filenamePrefix, '_genomicLoc.txt', sep='');
    
    write.table(loci, file=filenameGenomicLoc, sep='\t', quote=FALSE, 
                    row.names=FALSE, col.names=FALSE)
    
    # then run the perl script to get the DNA sequence of the regions
    filenameDnaSeq = paste(tempFolder, '/', filenamePrefix, 
                '_genomicLoc_DNAseq.txt', sep='');
    
    DNAfastaFile00 = DNAfastaFile
    # uncompress the fasta file, not remove the compressed
    if(grepl('\\.gz$', DNAfastaFile)[1]) {
        R.utils::gunzip(DNAfastaFile, remove=FALSE)
        DNAfastaFile = sub('\\.gz$','', DNAfastaFile)
    }
    
    command = paste('"', perlExec, '" "', perlScript_getDNAseq, '" "', 
    filenameGenomicLoc, '" "',DNAfastaFile, '" "',filenameDnaSeq, '"', sep='')

    system(command)

    if(grepl('\\.gz$', DNAfastaFile00)[1]) { # if uncompressed the fasta file 
        if (file.exists(DNAfastaFile)) { invisible(file.remove(DNAfastaFile))}
        DNAfastaFile= DNAfastaFile00
    }

    lociExt = read.table(filenameDnaSeq, sep='\t', stringsAsFactors=FALSE)

    kk = rep('NA', nrow(loci_0))
    kk[as.numeric(rownames(loci))] = lociExt[,ncol(lociExt)]
    # append the DNA seq to the end of the row.
    lociExt_0 = cbind(loci_0, DNAseq = kk) 

    # removing the tempory files
    fn=filenameGenomicLoc; if (file.exists(fn)) invisible(file.remove(fn))
    fn=filenameDnaSeq; if (file.exists(fn)) invisible(file.remove(fn))

    lociExt_0
}
