
buildUpstreamSeqs2 = function (fastaRead, organism="sce", provenance="harmen") 
{
# apply to a Biostrings::readFASTA object
    norf = length(fastaRead)
    cat("starting transformation to DNAString...")
    bstringlist = lapply(fastaRead, function(z) new("DNAString", 
        z$seq))
    cat("done.")
    fastaDesc = sapply(fastaRead, function(x) x$desc)
    orf = gsub(" .*", "", fastaDesc)
    orfnames = substr(orf, 2, nchar(orf))
    chr = gsub(".*Chr ", "", fastaDesc)
    chromvec = gsub(" .*$", "", chr)
    seqs = new.env(hash = TRUE)
    chrom = new.env(hash = TRUE)
    revComp = new.env(hash = TRUE)
    type = new.env(hash = TRUE)
    for (i in 1:norf) {
        if (i%%100 == 0) 
            cat(i)
        assign(orfnames[i], bstringlist[[i]], seqs)
        assign(orfnames[i], chromvec[i], chrom)
    }
    new("upstreamSeqs", seqs = seqs, chrom = chrom, organism = organism, 
        provenance = provenance)
}

