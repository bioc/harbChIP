
setMethod("show", "upstreamSeqs", function (object) 
{
#    cat("upstreamSeqs instance, organism ", organism(object),  Oct 07 conflict
#       with annotate::organism, eventually need namespace
    cat("upstreamSeqs instance, organism ", object@organism, 
        "\n")
    cat("There are ", tmpl <- length(seqs(object)), " entries\n")
    cat("first keys: \n")
    print(keys(object)[1:min(5, tmpl)])
#    cat(sum(unlist(as.list(isRevcomp(object)))), " sequences are reverse complement\n")
})

setGeneric("keys", function(x) standardGeneric("keys"))
setMethod("keys", "upstreamSeqs", function(x) ls(seqs(x)))
setGeneric("seqs", function(x) standardGeneric("seqs"))
setMethod("seqs", "upstreamSeqs", function(x) x@seqs)
setGeneric("organism", function(x) standardGeneric("organism"))
setMethod("organism", "upstreamSeqs", function(x) x@organism)

getUpstream = function (orfs, upstrob) 
{
    mget(orfs, seqs(upstrob))
}

setGeneric("Nmers", function(n, orf, usobj) standardGeneric("Nmers"))
setMethod("Nmers", c("numeric", "character", "upstreamSeqs"), function(n, orf,
usobj) { 
 if (length(orf)>1) stop("need single orf name")
 if (length(n)>1) stop("n must be numeric length 1")
 bs = getUpstream(orf, usobj)[[1]]
 ## Views(bs, start=1:(Biostrings::nchar(bs)-n+1),
 ## end=n:(Biostrings::nchar(bs)))  }) 
 Views(bs, start=1:(nchar(bs)-n+1), end=n:(nchar(bs)))
})
allhex = function(orf, usobj) Nmers(6, orf, usobj)
