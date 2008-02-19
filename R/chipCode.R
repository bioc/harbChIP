
#library(harbChIP)
#data(harbChIP)
#data(sceUpstr)

#1 generate and investigate hexamer tables

#library(Biostrings)

#us1 = getUpstream("YAL001C", sceUpstr)[[1]]
#hus1 = views(us1, start=1:495, end=6:500)
#hus1
#sort(table(as.character(hus1)), decr=TRUE)[1:8]

# how many unique hexamers are found in upstream
# regions of YAL001C and YAL002W?

# how many unique hexamers are present in YAL002W but
# absent in YAL001C?  verify that "GGAATC" meets this
# condition

#2 test for independence of heptamer occupancy and
# TF binding intensity for a given heptamer and TF
# example GGCGCTA, SNT2

#myhep = "GGCGCTA"
#countPattern( myhep, getUpstream("YAL001C", sceUpstr)[[1]] )
chkAllUS = function(patt, upstr) {
 orfs = keys(sceUpstr)
 allu = lapply( orfs, function(x) getUpstream(x, sceUpstr)[[1]] )
 names(allu) = orfs
 occ = sapply(allu, function(x) countPattern(patt, x))
 names(occ) = orfs
 occ
}

#UShep = chkAllUS( myhep, sceUpstr )
#summary(UShep)
#BRAT.snt2 = exprs(harbChIP)[,"SNT2"]
#common = intersect(names(UShep), names(BRAT.snt2))
#UShep = UShep[common]
#BRAT.snt2 = BRAT.snt2[common]
#summary(BRAT.snt2)
#hasHep = common %in% names(UShep[UShep>0])
#hiBind = common %in% names(na.omit(BRAT.snt2[BRAT.snt2>2]))
#table(hasHep, hiBind)
#fisher.test(table(hasHep, hiBind))
   
#3 write the software to test for enrichment of a given pattern
# in intergenic regions to which a given TF binds strongly

chkMotif4TF = function(motif, TF, chset, upstr, bthresh=2, countthresh=0) {
 if (!(TF %in% sampleNames(chset))) stop("TF not found in sampleNames of chset")
 cat("generating motif counts ...\n")
 allus = chkAllUS( motif, upstr )
 cat("done.\n")
 bvec = exprs(chset)[,TF]
 common = intersect(names(allus), names(bvec))
 allus = allus[common]
 bvec = bvec[common]
 hasMotif = common %in% names(allus[allus>countthresh])
 hiBind = common %in% names(na.omit(bvec[bvec > bthresh]))
 tt = fisher.test( table(hasMotif, hiBind) )
 tab = table(hasMotif, hiBind)
 ca = match.call()
 list(call=ca, tab=tab, test=tt)
}
 
