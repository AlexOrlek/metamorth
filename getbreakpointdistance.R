args = commandArgs(trailingOnly=TRUE)
library(GenomicRanges)
library(gsubfn)
library(data.table)


###ortholog breakpoint function - takes reciprocal best hits file that has been split by pairwise comparison; outputs sample names and bpdiststats

#pre-requisite functions
makepairs<-function(x) mapply(c, head(x,-1), tail(x,-1), SIMPLIFY = FALSE)
BP<-function(x,y) { #this function works on signed permuations (numeric vectors with +/- indicated)                                         
  out1<-makepairs(x)[!(makepairs(x) %in% makepairs(y) | makepairs(x) %in% makepairs(rev(y*-1)))]
  out2<-makepairs(x)[unlist(lapply(makepairs(x), function(z) length(unique(sign(z)))))>1]
  all<-c(out1,out2)
  return(all[!duplicated(all)]) #deduplicate to aovid double counting breakpoints that occur in out1 and out2                               
}

#orthologBP function
orthologBP<-function(splitdf) {
  splitnames<-c(splitdf$sample1[1],splitdf$sample2[1])
  
  #extracting start/end from position
  s1pos<-as.data.frame(do.call(rbind,sapply(as.vector(splitdf$sample1position), function(x) strsplit(x,'-'))))
  colnames(s1pos)<-c('s1start','s1end')
  s2pos<-as.data.frame(do.call(rbind,sapply(as.vector(splitdf$sample2position), function(x) strsplit(x,'-'))))
  colnames(s2pos)<-c('s2start','s2end')
  
  #load into iranges; add strand|row as name; and order by position
  s1ir<-IRanges(start=(as.numeric(as.vector(s1pos$s1start))), end = (as.numeric(as.vector(s1pos$s1end))), names=paste(splitdf$sample1strand,rownames(splitdf),sep='|'))
  s1out<-names(sort(s1ir))
  s2ir<-IRanges(start=(as.numeric(as.vector(s2pos$s2start))), end = (as.numeric(as.vector(s2pos$s2end))), names=paste(splitdf$sample2strand,rownames(splitdf),sep='|'))
  s2out<-names(sort(s2ir))
  
  #convert to strand|row to signed permutation
  s1perm<-as.data.frame(do.call(rbind,sapply(as.vector(s1out), function(x) strsplit(x,'\\|'))),stringsAsFactors=F)
  colnames(s1perm)<-c('strand','index')
  s1perm<-as.numeric(s1perm$strand)*as.numeric(s1perm$index)
  s2perm<-as.data.frame(do.call(rbind,sapply(as.vector(s2out), function(x) strsplit(x,'\\|'))),stringsAsFactors=F)
  colnames(s2perm)<-c('strand','index')
  s2perm<-as.numeric(s2perm$strand)*as.numeric(s2perm$index)
  
  #calculate breakpoint distance from signed permutations
  breakpoints<-length(BP(s1perm,s2perm))
  genes<-length(s1perm)
  stopifnot(genes==length(s2perm))
  bpdist<-breakpoints/genes
  
  #output sample1 sample2 bpdist
  return(c(splitnames,bpdist,breakpoints,genes))
}

###


#read in reciprocal best hits file (either calculated by metamorth or user-provided)
if (as.character(args[3])=='metamorth') {
   report<-fread(gsubfn('%1',list('%1'=args[1]),'%1/blast/allalignments_RBH.tsv'),sep='\t',header=TRUE)
} else {
   report<-fread(gsubfn('%1',list('%1'=args[4]),'%1'),sep='\t',header=FALSE)
   #order margin 1                                                                                                                                                                                             
   report<-as.data.frame(t(apply(report,1,sort)),stringsAsFactors = FALSE)
   colnames(report)<-c('sample1','sample2')

   #ordering below doesn't work (not sure why); solved by splitting data prior to ordering                                                                                                                     
   split1<-as.data.frame(separate(report,sample1,into=c("sample1", "sample1protein","sample1strand","sample1position"),sep="\\|"))[,c(1:4)]
   split2<-as.data.frame(separate(report,sample2,into=c("sample2", "sample2protein","sample2strand","sample2position"),sep="\\|"))[,c(2:5)]
   report<-data.frame(split1,split2)

   #order by column 1 and column 2                                                                                                                                                                             
   report2<-report[with(report, order(report$sample1,report$sample2)),]

   #remove non-duplicate rows (non-reciprocal hits)                                                                                                                                                            
   report3<-report2[duplicated(report2) | duplicated(report2, fromLast=TRUE),]

   #remove duplicate reciprocal                                                                                                                                                                                
   report4<-report3[!duplicated(report3),]
   report<-report4
   write.table(report, file=gsubfn('%1',list('%1'=args[1]),'%1/blast/allalignments_RBH.tsv'), sep='\t', quote=F, col.names=TRUE, row.names=FALSE)

}

#get included samples
samples<-as.data.frame(sort(unique(c(report$sample1,report$sample2))))
write.table(samples, file=gsubfn('%1',list('%1'=args[1]),'%1/included.txt'), sep='\t', quote=F, col.names=FALSE, row.names=FALSE)

#split data by sample1/sample2 columns (1 split per pairwise comparison)
splitdf<-split(report, list(report$sample1,report$sample2))
includedindices<-lapply(splitdf,nrow)>0 #prevents bug due to empty list element
splitdf<-splitdf[includedindices]           

#apply breakpoint calculation function to splitdf
library('foreach')
library('doParallel')

cl<-makeCluster(as.integer(args[2]))
registerDoParallel(cl)
clusterExport(cl,c("makepairs","BP","IRanges"))

finaldf<-parLapply(cl, splitdf,orthologBP)

stopCluster(cl)

#convert to dataframe, order, and write to file

finaldf<-as.data.frame(do.call(rbind,finaldf))
colnames(finaldf)<-c('Sample1','Sample2','Breakpoint_distance','Breakpoints','Genes')

finaldf<-finaldf[with(finaldf, order(finaldf$Sample1,finaldf$Sample2)),]

write.table(finaldf, file=gsubfn('%1',list('%1'=args[1]),'%1/output/breakpointdistance.tsv'), sep='\t', quote=F, col.names=TRUE, row.names=FALSE)

