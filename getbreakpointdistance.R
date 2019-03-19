args = commandArgs(trailingOnly=TRUE)
library(GenomicRanges)
library(gsubfn)
library(data.table)


###ortholog breakpoint function - takes reciprocal best hits file that has been split by pairwise comparison; outputs genome names and bpdiststats

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
  splitnames<-c(splitdf$qseqid[1],splitdf$sseqid[1])
  
  #extracting start/end from position
  s1pos<-as.data.frame(do.call(rbind,sapply(as.vector(splitdf$qseqprotposition), function(x) strsplit(x,'-'))))
  colnames(s1pos)<-c('s1start','s1end')
  s2pos<-as.data.frame(do.call(rbind,sapply(as.vector(splitdf$sseqprotposition), function(x) strsplit(x,'-'))))
  colnames(s2pos)<-c('s2start','s2end')
  
  #load into iranges; add strand|row as name; and order by position
  s1ir<-IRanges(start=(as.numeric(as.vector(s1pos$s1start))), end = (as.numeric(as.vector(s1pos$s1end))), names=paste(splitdf$qseqprotstrand,rownames(splitdf),sep='|'))
  s1out<-names(sort(s1ir))
  s2ir<-IRanges(start=(as.numeric(as.vector(s2pos$s2start))), end = (as.numeric(as.vector(s2pos$s2end))), names=paste(splitdf$sseqprotstrand,rownames(splitdf),sep='|'))
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
  
  #output qseqid sseqid bpdist
  return(c(splitnames,bpdist,breakpoints,genes))
}


reorderdf<-function(df) { #sort query/subject genomes alphabetically'; append NA accordingly
  ordervec<-order(df[c(1,2)])
  df[c(1,2)]<-df[ordervec]
  return(df)
}


###

#read in (reciprocal) best hits file (either calculated by metamorth or user-provided)
if (as.character(args[3])=='metamorth') {
   report<-fread(gsubfn('%1',list('%1'=args[1]),'%1/blast/allalignments_RBH.tsv'),select=1:8,sep='\t',header=TRUE)
} else { #read in best hits file; run getreciprocalhits code in case hits haven't been filtered to select reciprocal-only hits
   report<-fread(gsubfn('%1',list('%1'=args[4]),'%1'),select=c(1,2),sep='\t',header=TRUE) #header=TRUE/no header provided is less problematic than header=FALSE/header provided
   #order margin 1
   report<-as.data.frame(t(apply(report,1,reorderdf)),stringsAsFactors = FALSE)
   colnames(report)<-c('qseqid','sseqid')

   #ordering below doesn't work (not sure why); solved by splitting data prior to ordering                                                                                                                     
   split1<-as.data.frame(separate(report,qseqid,into=c("qseqid", "qseqprot","qseqprotstrand","qseqprotposition"),sep="\\|"))[,c(1:4)]
   split2<-as.data.frame(separate(report,sseqid,into=c("sseqid", "sseqprot","sseqprotstrand","sseqprotposition"),sep="\\|"))[,c(2:5)]
   report<-data.frame(split1,split2)

   #order columnwise; qseqid hit / sseqid hit                                                                                                                                                                
   report2<-report[with(report, order(report$qseqid,report$sseqid,report$qseqprotposition,report$sseqprotposition)),]

   #remove non-duplicate rows (non-reciprocal hits)                                                                                                                                                            
   report3<-report2[duplicated(report2) | duplicated(report2, fromLast=TRUE),]

   #remove duplicate reciprocal                                                                                                                                      
   report4<-report3[!duplicated(report3),]
   report<-report4
   write.table(report, file=gsubfn('%1',list('%1'=args[1]),'%1/blast/allalignments_RBH.tsv'), sep='\t', quote=F, col.names=TRUE, row.names=FALSE)

}

#get included genomes
genomes<-as.data.frame(sort(unique(c(report$qseqid,report$sseqid))))
write.table(genomes, file=gsubfn('%1',list('%1'=args[1]),'%1/included.txt'), sep='\t', quote=F, col.names=FALSE, row.names=FALSE)

#split data by qseqid/sseqid columns (1 split per pairwise comparison)
splitdf<-split(report, list(report$qseqid,report$sseqid))
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
colnames(finaldf)<-c('qseqid','sseqid','Breakpoint_distance','Breakpoints','Genes')

finaldf<-finaldf[with(finaldf, order(finaldf$qseqid,finaldf$sseqid)),]

write.table(finaldf, file=gsubfn('%1',list('%1'=args[1]),'%1/output/breakpointdistance.tsv'), sep='\t', quote=F, col.names=TRUE, row.names=FALSE)





#OLD CODE

#this is unexesssary - don't need NA column
# reorderdf<-function(df) { #sort query/subject genomes alphabetically'; append NA accordingly
#   ordervec<-order(df[c(1,2)])
#   df[c(1,2)]<-df[ordervec]
#   if (ordervec[1]==2) { #if alphabetic re-ordering has taken place; append NA (last column created for NAs)                                                                                                       
#     df[length(df)]<-NA
#   }  
#   return(df)
# }
