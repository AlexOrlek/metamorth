args = commandArgs(trailingOnly=TRUE)
library(gsubfn)
library(data.table)
library(tidyr)


reorderdf<-function(df) { #sort query/subject genomes alphabetically'; assign qcov=NA accordingly
  ordervec<-order(df[c(1,2)])
  df[c(1,2)]<-df[ordervec]
  df[c(7,8)]<-df[ordervec+6]
  if (ordervec[1]==2) { #if alphabetic re-ordering has taken place; assign qcov to NA
    df[6]<-NA
  }
  return(df)
}


report<-fread(gsubfn('%1',list('%1'=args[1]),'%1/blast/allalignments.tsv'),sep='\t',header=TRUE) #select=c(1,2)

#order margin 1 (rowwise) qseqid/sseqid; replace qcovhsp with NA accordingly (if alphabetic sort has occurred so query coverage is unknown); swap qlen/slen if alphabetic sort has occurred
report<-as.data.frame(t(apply(report,1,reorderdf)))
colnames(report)<-c('sample1','sample2','pident','evalue','bitscore','qcovhsp','qlen','slen')

#ordering below doesn't work (not sure why); solved by splitting data prior to ordering
split1<-as.data.frame(separate(report,sample1,into=c("sample1", "sample1protein","sample1strand","sample1position"),sep="\\|"))[,c(1:4)]
split2<-as.data.frame(separate(report,sample2,into=c("sample2", "sample2protein","sample2strand","sample2position"),sep="\\|"))[,c(2:5)]
report<-data.frame(split1,split2,report[,3:ncol(report)])

#order columnwise; sample1 hit / sample2 hit
report2<-report[with(report, order(report$sample1,report$sample2,report$sample1position,report$sample2position)),]

#remove non-duplicate rows (non-reciprocal hits); only order based on columns 1:8 not evalue/bitscore etc (may depend on blast directionality)
report3<-report2[duplicated(report2[,1:8]) | duplicated(report2[,1:8], fromLast=TRUE),]

#remove duplicate reciprocal - remove NAs since for each duplicate there will be one with qcov and one with qcov NA
report4<-report3[!is.na(report3$qcovhsp),]

write.table(report4, file=gsubfn('%1',list('%1'=args[1]),'%1/blast/allalignments_RBH.tsv'), sep='\t', quote=F, col.names=TRUE, row.names=FALSE)



#OLD CODE
#colnames(report)<-c('sample1','sample2','pid','qcov','qlen','slen')
#report<-report[,c(1,2)]

#report[,c(1,2)]<-as.data.frame(t(apply(report[,c(1,2)],1,sort)),stringsAsFactors = FALSE)
#colnames(report)<-c('sample1','sample2','pident','evalue','bitscore','qcovhsp','qlen','slen')

#report2<-report[with(report, order(report$sample1,report$sample2)),]

#report3<-report2[duplicated(report2) | duplicated(report2, fromLast=TRUE),]

#report4<-report3[!duplicated(report3),]