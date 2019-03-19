args = commandArgs(trailingOnly=TRUE)
library(gsubfn)
library(data.table)
library(tidyr)


reorderdf<-function(df) { #sort query/subject genomes alphabetically'; assign qcov=NA accordingly
  ordervec<-order(df[c(1,2)])
  df[c(1,2)]<-df[ordervec]
  #df[c(7,8)]<-df[ordervec+6] #swap qlen/slen if alphabetic sort has occurred
  if (ordervec[1]==2) { #if alphabetic re-ordering has taken place; assign qcov to NA
    df[13]<-NA
  }
  return(df)
}


report<-fread(gsubfn('%1',list('%1'=args[1]),'%1/blast/allalignments.tsv'),sep='\t',header=TRUE) #select=c(1,2)

#order margin 1 (rowwise) qseqid/sseqid; replace qcovhsp with NA accordingly (if alphabetic sort has occurred so query coverage is unknown)
report<-as.data.frame(t(apply(report,1,reorderdf)))
colnames(report)<-c('qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','qcovhsp','qlen','slen')

#ordering below doesn't work (not sure why); solved by splitting data prior to ordering
split1<-as.data.frame(separate(report,qseqid,into=c("qseqid", "qseqprot","qseqprotstrand","qseqprotposition"),sep="\\|"))[,c(1:4)]
split2<-as.data.frame(separate(report,sseqid,into=c("sseqid", "sseqprot","sseqprotstrand","sseqprotposition"),sep="\\|"))[,c(2:5)]
report<-data.frame(split1,split2,report[,3:ncol(report)])

#order columnwise; qseq hit / sseq hit
report2<-report[with(report, order(report$qseqid,report$sseqid,report$qseqprotposition,report$sseqprotposition)),]

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