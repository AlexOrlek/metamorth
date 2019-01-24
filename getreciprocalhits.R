args = commandArgs(trailingOnly=TRUE)
library(gsubfn)
library(data.table)
library(tidyr)

report<-fread(gsubfn('%1',list('%1'=args[1]),'%1/blast/allalignments.tsv'),select=c(1,2),sep='\t',header=FALSE)
#colnames(report)<-c('sample1','sample2','pid','qcov','qlen','slen')
#report<-report[,c(1,2)]

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

write.table(report4, file=gsubfn('%1',list('%1'=args[1]),'%1/blast/allalignments_RBH.tsv'), sep='\t', quote=F, col.names=TRUE, row.names=FALSE)