library(metagenomeSeq)
args = commandArgs(trailingOnly=TRUE)
data <- read.table(file = args[1],
		   header=TRUE,sep= "\t",row.names = "SAMPLES",check.names = FALSE)

data.metagenomeSeq = newMRexperiment(t(data),featureData=NULL, libSize=NULL, normFactors=NULL)
p = cumNormStat(data.metagenomeSeq)
data.cumnorm = cumNorm(data.metagenomeSeq, p=p)
data.CSS = t(MRcounts(data.cumnorm, norm=TRUE, log=TRUE))
SAMPLES <- row.names(data.CSS)
AD <-as.data.frame(SAMPLES)
SUM<-cbind(AD,data.CSS)
write.table(SUM,file= args[2],sep='\t',row.names = FALSE, quote= FALSE)
