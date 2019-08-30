args = commandArgs(trailingOnly=TRUE)
TSS.divide = function(x){
	         log((x+1)/(sum(x)+length(x)))
}
TSS.divide2 = function(x){
	                 (x)/(sum(x))
}
MinMAX = function(x){
		(x-min(x))/(max(x)-min(x))

}
## Read OTU 
df_otu <-read.csv(file=args[1],header = T,sep ='\t',check.names=FALSE)
index =length(df_otu)#-2
SAMPLES <- df_otu[,1]
SAMPLES <- as.list(SAMPLES)
##SAMPLES <- df_otu$SAMPLES
#print(SAMPLES)
#index2=index +1 
#Prog <- df_otu[,c(index2: length(df_otu))]
df_otu<-df_otu[,c(2:index)]
TSS =t(apply(df_otu, 1, TSS.divide))
TSS =apply(TSS, 2, MinMAX)
TSS2 = t(apply(df_otu, 1, TSS.divide2))
#AD<-cbind(SAMPLES,TSS)
AD <- cbind(t(as.data.frame(SAMPLES)), TSS)
colnames(AD)[1] <- c("SAMPLES")
write.table(AD,file = args[2],sep= '\t',row.names = FALSE,quote=FALSE)#,Prog),file = args[2],sep= '\t',row.names = FALSE,quote=FALSE)
AD2 <-cbind(t(as.data.frame(SAMPLES)), TSS2)
colnames(AD2)[1] <- c("SAMPLES")
write.table(AD2,file = args[3],sep= '\t',row.names = FALSE,quote=FALSE)
