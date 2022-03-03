setwd("/mnt/raid5/ChIPseqArchive/ChoMet/")
#install.packages("genome/BSGenome.Pf3D7.PlasmoDB.37_1.0.tar.gz")
require(Rsamtools);require(plyr);library(reshape2); library(dplyr); library(ggplot2); library(ggsignif); library(ggpubr);library(Rmisc);library(tidyr);
library(InPAS); library(BSGenome.Pf3D7.PlasmoDB.37); library(remotes); library(GenomicRanges);library(Gviz);library(rtracklayer);library(Rsamtools);library(GenomeInfoDb)

genome=Pf3D7

#############Loading proper experimental pairs as Granges for analysis 
#Load sample pairs that can be compared
X1.S16.X1.NF54.MC.K9.B_SPMR_FE_comp.bdg<-import.bedGraph("./singleSamplePeaks/SPMR/X1.S16.X1.NF54.MC.K9.B_SPMR_FE_comp.bdg", genome="Pf3D7")
X1.S19.X1.NF54._.K9.B_SPMR_FE_comp.bdg<-import.bedGraph("./singleSamplePeaks/SPMR/X1.S19.X1.NF54._.K9.B_SPMR_FE_comp.bdg", genome="Pf3D7")
X1.S24.X1.NF54.MC.K9.C_SPMR_FE_comp.bdg <-import.bedGraph("./singleSamplePeaks/SPMR/X1.S24.X1.NF54.MC.K9.C_SPMR_FE_comp.bdg", genome="Pf3D7")
X1.S27.X1.NF54._.K9.C_SPMR_FE_comp.bdg<-import.bedGraph("./singleSamplePeaks/SPMR/X1.S27.X1.NF54._.K9.C_SPMR_FE_comp.bdg", genome="Pf3D7")

#X1.S16.X1.NF54.MC.K9.B_SPMR_FE_comp.bdg<-X1.S16.X1.NF54.MC.K9.B_SPMR_FE_comp.bdg[]

seqlevels(X1.S16.X1.NF54.MC.K9.B_SPMR_FE_comp.bdg)<-seqlevels(genome)
seqlevels(X1.S19.X1.NF54._.K9.B_SPMR_FE_comp.bdg)<-seqlevels(genome)
seqlevels(X1.S24.X1.NF54.MC.K9.C_SPMR_FE_comp.bdg)<-seqlevels(genome)
seqlevels(X1.S27.X1.NF54._.K9.C_SPMR_FE_comp.bdg)<-seqlevels(genome)

seqlengths(X1.S16.X1.NF54.MC.K9.B_SPMR_FE_comp.bdg)<- seqlengths(genome)
seqlengths(X1.S19.X1.NF54._.K9.B_SPMR_FE_comp.bdg)<- seqlengths(genome)
seqlengths(X1.S24.X1.NF54.MC.K9.C_SPMR_FE_comp.bdg)<- seqlengths(genome)
seqlengths(X1.S27.X1.NF54._.K9.C_SPMR_FE_comp.bdg)<- seqlengths(genome)

K9.MC.pairs <- bindAsGRanges(rep1_MC=mcolAsRleList(X1.S16.X1.NF54.MC.K9.B_SPMR_FE_comp.bdg,"score"), rep1_=mcolAsRleList(X1.S19.X1.NF54._.K9.B_SPMR_FE_comp.bdg,"score"),rep2_MC=mcolAsRleList(X1.S24.X1.NF54.MC.K9.C_SPMR_FE_comp.bdg,"score"),rep2_=mcolAsRleList(X1.S27.X1.NF54._.K9.C_SPMR_FE_comp.bdg,"score"))
K9.MC.pairs$rep1K9delta<-round(K9.MC.pairs$rep1_-K9.MC.pairs$rep1_MC,4)
K9.MC.pairs$rep2K9delta<-round(K9.MC.pairs$rep2_-K9.MC.pairs$rep2_MC,4)
K9.MC.pairs<-K9.MC.pairs[seqnames(K9.MC.pairs) %in% paste0("Pf3D7_",sprintf("%02.0f",1:14), "_v3")]
seqlevels(K9.MC.pairs)<- paste0("Pf3D7_",sprintf("%02.0f",1:14), "_v3")
genome(K9.MC.pairs)<-"3D7"

#write.table(data.frame(bindAsGRanges(mcolAsRleList(K9.MC.pairs,"rep1K9delta")))[,c("seqnames","start","end","V1")],file=paste0("./byTreatAbPeaks/deltaK9.MC_rep1.bdg"), quote=F, sep="\t", col.names = F, row.names = F)
#write.table(data.frame(bindAsGRanges(mcolAsRleList(K9.MC.pairs,"rep2K9delta")))[,c("seqnames","start","end","V1")],file=paste0("./byTreatAbPeaks/deltaK9.MC_rep2.bdg"), quote=F, sep="\t", col.names = F, row.names = F)
#rm(X1.S16.X1.NF54.MC.K9.B_SPMR_FE_comp.bdg,X1.S19.X1.NF54._.K9.B_SPMR_FE_comp.bdg,X1.S24.X1.NF54.MC.K9.C_SPMR_FE_comp.bdg,X1.S27.X1.NF54._.K9.C_SPMR_FE_comp.bdg)



###K4 delta MC pairs
X1.S4.X1.NF54.MC.K4.A_SPMR_FE_comp.bdg <-import.bedGraph("./singleSamplePeaks/SPMR/X1.S4.X1.NF54.MC.K4.A_SPMR_FE_comp.bdg", genome="3D7")
seqlevels(X1.S4.X1.NF54.MC.K4.A_SPMR_FE_comp.bdg)<-seqlevels(genome)
seqlengths(X1.S4.X1.NF54.MC.K4.A_SPMR_FE_comp.bdg)<- seqlengths(genome)

X1.S8.X1.NF54._.K4.A_SPMR_FE_comp.bdg<-import.bedGraph("./singleSamplePeaks/SPMR/X1.S8.X1.NF54._.K4.A_SPMR_FE_comp.bdg", genome="3D7")
seqlevels(X1.S8.X1.NF54._.K4.A_SPMR_FE_comp.bdg)<-seqlevels(genome)
seqlengths(X1.S8.X1.NF54._.K4.A_SPMR_FE_comp.bdg)<- seqlengths(genome)

X1.S15.X1.NF54.MC.K4.B_SPMR_FE_comp.bdg <-import.bedGraph("./singleSamplePeaks/SPMR/X1.S15.X1.NF54.MC.K4.B_SPMR_FE_comp.bdg", genome="3D7")
seqlevels(X1.S15.X1.NF54.MC.K4.B_SPMR_FE_comp.bdg)<-seqlevels(genome)
seqlengths(X1.S15.X1.NF54.MC.K4.B_SPMR_FE_comp.bdg)<- seqlengths(genome)

X1.S18.X1.NF54._.K4.B_SPMR_FE_comp.bdg<-import.bedGraph("./singleSamplePeaks/SPMR/X1.S18.X1.NF54._.K4.B_SPMR_FE_comp.bdg", genome="3D7")
seqlevels(X1.S18.X1.NF54._.K4.B_SPMR_FE_comp.bdg)<-seqlevels(genome)
seqlengths(X1.S18.X1.NF54._.K4.B_SPMR_FE_comp.bdg)<- seqlengths(genome)

X1.S23.X1.NF54.MC.K4.C_SPMR_FE_comp.bdg <-import.bedGraph("./singleSamplePeaks/SPMR/X1.S23.X1.NF54.MC.K4.C_SPMR_FE_comp.bdg", genome="3D7")
seqlevels(X1.S23.X1.NF54.MC.K4.C_SPMR_FE_comp.bdg)<-seqlevels(genome)
seqlengths(X1.S23.X1.NF54.MC.K4.C_SPMR_FE_comp.bdg)<- seqlengths(genome)

X1.S26.X1.NF54._.K4.C_SPMR_FE_comp.bdg<-import.bedGraph("./singleSamplePeaks/SPMR/X1.S26.X1.NF54._.K4.C_SPMR_FE_comp.bdg", genome="3D7")
seqlevels(X1.S26.X1.NF54._.K4.C_SPMR_FE_comp.bdg)<-seqlevels(genome)
seqlengths(X1.S26.X1.NF54._.K4.C_SPMR_FE_comp.bdg)<- seqlengths(genome)


K4.MC.pairs <- bindAsGRanges(
  rep1_MC=mcolAsRleList(X1.S4.X1.NF54.MC.K4.A_SPMR_FE_comp.bdg,"score"), 
  rep1_=mcolAsRleList(X1.S8.X1.NF54._.K4.A_SPMR_FE_comp.bdg,"score"),
  rep2_MC=mcolAsRleList(X1.S15.X1.NF54.MC.K4.B_SPMR_FE_comp.bdg,"score"),
  rep2_=mcolAsRleList(X1.S18.X1.NF54._.K4.B_SPMR_FE_comp.bdg,"score"),
  rep3_MC=mcolAsRleList(X1.S23.X1.NF54.MC.K4.C_SPMR_FE_comp.bdg,"score"),
  rep3_=mcolAsRleList(X1.S26.X1.NF54._.K4.C_SPMR_FE_comp.bdg,"score")
)

K4.MC.pairs$rep1K4delta<-round(K4.MC.pairs$rep1_-K4.MC.pairs$rep1_MC,4)
K4.MC.pairs$rep2K4delta<-round(K4.MC.pairs$rep2_-K4.MC.pairs$rep2_MC,4)
K4.MC.pairs$rep3K4delta<-round(K4.MC.pairs$rep3_-K4.MC.pairs$rep3_MC,4)
genome(K4.MC.pairs)<-"3D7"
K4.MC.pairs<-K4.MC.pairs[seqnames(K4.MC.pairs) %in% paste0("Pf3D7_",sprintf("%02.0f",1:14), "_v3")]
seqlevels(K4.MC.pairs)<- paste0("Pf3D7_",sprintf("%02.0f",1:14), "_v3")

#write.table(data.frame(bindAsGRanges(mcolAsRleList(K4.MC.pairs,"rep1K4delta")))[,c("seqnames","start","end","V1")],file=paste0("./byTreatAbPeaks/deltaK4.MC_rep1.bdg"), quote=F, sep="\t", col.names = F, row.names = F)
#write.table(data.frame(bindAsGRanges(mcolAsRleList(K4.MC.pairs,"rep2K4delta")))[,c("seqnames","start","end","V1")],file=paste0("./byTreatAbPeaks/deltaK4.MC_rep2.bdg"), quote=F, sep="\t", col.names = F, row.names = F)
#write.table(data.frame(bindAsGRanges(mcolAsRleList(K4.MC.pairs,"rep3K4delta")))[,c("seqnames","start","end","V1")],file=paste0("./byTreatAbPeaks/deltaK4.MC_rep3.bdg"), quote=F, sep="\t", col.names = F, row.names = F)
#rm(X1.S4.X1.NF54.MC.K4.A_SPMR_FE_comp.bdg,X1.S8.X1.NF54._.K4.A_SPMR_FE_comp.bdg,X1.S15.X1.NF54.MC.K4.B_SPMR_FE_comp.bdg,X1.S18.X1.NF54._.K4.B_SPMR_FE_comp.bdg,X1.S23.X1.NF54.MC.K4.C_SPMR_FE_comp.bdg,X1.S26.X1.NF54._.K4.C_SPMR_FE_comp.bdg)

#####Reading in MCvsMCDZA pairs, save as Grange, calculate sum of deltas in K9s and K4

X3.S45.X2.NF54.MC.K9.B_SPMR_FE_comp.bdg<-import.bedGraph("./singleSamplePeaks/SPMR/X3.S45.X2.NF54.MC.K9.B_SPMR_FE_comp.bdg", genome="3D7")
X2.S58.X2.NF54.MC.K9.C_SPMR_FE_comp.bdg<-import.bedGraph("./singleSamplePeaks/SPMR/X2.S58.X2.NF54.MC.K9.C_SPMR_FE_comp.bdg", genome="3D7")
X3.S2.X3.NF54.MC.K9.A_SPMR_FE_comp.bdg<-import.bedGraph("./singleSamplePeaks/SPMR/X3.S2.X3.NF54.MC.K9.A_SPMR_FE_comp.bdg", genome="3D7")
X3.S46.X2.NF54.MCDZA.K9.B_SPMR_FE_comp.bdg<-import.bedGraph("./singleSamplePeaks/SPMR/X3.S46.X2.NF54.MCDZA.K9.B_SPMR_FE_comp.bdg", genome="3D7")
X2.S59.X2.NF54.MCDZA.K9.C_SPMR_FE_comp.bdg <-import.bedGraph("./singleSamplePeaks/SPMR/X2.S59.X2.NF54.MCDZA.K9.C_SPMR_FE_comp.bdg", genome="3D7")
X3.S6.X3.NF54.MCDZA.K9.A_SPMR_FE_comp.bdg<-import.bedGraph("./singleSamplePeaks/SPMR/X3.S6.X3.NF54.MCDZA.K9.A_SPMR_FE_comp.bdg", genome="3D7")

seqlevels(X3.S45.X2.NF54.MC.K9.B_SPMR_FE_comp.bdg)<-seqlevels(genome)
seqlevels(X2.S58.X2.NF54.MC.K9.C_SPMR_FE_comp.bdg)<-seqlevels(genome)
seqlevels(X3.S2.X3.NF54.MC.K9.A_SPMR_FE_comp.bdg)<-seqlevels(genome)
seqlevels(X3.S46.X2.NF54.MCDZA.K9.B_SPMR_FE_comp.bdg)<-seqlevels(genome)
seqlevels(X2.S59.X2.NF54.MCDZA.K9.C_SPMR_FE_comp.bdg)<-seqlevels(genome)
seqlevels(X3.S6.X3.NF54.MCDZA.K9.A_SPMR_FE_comp.bdg)<-seqlevels(genome)

seqlengths(X3.S45.X2.NF54.MC.K9.B_SPMR_FE_comp.bdg)<- seqlengths(genome)
seqlengths(X2.S58.X2.NF54.MC.K9.C_SPMR_FE_comp.bdg)<- seqlengths(genome)
seqlengths(X3.S2.X3.NF54.MC.K9.A_SPMR_FE_comp.bdg)<- seqlengths(genome)
seqlengths(X3.S46.X2.NF54.MCDZA.K9.B_SPMR_FE_comp.bdg)<- seqlengths(genome)
seqlengths(X2.S59.X2.NF54.MCDZA.K9.C_SPMR_FE_comp.bdg)<- seqlengths(genome)
seqlengths(X3.S6.X3.NF54.MCDZA.K9.A_SPMR_FE_comp.bdg)<- seqlengths(genome)

K9.MCDZA.pairs <- bindAsGRanges(rep1_MC=mcolAsRleList(X3.S45.X2.NF54.MC.K9.B_SPMR_FE_comp.bdg,"score"), rep1_MCDZA=mcolAsRleList(X3.S46.X2.NF54.MCDZA.K9.B_SPMR_FE_comp.bdg,"score"),rep2_MC=mcolAsRleList(X2.S58.X2.NF54.MC.K9.C_SPMR_FE_comp.bdg,"score"),rep2_MCDZA=mcolAsRleList(X2.S59.X2.NF54.MCDZA.K9.C_SPMR_FE_comp.bdg,"score"),rep3_MC=mcolAsRleList(X3.S2.X3.NF54.MC.K9.A_SPMR_FE_comp.bdg,"score"),rep3_MCDZA=mcolAsRleList(X3.S6.X3.NF54.MCDZA.K9.A_SPMR_FE_comp.bdg,"score"))
K9.MCDZA.pairs$rep1K9delta<-round(K9.MCDZA.pairs$rep1_MCDZA-K9.MCDZA.pairs$rep1_MC,4)
K9.MCDZA.pairs$rep2K9delta<-round(K9.MCDZA.pairs$rep2_MCDZA-K9.MCDZA.pairs$rep2_MC,4)
K9.MCDZA.pairs$rep3K9delta<-round(K9.MCDZA.pairs$rep3_MCDZA-K9.MCDZA.pairs$rep3_MC,4)
K9.MCDZA.pairs<-K9.MCDZA.pairs[seqnames(K9.MCDZA.pairs) %in% paste0("Pf3D7_",sprintf("%02.0f",1:14), "_v3")]
seqlevels(K9.MCDZA.pairs)<- paste0("Pf3D7_",sprintf("%02.0f",1:14), "_v3")
genome(K9.MCDZA.pairs)<-"3D7"

#write.table(data.frame(bindAsGRanges(mcolAsRleList(K9.MCDZA.pairs,"rep1K9delta")))[,c("seqnames","start","end","V1")],file=paste0("./byTreatAbPeaks/deltaK9.MCDZA_rep1.bdg"), quote=F, sep="\t", col.names = F, row.names = F)
#write.table(data.frame(bindAsGRanges(mcolAsRleList(K9.MCDZA.pairs,"rep2K9delta")))[,c("seqnames","start","end","V1")],file=paste0("./byTreatAbPeaks/deltaK9.MCDZA_rep2.bdg"), quote=F, sep="\t", col.names = F, row.names = F)
#write.table(data.frame(bindAsGRanges(mcolAsRleList(K9.MCDZA.pairs,"rep3K9delta")))[,c("seqnames","start","end","V1")],file=paste0("./byTreatAbPeaks/deltaK9.MCDZA_rep3.bdg"), quote=F, sep="\t", col.names = F, row.names = F)
#rm(X3.S45.X2.NF54.MC.K9.B_SPMR_FE_comp.bdg,X2.S58.X2.NF54.MC.K9.C_SPMR_FE_comp.bdg,X3.S2.X3.NF54.MC.K9.A_SPMR_FE_comp.bdg,X3.S46.X2.NF54.MCDZA.K9.B_SPMR_FE_comp.bdg,X2.S59.X2.NF54.MCDZA.K9.C_SPMR_FE_comp.bdg,X3.S6.X3.NF54.MCDZA.K9.A_SPMR_FE_comp.bdg)

#K4
X3.S47.X2.NF54.MC.K4.B_SPMR_FE_comp.bdg<-import.bedGraph("./singleSamplePeaks/SPMR/X3.S47.X2.NF54.MC.K4.B_SPMR_FE_comp.bdg", genome="3D7")
seqlevels(X3.S47.X2.NF54.MC.K4.B_SPMR_FE_comp.bdg)<-seqlevels(genome)
seqlengths(X3.S47.X2.NF54.MC.K4.B_SPMR_FE_comp.bdg)<- seqlengths(genome)
X2.S47.X2.NF54.MCDZA.K4.B_SPMR_FE_comp.bdg <-import.bedGraph("./singleSamplePeaks/SPMR/X2.S47.X2.NF54.MCDZA.K4.B_SPMR_FE_comp.bdg", genome="3D7")
seqlevels(X2.S47.X2.NF54.MCDZA.K4.B_SPMR_FE_comp.bdg)<-seqlevels(genome)
seqlengths(X2.S47.X2.NF54.MCDZA.K4.B_SPMR_FE_comp.bdg)<- seqlengths(genome)

X2.S61.X2.NF54.MCDZA.K4.C_SPMR_FE_comp.bdg <-import.bedGraph("./singleSamplePeaks/SPMR/X2.S61.X2.NF54.MCDZA.K4.C_SPMR_FE_comp.bdg", genome="3D7")
seqlevels(X2.S61.X2.NF54.MCDZA.K4.C_SPMR_FE_comp.bdg)<-seqlevels(genome)
seqlengths(X2.S61.X2.NF54.MCDZA.K4.C_SPMR_FE_comp.bdg)<- seqlengths(genome)
X2.S60.X2.NF54.MC.K4.C_SPMR_FE_comp.bdg <-import.bedGraph("./singleSamplePeaks/SPMR/X2.S60.X2.NF54.MC.K4.C_SPMR_FE_comp.bdg", genome="3D7")
seqlevels(X2.S60.X2.NF54.MC.K4.C_SPMR_FE_comp.bdg)<-seqlevels(genome)
seqlengths(X2.S60.X2.NF54.MC.K4.C_SPMR_FE_comp.bdg)<- seqlengths(genome)

X3.S3.X3.NF54.MC.K4.A_SPMR_FE_comp.bdg<-import.bedGraph("./singleSamplePeaks/SPMR/X3.S3.X3.NF54.MC.K4.A_SPMR_FE_comp.bdg", genome="3D7")
seqlevels(X3.S3.X3.NF54.MC.K4.A_SPMR_FE_comp.bdg)<-seqlevels(genome)
seqlengths(X3.S3.X3.NF54.MC.K4.A_SPMR_FE_comp.bdg)<- seqlengths(genome)
X3.S7.X3.NF54.MCDZA.K4.A_SPMR_FE_comp.bdg<-import.bedGraph("./singleSamplePeaks/SPMR/X3.S7.X3.NF54.MCDZA.K4.A_SPMR_FE_comp.bdg", genome="3D7")
seqlevels(X3.S7.X3.NF54.MCDZA.K4.A_SPMR_FE_comp.bdg)<-seqlevels(genome)
seqlengths(X3.S7.X3.NF54.MCDZA.K4.A_SPMR_FE_comp.bdg)<- seqlengths(genome)


K4.MCDZA.pairs <- bindAsGRanges(rep1_MC=mcolAsRleList(X3.S47.X2.NF54.MC.K4.B_SPMR_FE_comp.bdg,"score"), rep1_MCDZA=mcolAsRleList(X2.S47.X2.NF54.MCDZA.K4.B_SPMR_FE_comp.bdg,"score"),rep2_MC=mcolAsRleList(X2.S60.X2.NF54.MC.K4.C_SPMR_FE_comp.bdg,"score"),rep2_MCDZA=mcolAsRleList(X2.S61.X2.NF54.MCDZA.K4.C_SPMR_FE_comp.bdg,"score"),rep3_MC=mcolAsRleList(X3.S3.X3.NF54.MC.K4.A_SPMR_FE_comp.bdg,"score"),rep3_MCDZA=mcolAsRleList(X3.S7.X3.NF54.MCDZA.K4.A_SPMR_FE_comp.bdg,"score"))
K4.MCDZA.pairs$rep1K4delta<-round(K4.MCDZA.pairs$rep1_MCDZA-K4.MCDZA.pairs$rep1_MC,4)
K4.MCDZA.pairs$rep2K4delta<-round(K4.MCDZA.pairs$rep2_MCDZA-K4.MCDZA.pairs$rep2_MC,4)
K4.MCDZA.pairs$rep3K4delta<-round(K4.MCDZA.pairs$rep3_MCDZA-K4.MCDZA.pairs$rep3_MC,4)
K4.MCDZA.pairs<-K4.MCDZA.pairs[seqnames(K4.MCDZA.pairs) %in% paste0("Pf3D7_",sprintf("%02.0f",1:14), "_v3")]
seqlevels(K4.MCDZA.pairs)<- paste0("Pf3D7_",sprintf("%02.0f",1:14), "_v3")
genome(K4.MCDZA.pairs)<-"3D7"

save(K4.MC.pairs,K9.MC.pairs,K4.MCDZA.pairs,K9.MCDZA.pairs,file="K4me3.K9me3.coverage.RData")

###Quantify % difference in coverage
K9mergedPeaks<- data.frame(import.bed("./K9.mergedPeaks.bed"))
K9mergedPeaks<-sort(makeGRangesFromDataFrame(K9mergedPeaks,keep.extra.columns = T))

K4mergedPeaks<- data.frame(import.bed("./byTreatAbPeaks/K4.mergedPeaks.bed"))
K4mergedPeaks<-sort(makeGRangesFromDataFrame(K4mergedPeaks,keep.extra.columns = T))

options(ucscChromosomeNames=FALSE)

K9peaks.MC.pairs<-subsetByOverlaps(K9.MC.pairs,K9mergedPeaks)
K9peaks.MCDZA.pairs<-subsetByOverlaps(K9.MCDZA.pairs,K9mergedPeaks)
K4peaks.MC.pairs<-subsetByOverlaps(K4.MC.pairs,K4mergedPeaks)
K4peaks.MCDZA.pairs<-subsetByOverlaps(K4.MCDZA.pairs,K4mergedPeaks)


plotme<-
  data.frame(
  mod=c(rep("H3K4me3",5),rep("H3K9me3",5)),
  treat=c(rep("DZA",3),rep("MC",2),rep("DZA",3),rep("MC",2)),
  rep=c(LETTERS[1:3], LETTERS[1:2], LETTERS[1:3], LETTERS[1:2]),
  pctdiff=c(
          sum(K4peaks.MCDZA.pairs$rep1K4delta,na.rm=T)/sum(K4peaks.MCDZA.pairs$rep1_MC,na.rm=T),
          sum(K4peaks.MCDZA.pairs$rep2K4delta,na.rm=T)/sum(K4peaks.MCDZA.pairs$rep2_MC,na.rm=T),
          sum(K4peaks.MCDZA.pairs$rep3K4delta,na.rm=T)/sum(K4peaks.MCDZA.pairs$rep3_MC,na.rm=T),
          sum(K4peaks.MC.pairs$rep1K4delta,na.rm=T)/sum(K4peaks.MCDZA.pairs$rep1_MC,na.rm=T),
          sum(K4peaks.MC.pairs$rep2K4delta,na.rm=T)/sum(K4peaks.MCDZA.pairs$rep2_MC,na.rm=T),
          sum(K9peaks.MCDZA.pairs$rep1K9delta,na.rm=T)/sum(K9peaks.MCDZA.pairs$rep1_MC,na.rm=T),
          sum(K9peaks.MCDZA.pairs$rep2K9delta,na.rm=T)/sum(K9peaks.MCDZA.pairs$rep2_MC,na.rm=T),
          sum(K9peaks.MCDZA.pairs$rep3K9delta,na.rm=T)/sum(K9peaks.MCDZA.pairs$rep3_MC,na.rm=T),
          sum(K9peaks.MC.pairs$rep1K9delta,na.rm=T)/sum(K9peaks.MCDZA.pairs$rep1_MC,na.rm=T),
          sum(K9peaks.MC.pairs$rep2K9delta,na.rm=T)/sum(K9peaks.MCDZA.pairs$rep2_MC,na.rm=T))
  )

ddply(plotme,.(treat,mod),summarize,mean=mean(pctdiff)*100)

ggplot(plotme,aes(x=treat,y=pctdiff*100,fill=mod))+
  stat_summary(fun = "mean", geom = "bar", position="dodge")+
  stat_summary(fun.data= mean_cl_boot, geom="errorbar", position=position_dodge(0.95,preserve = "single"),width=0.3)+
  geom_point(position=position_dodge(0.95))+
  geom_hline(yintercept = 0,color="black")+
  scale_fill_manual(values=c("blue","red"))+
  scale_y_continuous(breaks=c(10,-10,-20,-30,-40))+
  theme_classic()+
  labs(y="% difference in coverage")+
  theme(
    axis.text.y = element_text(color="black"),
    axis.line.x = element_line(size=0.75),
    axis.line.y = element_line(size=0.75), 
    axis.ticks.length = unit(.35,"cm")
  )

###Quantify % difference in coverage for ap2g
plotme1<-
  data.frame(
    mod=c(rep("H3K4me3",5),rep("H3K9me3",5)),
    treat=c(rep("DZA",3),rep("MC",2),rep("DZA",3),rep("MC",2)),
    rep=c(LETTERS[1:3], LETTERS[1:2], LETTERS[1:3], LETTERS[1:2]),
    pctdiff=c(
      sum(K4.MCDZA.pairs$rep1K4delta,na.rm=T)/sum(K4.MCDZA.pairs$rep1_MC,na.rm=T),
      sum(K4.MCDZA.pairs$rep2K4delta,na.rm=T)/sum(K4.MCDZA.pairs$rep2_MC,na.rm=T),
      sum(K4.MCDZA.pairs$rep3K4delta,na.rm=T)/sum(K4.MCDZA.pairs$rep3_MC,na.rm=T),
      sum(K4.MC.pairs$rep1K4delta,na.rm=T)/sum(K4.MC.pairs$rep1_MC,na.rm=T),
      sum(K4.MC.pairs$rep2K4delta,na.rm=T)/sum(K4.MC.pairs$rep2_MC,na.rm=T),
      sum(K9.MCDZA.pairs$rep1K9delta,na.rm=T)/sum(K9.MCDZA.pairs$rep1_MC,na.rm=T),
      sum(K9.MCDZA.pairs$rep2K9delta,na.rm=T)/sum(K9.MCDZA.pairs$rep2_MC,na.rm=T),
      sum(K9.MCDZA.pairs$rep3K9delta,na.rm=T)/sum(K9.MCDZA.pairs$rep3_MC,na.rm=T),
      sum(K9.MC.pairs$rep1K9delta,na.rm=T)/sum(K9.MC.pairs$rep1_MC,na.rm=T),
      sum(K9.MC.pairs$rep2K9delta,na.rm=T)/sum(K9.MC.pairs$rep2_MC,na.rm=T))
  )

ggplot(plotme1,aes(x=treat,y=pctdiff*100,fill=mod))+
  stat_summary(fun = "mean", geom = "bar", position="dodge")+
  stat_summary(fun.data= mean_cl_boot, geom="errorbar", position=position_dodge(0.95,preserve = "single"),width=0.3)+
  geom_point(position=position_dodge(0.95))+
  geom_hline(yintercept = 0,color="black")+
  scale_fill_manual(values=c("black","red"))+
  theme_classic()+
  labs(y="% difference in coverage")+
  theme(
    axis.text.y = element_text(color="black"),
    axis.line.x = element_line(size=0.75),
    axis.line.y = element_line(size=0.75), 
    axis.ticks.length = unit(.35,"cm")
  )

#### ap2-g locus
ap2g<-K9mergedPeaks[52]
K4ap2g<-K4mergedPeaks[1822:1823]
K9.MC.pairs.AP2G<- subsetByOverlaps(K9.MC.pairs,ap2g)
K4.MC.pairs.AP2G<-subsetByOverlaps(K4.MC.pairs,K4ap2g)
K9.MCDZA.pairs.AP2G<-subsetByOverlaps(K9.MCDZA.pairs,ap2g)
K4.MCDZA.pairs.AP2G<-subsetByOverlaps(K4.MCDZA.pairs,K4ap2g)

plotme2<-
  data.frame(
    mod=c(rep("H3K4me3",5),rep("H3K9me3",5)),
    treat=c(rep("DZA",3),rep("MC",2),rep("DZA",3),rep("MC",2)),
    rep=c(LETTERS[1:3], LETTERS[1:2], LETTERS[1:3], LETTERS[1:2]),
    pctdiff=c(
      sum(K4.MCDZA.pairs.AP2G$rep1K4delta,na.rm=T)/sum(K4.MCDZA.pairs.AP2G$rep1_MC,na.rm=T),
      sum(K4.MCDZA.pairs.AP2G$rep2K4delta,na.rm=T)/sum(K4.MCDZA.pairs.AP2G$rep2_MC,na.rm=T),
      sum(K4.MCDZA.pairs.AP2G$rep3K4delta,na.rm=T)/sum(K4.MCDZA.pairs.AP2G$rep3_MC,na.rm=T),
      sum(K4.MC.pairs.AP2G$rep1K4delta,na.rm=T)/sum(K4.MC.pairs.AP2G$rep1_MC,na.rm=T),
      sum(K4.MC.pairs.AP2G$rep2K4delta,na.rm=T)/sum(K4.MC.pairs.AP2G$rep2_MC,na.rm=T),
      sum(K9.MCDZA.pairs.AP2G$rep1K9delta,na.rm=T)/sum(K9.MCDZA.pairs.AP2G$rep1_MC,na.rm=T),
      sum(K9.MCDZA.pairs.AP2G$rep2K9delta,na.rm=T)/sum(K9.MCDZA.pairs.AP2G$rep2_MC,na.rm=T),
      sum(K9.MCDZA.pairs.AP2G$rep3K9delta,na.rm=T)/sum(K9.MCDZA.pairs.AP2G$rep3_MC,na.rm=T),
      sum(K9.MC.pairs.AP2G$rep1K9delta,na.rm=T)/sum(K9.MC.pairs.AP2G$rep1_MC,na.rm=T),
      sum(K9.MC.pairs.AP2G$rep2K9delta,na.rm=T)/sum(K9.MC.pairs.AP2G$rep2_MC,na.rm=T))
  )

ggplot(plotme2,aes(x=treat,y=pctdiff*100,fill=mod))+
  stat_summary(fun = "mean", geom = "bar", position="dodge")+
  stat_summary(fun.data= mean_cl_boot, geom="errorbar", position=position_dodge(0.95,preserve = "single"),width=0.3)+
  geom_point(position=position_dodge(0.95))+
  geom_hline(yintercept = 0,color="black")+
  scale_fill_manual(values=c("black","red"))+
  theme_classic()+
  labs(y="% difference in coverage")+
  theme(
        axis.text.y = element_text(color="black"),
        axis.line.x = element_line(size=0.75),
        axis.line.y = element_line(size=0.75), 
        axis.ticks.length = unit(.35,"cm")
  )

#write.table(data.frame(bindAsGRanges(mcolAsRleList(K4.MCDZA.pairs,"rep1K4delta")))[,c("seqnames","start","end","V1")],file=paste0("./byTreatAbPeaks/deltaK4.MCDZA_rep1.bdg"), quote=F, sep="\t", col.names = F, row.names = F)
#write.table(data.frame(bindAsGRanges(mcolAsRleList(K4.MCDZA.pairs,"rep2K4delta")))[,c("seqnames","start","end","V1")],file=paste0("./byTreatAbPeaks/deltaK4.MCDZA_rep2.bdg"), quote=F, sep="\t", col.names = F, row.names = F)
#write.table(data.frame(bindAsGRanges(mcolAsRleList(K4.MCDZA.pairs,"rep3K4delta")))[,c("seqnames","start","end","V1")],file=paste0("./byTreatAbPeaks/deltaK4.MCDZA_rep3.bdg"), quote=F, sep="\t", col.names = F, row.names = F)
rm(X3.S47.X2.NF54.MC.K4.B_SPMR_FE_comp.bdg,X2.S60.X2.NF54.MC.K4.C_SPMR_FE_comp.bdg,X3.S3.X3.NF54.MC.K4.A_SPMR_FE_comp.bdg,X2.S47.X2.NF54.MCDZA.K4.B_SPMR_FE_comp.bdg,X2.S61.X2.NF54.MCDZA.K4.C_SPMR_FE_comp.bdg,X3.S7.X3.NF54.MCDZA.K4.A_SPMR_FE_comp.bdg)

#save(list=c("K9.MC.pairs", "K4.MC.pairs","K9.MCDZA.pairs","K4.MCDZA.pairs"), file="K4.K9.pairs_deltas.RData")
load(file="K4.K9.pairs_deltas.RData")

#####  Analysis  of differential enrichment in R #####
plotme<-
  data.frame(
    mod=c(rep("H3K4me3",5),rep("H3K9me3",5)),
    treat=c(rep("DZA",3),rep("MC",2),rep("DZA",3),rep("MC",2)),
    rep=c(LETTERS[1:3], LETTERS[1:2], LETTERS[1:3], LETTERS[1:2]),
    pctdiff=c(
      sum(K4peaks.MCDZA.pairs$rep1_MCDZA,na.rm=T)/sum(K4peaks.MCDZA.pairs$rep1_MC,na.rm=T)-1,
      sum(K4peaks.MCDZA.pairs$rep2_MCDZA,na.rm=T)/sum(K4peaks.MCDZA.pairs$rep2_MC,na.rm=T)-1,
      sum(K4peaks.MCDZA.pairs$rep3_MCDZA,na.rm=T)/sum(K4peaks.MCDZA.pairs$rep3_MC,na.rm=T)-1,
      sum(K4peaks.MC.pairs$rep1_,na.rm=T)/sum(K4peaks.MC.pairs$rep1_MC,na.rm=T)-1,
      sum(K4peaks.MC.pairs$rep2_,na.rm=T)/sum(K4peaks.MC.pairs$rep2_MC,na.rm=T)-1,
      sum(K9peaks.MCDZA.pairs$rep1_MCDZA,na.rm=T)/sum(K9peaks.MCDZA.pairs$rep1_MC,na.rm=T)-1,
      sum(K9peaks.MCDZA.pairs$rep2_MCDZA,na.rm=T)/sum(K9peaks.MCDZA.pairs$rep2_MC,na.rm=T)-1,
      sum(K9peaks.MCDZA.pairs$rep3_MCDZA,na.rm=T)/sum(K9peaks.MCDZA.pairs$rep3_MC,na.rm=T)-1,
      sum(K9peaks.MC.pairs$rep1_,na.rm=T)/sum(K9peaks.MC.pairs$rep1_MC,na.rm=T)-1,
      sum(K9peaks.MC.pairs$rep2_,na.rm=T)/sum(K9peaks.MC.pairs$rep2_MC,na.rm=T)-1)
  )

ggplot(plotme[plotme$treat=="MC",],aes(x=treat,y=pctdiff*100,fill=mod))+
  stat_summary(fun = "mean", geom = "bar", position="dodge")+
  geom_hline(yintercept = 0,color="black")+
  scale_fill_manual(values=c("blue","red"))+
  scale_y_continuous(breaks=c(10,5,0,-5,-10,-15,-20,-25,-30), limits=c(-25,15))+
  theme_classic()+
  labs(x="",y="")+
  theme(
    axis.text.y = element_text(color="black"),
    axis.line.x = element_line(size=0.75),
    axis.line.y = element_line(size=0.75), 
    axis.ticks.length = unit(.35,"cm")
  )

#genome-wide H3K9me3 _ vs MC
test<-t.test(c(data.frame(K9peaks.MC.pairs)$rep1_,data.frame(K9peaks.MC.pairs)$rep2_),c(data.frame(K9peaks.MC.pairs)$rep1_MC,data.frame(K9peaks.MC.pairs)$rep2_MC),paired=F)
test
round(test$estimate[1]/test$estimate[2]*100-100,1)
#genome-wide H3K4me3 _ vs MC
test<-t.test(c(data.frame(K4peaks.MC.pairs)$rep1_,data.frame(K4peaks.MC.pairs)$rep2_),c(data.frame(K4peaks.MC.pairs)$rep1_MC,data.frame(K4peaks.MC.pairs)$rep2_MC),paired=F)
test
round(test$estimate[1]/test$estimate[2]*100-100,1)

ggplot(plotme[plotme$treat=="DZA",],aes(x=treat,y=pctdiff*100,fill=mod))+
  stat_summary(fun = "mean", geom = "bar", position="dodge")+
  geom_hline(yintercept = 0,color="black")+
  scale_fill_manual(values=c("blue","red"))+
  scale_y_continuous(breaks=c(10,5,0,-5,-10,-15,-20,-25,-30), limits=c(-45,0))+
  theme_classic()+
  labs(x="",y="")+
  theme(
    axis.text.y = element_text(color="black"),
    axis.line.x = element_line(size=0.75),
    axis.line.y = element_line(size=0.75), 
    axis.ticks.length = unit(.35,"cm")
  )

#genome-wide H3K9me3 MCDZA vs MC
test<-t.test(c(data.frame(K9peaks.MCDZA.pairs)$rep1_MCDZA,data.frame(K9peaks.MCDZA.pairs)$rep2_MCDZA,data.frame(K9peaks.MCDZA.pairs)$rep3_MCDZA),c(data.frame(K9peaks.MCDZA.pairs)$rep1_MC,data.frame(K9peaks.MCDZA.pairs)$rep2_MC,data.frame(K9peaks.MCDZA.pairs)$rep3_MC),paired=F)
test
round(test$estimate[1]/test$estimate[2]*100-100,1)
#genome-wide H3K4me3 MCDZA vs MC
test<-t.test(c(data.frame(K4peaks.MCDZA.pairs)$rep1_MCDZA,data.frame(K4peaks.MCDZA.pairs)$rep2_MCDZA,data.frame(K4peaks.MCDZA.pairs)$rep3_MCDZA),c(data.frame(K4peaks.MCDZA.pairs)$rep1_MC,data.frame(K4peaks.MCDZA.pairs)$rep2_MC,data.frame(K4peaks.MCDZA.pairs)$rep3_MC),paired=F)
test
round(test$estimate[1]/test$estimate[2]*100-100,1)

# mean of ((sum of rep1_)/(sum of rep1_MC) and (sum of rep2_)/(sum of rep_2MC))
# 
# mean(c(sum(data.frame(K9peaks.MCDZA.pairs)[,7],na.rm=T)/sum(data.frame(K9peaks.MCDZA.pairs)[,6],na.rm=T),
#        sum(data.frame(K9peaks.MCDZA.pairs)[,9],na.rm=T)/sum(data.frame(K9peaks.MCDZA.pairs)[,8],na.rm=T),
#        sum(data.frame(K9peaks.MCDZA.pairs)[,11],na.rm=T)/sum(data.frame(K9peaks.MCDZA.pairs)[,10],na.rm=T)))

# plotme<-rbind(
#   data.frame(treat="DZA", pos="K9", genome=colSums(data.frame(K9.MCDZA.pairs)[,12:14],na.rm=T), peaks=colSums(data.frame(subsetByOverlaps(K9.MCDZA.pairs,K9mergedPeaks))[,12:14]), peak.base=colSums(data.frame(subsetByOverlaps(K9.MCDZA.pairs,K9mergedPeaks))[,12:14])/K9.length),
#   data.frame(treat="DZA", pos="K4",genome=colSums(data.frame(K4.MCDZA.pairs)[,12:14],na.rm=T), peaks=colSums(data.frame(subsetByOverlaps(K4.MCDZA.pairs,K4mergedPeaks))[,12:14]), peak.base=colSums(data.frame(subsetByOverlaps(K4.MCDZA.pairs,K4mergedPeaks))[,12:14])/K4.length),
#   data.frame(treat="MC", pos="K9", genome=colSums(data.frame(K9.MC.pairs)[,10:11],na.rm=T), peaks=colSums(data.frame(subsetByOverlaps(K9.MC.pairs,K9mergedPeaks))[,10:11]), peak.base=colSums(data.frame(subsetByOverlaps(K9.MC.pairs,K9mergedPeaks))[,10:11])/K9.length),
#   data.frame(treat="MC", pos="K4", genome=colSums(data.frame(K4.MC.pairs)[,12:14],na.rm=T), peaks=colSums(data.frame(subsetByOverlaps(K4.MC.pairs,K4mergedPeaks))[,12:14]), peak.base=colSums(data.frame(subsetByOverlaps(K4.MC.pairs,K4mergedPeaks))[,12:14])/K4.length)
# )
# plotme$treat<-factor(plotme$treat,levels=c("MC","DZA"))

# ggplot(plotme, aes(x=treat,y=peak.base,fill=pos))+
#   scale_fill_manual(values=c("blue","red"))+
#   stat_summary(fun = "mean", geom = "bar", position="dodge")+
#   stat_summary(fun.data= mean_cl_boot, geom="errorbar", position=position_dodge(0.95,preserve = "single"),width=0.3)+
#   theme_classic()

#### ap2-g locus
ap2g<-K9mergedPeaks[52]
#K4ap2g<-K4mergedPeaks[1822:1823] #??

K9.MC.pairs.AP2G<- subsetByOverlaps(K9.MC.pairs,K9mergedPeaks[52])
K4.MC.pairs.AP2G<-subsetByOverlaps(K4.MC.pairs,K9mergedPeaks[52])
K9.MCDZA.pairs.AP2G<-subsetByOverlaps(K9.MCDZA.pairs,K9mergedPeaks[52])
K4.MCDZA.pairs.AP2G<-subsetByOverlaps(K4.MCDZA.pairs,K9mergedPeaks[52])

plotme<-
  data.frame(
    mod=c(rep("H3K4me3",5),rep("H3K9me3",5)),
    treat=c(rep("DZA",3),rep("MC",2),rep("DZA",3),rep("MC",2)),
    rep=c(LETTERS[1:3], LETTERS[1:2], LETTERS[1:3], LETTERS[1:2]),
    pctdiff=c(
      sum(K9.MC.pairs.AP2G$rep1_,na.rm=T)/sum(K9.MC.pairs.AP2G$rep1_MC,na.rm=T)-1,
      sum(K9.MC.pairs.AP2G$rep2_,na.rm=T)/sum(K9.MC.pairs.AP2G$rep2_MC,na.rm=T)-1,
      sum(K4.MCDZA.pairs.AP2G$rep1_MCDZA,na.rm=T)/sum(K4.MCDZA.pairs.AP2G$rep1_MC,na.rm=T)-1,
      sum(K4.MCDZA.pairs.AP2G$rep2_MCDZA,na.rm=T)/sum(K4.MCDZA.pairs.AP2G$rep2_MC,na.rm=T)-1,
      sum(K4.MCDZA.pairs.AP2G$rep3_MCDZA,na.rm=T)/sum(K4.MCDZA.pairs.AP2G$rep3_MC,na.rm=T)-1,
      sum(K4.MC.pairs.AP2G$rep1_,na.rm=T)/sum(K4.MC.pairs.AP2G$rep1_MC,na.rm=T)-1,
      sum(K4.MC.pairs.AP2G$rep2_,na.rm=T)/sum(K4.MC.pairs.AP2G$rep2_MC,na.rm=T)-1,
      sum(K9.MCDZA.pairs.AP2G$rep1_MCDZA,na.rm=T)/sum(K9.MCDZA.pairs.AP2G$rep1_MC,na.rm=T)-1,
      sum(K9.MCDZA.pairs.AP2G$rep2_MCDZA,na.rm=T)/sum(K9.MCDZA.pairs.AP2G$rep2_MC,na.rm=T)-1,
      sum(K9.MCDZA.pairs.AP2G$rep3_MCDZA,na.rm=T)/sum(K9.MCDZA.pairs.AP2G$rep3_MC,na.rm=T)-1,

  ))

ggplot(plotme[plotme$treat=="MC",],aes(x=treat,y=pctdiff*100,fill=mod,label=pctdiff*100))+
  stat_summary(fun = "mean", geom = "bar", position="dodge")+
  geom_hline(yintercept = 0,color="black")+
  scale_fill_manual(values=c("blue","red"))+
  scale_y_continuous(breaks=c(0,-10,-20,-30,-40), limits=c(-50,0))+
  theme_classic()+
  labs(x="",y="")+
  theme(
    axis.text.y = element_text(color="black"),
    axis.line.x = element_line(size=0.75),
    axis.line.y = element_line(size=0.75), 
    axis.ticks.length = unit(.35,"cm")
  )

#ap2-g locus H3K9me3 _ vs MC
test<-t.test(c(data.frame(K9.MC.pairs.AP2G)$rep1_,data.frame(K9.MC.pairs.AP2G)$rep2_),c(data.frame(K9.MC.pairs.AP2G)$rep1_MC,data.frame(K9.MC.pairs.AP2G)$rep2_MC),paired=F)
test
round(test$estimate[1]/test$estimate[2]*100-100,1)
#ap2-g locus H3K4me3 _ vs MC
test<-t.test(c(data.frame(K4.MC.pairs.AP2G)$rep1_,data.frame(K4.MC.pairs.AP2G)$rep2_),c(data.frame(K4.MC.pairs.AP2G)$rep1_MC,data.frame(K4.MC.pairs.AP2G)$rep2_MC),paired=F)
test
round(test$estimate[1]/test$estimate[2]*100-100,1)

ggplot(plotme[plotme$treat=="DZA",],aes(x=treat,y=pctdiff*100,fill=mod,col=mod,label=pctdiff*100))+
  stat_summary(fun = "mean", geom = "bar", position="dodge")+
  geom_hline(yintercept = 0,color="black")+
  scale_fill_manual(values=c("blue","red"))+
  geom_text()+
  scale_y_continuous(breaks=c(0,-10,-20,-30,-40), limits=c(-50,20))+
  theme_classic()+
  labs(x="",y="")+
  theme(
    axis.text.y = element_text(color="black"),
    axis.line.x = element_line(size=0.75),
    axis.line.y = element_line(size=0.75), 
    axis.ticks.length = unit(.35,"cm")
  )

#ap2-g locus H3K9me3 MCDZA vs MC
test<-t.test(c(data.frame(K9.MCDZA.pairs.AP2G)$rep1_MCDZA,data.frame(K9.MCDZA.pairs.AP2G)$rep2_MCDZA,data.frame(K9.MCDZA.pairs.AP2G)$rep3_MCDZA),c(data.frame(K9.MCDZA.pairs.AP2G)$rep1_MC,data.frame(K9.MCDZA.pairs.AP2G)$rep2_MC,data.frame(K9.MCDZA.pairs.AP2G)$rep3_MC),paired=F)
test
round(test$estimate[1]/test$estimate[2]*100-100,1)
#ap2-g locus H3K4me3 MCDZA vs MC
test<-t.test(c(data.frame(K4.MCDZA.pairs.AP2G)$rep1_MCDZA,data.frame(K4.MCDZA.pairs.AP2G)$rep2_MCDZA,data.frame(K4.MCDZA.pairs.AP2G)$rep3_MCDZA),c(data.frame(K4.MCDZA.pairs.AP2G)$rep1_MC,data.frame(K4.MCDZA.pairs.AP2G)$rep2_MC,data.frame(K4.MCDZA.pairs.AP2G)$rep3_MC),paired=F)
test
round(test$estimate[1]/test$estimate[2]*100-100,1)
