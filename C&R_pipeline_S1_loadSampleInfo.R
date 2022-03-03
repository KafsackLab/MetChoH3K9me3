setwd("/mnt/raid5/ChIPseqArchive/ChoMet/")
require(Rsamtools);require(plyr);
library(reshape2); library(dplyr); library(ggplot2); library(ggsignif); library(ggpubr);library(Rmisc);library(tidyr)
library(InPAS); library(BSGenome.Pf3D7.PlasmoDB.37); library(remotes); library(GenomicRanges)

####Generating necessary folders for processed samples (if doesn't already exist)
if(!dir.exists("./bams/"))dir.create("./bams")
if(!dir.exists("./paired.reads/"))dir.create("./paired.reads")
if(!dir.exists("./singleSamplePeaks/"))dir.create("./singleSamplePeaks")
if(!dir.exists("./singleSamplePeaks/tdfs/"))dir.create("./singleSamplePeaks/tdfs/")
if(!dir.exists("./singleSamplePeaks/SPMR"))dir.create("./singleSamplePeaks/SPMR")
if(!dir.exists("./byTreatAbPeaks/"))dir.create("./byTreatAbPeaks")
if(!dir.exists("./byTreatAbPeaks/tdfs/"))dir.create("./byTreatAbPeaks/tdfs/")
if(!dir.exists("./byExpTreatAbPeaks/"))dir.create("./byExpTreatAbPeaks")
if(!dir.exists("./byExpTreatAbPeaks/tdfs/"))dir.create("./byExpTreatAbPeaks/tdfs/")


#### Load sample info into R#### 
SampleInfo<-read.delim("AllSamples.tab")
SampleInfo$sampleID<-paste0("S",SampleInfo$sampleID)
SampleInfo$ID<-paste0(SampleInfo$exp,".",SampleInfo$sampleID)

SampleInfo[SampleInfo$exp=="X2" & SampleInfo$rep=="B",]
SampleInfo[SampleInfo$exp=="X3" & SampleInfo$rep=="X2B",] # X3.S46 == X2.S45, the other three were not sequenced in X2

FQ.files<-data.frame(rbind(
  data.frame(exp="X1",sampleID=NA,read=1,old.filename=paste0("../20200821_CutandRun/reads/",list.files(path="../20200821_CutandRun/reads/",pattern=".fastq.gz"))),
  data.frame(exp="X2",sampleID=NA,read=1,old.filename=paste0("../20201202_CutandRun/reads/",list.files(path="../20201202_CutandRun/reads/",pattern=".fastq.gz"))),
  data.frame(exp="X3",sampleID=NA,read=1,old.filename=paste0("../20210412_CutandRun/reads/",list.files(path="../20210412_CutandRun/reads/",pattern=".fastq.gz")))
))

FQ.files$read<-substring(FQ.files$old.filename,regexpr("_R\\d_",FQ.files$old.filename,perl=TRUE)+1,regexpr("_R\\d_",FQ.files$old.filename,perl=TRUE)+2)
FQ.files$sampleID<-substring(FQ.files$old.filename,regexpr("_S\\d+_",FQ.files$old.filename,perl=TRUE)+1,regexpr("_S\\d+_",FQ.files$old.filename,perl=TRUE)+3)
FQ.files$sampleID<-sub("_","",FQ.files$sampleID)
FQ.files<-FQ.files[FQ.files$sampleID!="S0",]
FQ.files$ID<-paste0(FQ.files$exp,".",FQ.files$sampleID)
FQ.files<-FQ.files[FQ.files$ID %in% SampleInfo$ID,]

setdiff(SampleInfo$ID,FQ.files$ID)
setdiff(FQ.files$ID,SampleInfo$ID)

FQ.files.r1<-FQ.files[FQ.files$read=="R1",]
FQ.files.r2<-FQ.files[FQ.files$read=="R2",]	

SampleInfo$r1<-FQ.files.r1$old.filename[match(SampleInfo$ID,FQ.files.r1$ID)]
SampleInfo$r2<-FQ.files.r2$old.filename[match(SampleInfo$ID,FQ.files.r2$ID)]

SampleInfo<-SampleInfo[SampleInfo$qual=="OK",] #remove bad quality libraries

SampleInfo$exp[SampleInfo$rep=="X2B"]<-"X2"  # for exp 2 libraries together with exp 3 
SampleInfo$rep[SampleInfo$rep=="X2B"]<-"B"

RunMe<-SampleInfo[SampleInfo$line =="NF54" & SampleInfo$treat!="LS",]

RunMe$sample<-paste0(RunMe$ID,".",RunMe$exp,".",RunMe$line,".",RunMe$treat,".",RunMe$Ab,".",RunMe$rep)
RunMe[RunMe$sample %in% RunMe$sample[duplicated(RunMe$sample)],]
RunMe$bam<-paste0("./bams/",RunMe$sample,".bam")
