setwd("/mnt/raid5/ChIPseqArchive/ChoMet/")
###DataTracks
options(ucscChromosomeNames=FALSE)
dataROI <- DataTrack(K9.MC.pairs,ucscChromosomeNames=FALSE,chromosome = "Pf3D7_12_v3") 
plotTracks(dataROI,type="l")

#####
K9MCtrackrep2 <- bindAsGRanges(rep2_MC=mcolAsRleList(K9.MC.pairs,"rep2_MC"))
dataROI <- DataTrack(K9MCtrackrep2,ucscChromosomeNames=FALSE,chromosome = "Pf3D7_12_v3") 
plotTracks(c(dataROI,genomeAxis),type="histogram",col.histogram="red",window=-1, windowSize = 5000)

K9_trackrep2 <- bindAsGRanges(rep2_=mcolAsRleList(K9.MC.pairs,"rep2_"))
dataROI <- DataTrack(K9_trackrep2,ucscChromosomeNames=FALSE,chromosome = "Pf3D7_12_v3") 
plotTracks(c(dataROI,genomeAxis),type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(0,6))

K9MCDZAtrackrep2_MC <- bindAsGRanges(rep1_MC=mcolAsRleList(K9.MCDZA.pairs,"rep2_MC"))
dataROI <- DataTrack(K9MCDZAtrackrep2_MC,ucscChromosomeNames=FALSE,chromosome = "Pf3D7_12_v3") 
plotTracks(c(dataROI,genomeAxis),type="histogram",col.histogram="red",window=-1, windowSize = 5000)

K9MCDZAtrackrep2_MCDZA <- bindAsGRanges(rep1_MCDZA=mcolAsRleList(K9.MCDZA.pairs,"rep1_MCDZA"))
dataROI <- DataTrack(K9MCDZAtrackrep2_MCDZA,ucscChromosomeNames=FALSE,chromosome = "Pf3D7_12_v3") 
plotTracks(c(dataROI,genomeAxis),type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(0,8))


K4MCtrackrep1 <- bindAsGRanges(rep1_MC=mcolAsRleList(K4.MC.pairs,"rep1_MC"))
dataROI <- DataTrack(K4MCtrackrep1,ucscChromosomeNames=FALSE,chromosome = "Pf3D7_12_v3") 
plotTracks(c(dataROI,genomeAxis),type="histogram",col.histogram="black",window=-1, windowSize = 5000, ylim=c(0,6))

K4_trackrep1 <- bindAsGRanges(rep1_=mcolAsRleList(K4.MC.pairs,"rep1_"))
dataROI <- DataTrack(K4_trackrep1,ucscChromosomeNames=FALSE,chromosome = "Pf3D7_12_v3") 
plotTracks(c(dataROI,genomeAxis),type="histogram",col.histogram="black",window=-1, windowSize = 5000, ylim=c(0,6))

K4MCDZAtrackrep2_MC <- bindAsGRanges(rep2_MC=mcolAsRleList(K4.MCDZA.pairs,"rep2_MC"))
dataROI <- DataTrack(K4MCDZAtrackrep2_MC,ucscChromosomeNames=FALSE,chromosome = "Pf3D7_12_v3") 
plotTracks(c(dataROI,genomeAxis),type="histogram",col.histogram="black",window=-1, windowSize = 5000,ylim=c(0,8))

K4MCDZAtrackrep2_MCDZA <- bindAsGRanges(rep2_MCDZA=mcolAsRleList(K4.MCDZA.pairs,"rep2_MCDZA"))
dataROI <- DataTrack(K4MCDZAtrackrep2_MCDZA,ucscChromosomeNames=FALSE,chromosome = "Pf3D7_12_v3") 
plotTracks(c(dataROI,genomeAxis),type="histogram",col.histogram="black",window=-1, windowSize = 5000,ylim=c(0,8))

#### deltaMC
options(ucscChromosomeNames=FALSE)
K9mergedPeaks<- data.frame(import.bed("./K9.mergedPeaks.bed"))
K9mergedPeaks<- droplevels(K9mergedPeaks,exclude="Pf_M76611")
K9mergedPeaks<- K9mergedPeaks[1:63,]
levels(K9mergedPeaks$seqnames)<-c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10","Chr11","Chr12","Chr13","Chr14")
K9mergedPeaks<-sort(makeGRangesFromDataFrame(K9mergedPeaks,keep.extra.columns = T))
rep2K9deltaMC<- bindAsGRanges(rep2_deltaMC=mcolAsRleList(K9.MC.pairs,"rep2K9delta"))
rep2K4deltaMC <- bindAsGRanges(rep2_deltaMC=mcolAsRleList(K4.MC.pairs,"rep2K4delta"))
seqlevels(rep2K9deltaMC)<-c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10","Chr11","Chr12","Chr13","Chr14")
seqlevels(rep2K4deltaMC)<-c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10","Chr11","Chr12","Chr13","Chr14")

#chromosome1
rep2K9deltaMC_chrom1 <- DataTrack(rep2K9deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr01",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-3,0.5)) 
rep2K4deltaMC_chrom1 <- DataTrack(rep2K4deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr01",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-3,0.5),size=2)
chrom1AnnoT<- AnnotationTrack(K9mergedPeaks,genome="51",chromosome= "Chr01",name="Chr01", fill="red")
plotTracks(c(rep2K4deltaMC_chrom1,rep2K9deltaMC_chrom1,chrom1AnnoT))

#chromosome2
rep2K9deltaMC_chrom2 <- DataTrack(rep2K9deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr02",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-1.5,0.5)) 
rep2K4deltaMC_chrom2 <- DataTrack(rep2K4deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr02",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-1.5,0.5),size=2)
chrom2AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr02",name="Chr02", fill="red")
plotTracks(c(rep2K4deltaMC_chrom2,rep2K9deltaMC_chrom2,chrom2AnnoT))


###chromosome3
rep2K9deltaMC_chrom3 <- DataTrack(rep2K9deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr03",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-1.5,0.5)) 
rep2K4deltaMC_chrom3 <- DataTrack(rep2K4deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr03",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-1.5,0.5),size=2)
chrom3AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr03",name="Chr03", fill="red")
plotTracks(c(rep2K4deltaMC_chrom3,rep2K9deltaMC_chrom3,chrom3AnnoT))

###chromosome4
rep2K9deltaMC_chrom4 <- DataTrack(rep2K9deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr04",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-1.5,0.5)) 
rep2K4deltaMC_chrom4 <- DataTrack(rep2K4deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr04",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-1.5,0.5),size=2)
chrom4AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr04",name="Chr04", fill="red")
plotTracks(c(rep2K4deltaMC_chrom4,rep2K9deltaMC_chrom4,chrom4AnnoT))

###chromosome5
rep2K9deltaMC_chrom5 <- DataTrack(rep2K9deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr05",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-1.7,0.5)) 
rep2K4deltaMC_chrom5 <- DataTrack(rep2K4deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr05",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-1.7,0.5),size=2)
chrom5AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr05",name="Chr05", fill="red")
plotTracks(c(rep2K4deltaMC_chrom5,rep2K9deltaMC_chrom5,chrom5AnnoT))

###chromosome6
rep2K9deltaMC_chrom6 <- DataTrack(rep2K9deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr06",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-2.3,0.5)) 
rep2K4deltaMC_chrom6 <- DataTrack(rep2K4deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr06",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-2.3,0.5),size=2)
chrom6AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr06",name="Chr06", fill="red")
plotTracks(c(rep2K4deltaMC_chrom6,rep2K9deltaMC_chrom6,chrom6AnnoT))

###chromosome7
rep2K9deltaMC_chrom7 <- DataTrack(rep2K9deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr07",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-2.3,0.5)) 
rep2K4deltaMC_chrom7 <- DataTrack(rep2K4deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr07",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-2.3,0.5),size=2)
chrom7AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr07",name="Chr07", fill="red")
plotTracks(c(rep2K4deltaMC_chrom7,rep2K9deltaMC_chrom7,chrom7AnnoT))

###chromosome8
rep2K9deltaMC_chrom8 <- DataTrack(rep2K9deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr08",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-1.8,0.5)) 
rep2K4deltaMC_chrom8 <- DataTrack(rep2K4deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr08",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-1.8,0.5),size=2)
chrom8AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr08",name="Pf3D7_08", fill="red")
plotTracks(c(rep2K4deltaMC_chrom8,rep2K9deltaMC_chrom8,chrom8AnnoT))

###chromosome9
rep2K9deltaMC_chrom9 <- DataTrack(rep2K9deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr09",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-2.5,0.5)) 
rep2K4deltaMC_chrom9 <- DataTrack(rep2K4deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr09",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-2.5,0.5),size=2)
chrom9AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr09",name="Chr09", fill="red")
plotTracks(c(rep2K4deltaMC_chrom9,rep2K9deltaMC_chrom9,chrom9AnnoT))

###chromosome10
rep2K9deltaMC_chrom10 <- DataTrack(rep2K9deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr10",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-2,0.5)) 
rep2K4deltaMC_chrom10 <- DataTrack(rep2K4deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr10",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-2,0.5),size=2)
chrom10AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr10",name="Chr10", fill="red")
plotTracks(c(rep2K4deltaMC_chrom10,rep2K9deltaMC_chrom10,chrom10AnnoT))

###chromosome11
rep2K9deltaMC_chrom11 <- DataTrack(rep2K9deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr11",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-1.2,0.5)) 
rep2K4deltaMC_chrom11 <- DataTrack(rep2K4deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr11",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-1.2,0.5),size=2)
chrom11AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr11",name="Chr11", fill="red")
plotTracks(c(rep2K4deltaMC_chrom11,rep2K9deltaMC_chrom11,chrom11AnnoT))

###chromosome12
rep2K9deltaMC_chrom12 <- DataTrack(rep2K9deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr12",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-2.7,0.5)) 
rep2K4deltaMC_chrom12 <- DataTrack(rep2K4deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr12",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-2.7,0.5),size=2)
chrom12AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr12",name="Pf3D7_12", fill="red")
#plotTracks(c(rep2K9deltaMC_chrom12,rep2K4deltaMC_chrom12,chrom12AnnoT,CDS.track,genomeAxis), from=900778,to=950599)
plotTracks(c(rep2K4deltaMC_chrom12,rep2K9deltaMC_chrom12,chrom12AnnoT))

seqlevels(rep2K9deltaMC)<- c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10","Chr11","Chr12","Chr13","Chr14")
seqlevels(rep2K4deltaMC)<- c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10","Chr11","Chr12","Chr13","Chr14")
genome(rep2K9deltaMC)<-"51"
genome(rep2K4deltaMC)<-"51"

rep2K9deltaMC_ch12<-rep2K9deltaMC[rep2K9deltaMC@seqnames == "Chr12",]
plotme<-data.frame(sort(subsetByOverlaps(gff[gff$feature%in% c("CDS","ncRNA"),],rep2K9deltaMC_ch12)))
plotme$STRAND<-"plus"
plotme$STRAND[plotme$strand=="-"]<-"minus"

CDS.track<-GeneRegionTrack(
  plotme,
  chr="Chr12",
  stacking = "pack",
  feature=plotme$STRAND,
  col="black",
  name="genes",
  rotation.title=0,
  transcriptAnnotation = "gene",
  plus="blue",
  minus="red",
  cex.group=0.75,
  cex.title=1,
  just.group="left",
  fontsize.group=10,
  fontcolor.group="black"
)

ap2g<-data.frame(K9mergedPeaks[K9mergedPeaks@ranges== "900778-917599",])
ap2g<-GeneRegionTrack(ap2g,chr="Chr12",fill="red")
genomeAxis <- GenomeAxisTrack(name="Chr12",littleTicks=F,fontsize=12,col="black",showTitle=T,rotation.title=0,fill.range="transparent",fontcolor="black",labelPos="below",cex.title=0.75,col.border.title="white",ticksAt=c(seq(890000, 920000, by=10000)),cex.title=1,font.size=10) 
## motifs :ticksAt= c(890000,892000,894000,896000,898000,900000,902000,903663,903670,903690,904000,904304,904375,904717,904857,904894,906000,908000,910000,915000,920000)
rep2K9deltaMC_chrom12 <- DataTrack(rep2K9deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr12",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-2.7,2.7),name="H3K9me3 \n(Low SAM - High SAM)",fontsize.legend=22) 
rep2K4deltaMC_chrom12 <- DataTrack(rep2K4deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr12",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-2.7,2.7),name="H3K4me3 \n(Low SAM - High SAM)",fontsize.legend=22)
motifs<-HighlightTrack(c(rep2K4deltaMC_chrom12,rep2K9deltaMC_chrom12,ap2g,CDS.track),range=gene,genome="51",start=903663,end=904935,inBackground=F,fill="gray50",col="black",alpha=0.3)
chrom12AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr12",name="Pf3D7_12", fill="red")
plotTracks(c(motifs,genomeAxis), from=887000,to=923000,just.group="below",background.panel = "white", background.title = "black",fontsize=14,sizes=c(5,5,0.3,2,1.5))
pdf(file = "./SAM_PF3D7_1222600_DiffTrack.pdf",width = 8 ,height=4)
dev.off() 

###chromosome13
rep2K9deltaMC_chrom13 <- DataTrack(rep2K9deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr13",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-1.5,0.5)) 
rep2K4deltaMC_chrom13 <- DataTrack(rep2K4deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr13",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-1.5,0.5),size=2)
chrom13AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr13",name="Pf3D7_13", fill="red")
plotTracks(c(rep2K4deltaMC_chrom13,rep2K9deltaMC_chrom13,chrom13AnnoT))

###chromosome14
rep2K9deltaMC_chrom14 <- DataTrack(rep2K9deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr14",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-1.5,0.5)) 
rep2K4deltaMC_chrom14 <- DataTrack(rep2K4deltaMC,ucscChromosomeNames=FALSE,chromosome = "Chr14",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-1.5,0.5),size=2)
chrom14AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr14",name="Chr14", fill="red")
plotTracks(c(rep2K4deltaMC_chrom14,rep2K9deltaMC_chrom14,chrom14AnnoT))

####Delta MCDZA
options(ucscChromosomeNames=FALSE)
rep1K9deltaMCDZA <- bindAsGRanges(rep2_deltaMCDZA=mcolAsRleList(K9.MCDZA.pairs,"rep1K9delta"))
rep2K4deltaMCDZA <- bindAsGRanges(rep2_deltaMC=mcolAsRleList(K4.MCDZA.pairs,"rep2K4delta"))
seqlevels(rep1K9deltaMCDZA)<-c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10","Chr11","Chr12","Chr13","Chr14")
seqlevels(rep2K4deltaMCDZA)<-c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10","Chr11","Chr12","Chr13","Chr14")


#chromosome1
rep1K9deltaMCDZA_chrom1 <- DataTrack(rep1K9deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr01",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-3,0.5)) 
rep2K4deltaMCDZA_chrom1 <- DataTrack(rep2K4deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr01",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-3,0.5),size=2)
chrom1AnnoT<- AnnotationTrack(K9mergedPeaks,genome="51",chromosome= "Chr01",name="Chr01", fill="red")
plotTracks(c(rep2K4deltaMCDZA_chrom1,rep1K9deltaMCDZA_chrom1,chrom1AnnoT))

#chromosome2
rep1K9deltaMCDZA_chrom2 <- DataTrack(rep1K9deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr02",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-2.5,0.5)) 
rep2K4deltaMCDZA_chrom2 <- DataTrack(rep2K4deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr02",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-2.5,0.5),size=2)
chrom2AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr02",name="Chr02", fill="red")
plotTracks(c(rep2K4deltaMCDZA_chrom2,rep1K9deltaMCDZA_chrom2,chrom2AnnoT))


###chromosome3
rep1K9deltaMCDZA_chrom3 <- DataTrack(rep1K9deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr03",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-2.5,0.5)) 
rep2K4deltaMCDZA_chrom3 <- DataTrack(rep2K4deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr03",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-2.5,0.5),size=2)
chrom3AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr03",name="Chr03", fill="red")
plotTracks(c(rep2K4deltaMCDZA_chrom3,rep1K9deltaMCDZA_chrom3,chrom3AnnoT))

###chromosome4
rep1K9deltaMCDZA_chrom4 <- DataTrack(rep1K9deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr04",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-3.5,0.5)) 
rep2K4deltaMCDZA_chrom4 <- DataTrack(rep2K4deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr04",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-3.5,0.5),size=2)
chrom4AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr04",name="Chr04", fill="red")
plotTracks(c(rep2K4deltaMCDZA_chrom4,rep1K9deltaMCDZA_chrom4,chrom4AnnoT))

###chromosome5
rep1K9deltaMCDZA_chrom5 <- DataTrack(rep1K9deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr05",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-2.5,0.5)) 
rep2K4deltaMCDZA_chrom5 <- DataTrack(rep2K4deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr05",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-2.5,0.5),size=2)
chrom5AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr05",name="Chr05", fill="red")
plotTracks(c(rep2K4deltaMCDZA_chrom5,rep1K9deltaMCDZA_chrom5,chrom5AnnoT))

###chromosome6
rep1K9deltaMCDZA_chrom6 <- DataTrack(rep1K9deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr06",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-2.3,0.5)) 
rep2K4deltaMCDZA_chrom6 <- DataTrack(rep2K4deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr06",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-2.3,0.5),size=2)
chrom6AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr06",name="Chr06", fill="red")
plotTracks(c(rep2K4deltaMCDZA_chrom6,rep1K9deltaMCDZA_chrom6,chrom6AnnoT))

###chromosome7
rep1K9deltaMCDZA_chrom7 <- DataTrack(rep1K9deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr07",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-2.8,0.5)) 
rep2K4deltaMCDZA_chrom7 <- DataTrack(rep2K4deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr07",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-2.8,0.5),size=2)
chrom7AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr07",name="Chr07", fill="red")
plotTracks(c(rep2K4deltaMCDZA_chrom7,rep1K9deltaMCDZA_chrom7,chrom7AnnoT))

###chromosome8
rep1K9deltaMCDZA_chrom8 <- DataTrack(rep1K9deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr08",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-3,0.5)) 
rep2K4deltaMCDZA_chrom8 <- DataTrack(rep2K4deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr08",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-3,0.5),size=2)
chrom8AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr08",name="Pf3D7_08", fill="red")
plotTracks(c(rep2K4deltaMCDZA_chrom8,rep1K9deltaMCDZA_chrom8,chrom8AnnoT))

###chromosome9
rep1K9deltaMCDZA_chrom9 <- DataTrack(rep1K9deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr09",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-2.5,0.5)) 
rep2K4deltaMCDZA_chrom9 <- DataTrack(rep2K4deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr09",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-2.5,0.5),size=2)
chrom9AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr09",name="Chr09", fill="red")
plotTracks(c(rep2K4deltaMCDZA_chrom9,rep1K9deltaMCDZA_chrom9,chrom9AnnoT))

###chromosome10
rep1K9deltaMCDZA_chrom10 <- DataTrack(rep1K9deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr10",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-3,0.5)) 
rep2K4deltaMCDZA_chrom10 <- DataTrack(rep2K4deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr10",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-3,0.5),size=2)
chrom10AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr10",name="Chr10", fill="red")
plotTracks(c(rep2K4deltaMCDZA_chrom10,rep1K9deltaMCDZA_chrom10,chrom10AnnoT))

###chromosome11
rep1K9deltaMCDZA_chrom11 <- DataTrack(rep1K9deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr11",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-2.5,0.5)) 
rep2K4deltaMCDZA_chrom11 <- DataTrack(rep2K4deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr11",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-2.5,0.5),size=2)
chrom11AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr11",name="Chr11", fill="red")
plotTracks(c(rep2K4deltaMCDZA_chrom11,rep1K9deltaMCDZA_chrom11,chrom11AnnoT))

###chromosome12
rep1K9deltaMCDZA_chrom12<-rep1K9deltaMCDZA[rep1K9deltaMCDZA@seqnames == "Chr12",]
plotme<-data.frame(sort(subsetByOverlaps(gff[gff$feature%in% c("CDS","ncRNA"),],rep1K9deltaMCDZA_chrom12)))
plotme$STRAND<-"plus"
plotme$STRAND[plotme$strand=="-"]<-"minus"

CDS.track<-GeneRegionTrack(
  plotme,
  chr="Chr12",
  stacking = "pack",
  feature=plotme$STRAND,
  col="black",
  name="genes",
  rotation.title=0,
  transcriptAnnotation = "gene",
  plus="blue",
  minus="red",
  cex.group=0.75,
  cex.title=1,
  just.group="left",
  fontsize.group=10,
  fontcolor.group="black"
)

ap2g<-data.frame(K9mergedPeaks[K9mergedPeaks@ranges== "900778-917599",])
ap2g<-GeneRegionTrack(ap2g,chr="Chr12",fill="red")
genomeAxis <- GenomeAxisTrack(name="Chr12",littleTicks=F,fontsize=12,col="black",showTitle=T,rotation.title=0,fill.range="transparent",fontcolor="black",labelPos="below",cex.title=0.75,col.border.title="white",ticksAt=c(seq(890000, 920000, by=10000)),cex.title=1,font.size=10) 
## motifs :ticksAt= c(890000,892000,894000,896000,898000,900000,902000,903663,903670,903690,904000,904304,904375,904717,904857,904894,906000,908000,910000,915000,920000)
rep1K9deltaMCDZA_chrom12 <- DataTrack(rep1K9deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr12",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-1,1),name="H3K9me3 \n(High SAH - Low SAH)",fontsize.legend=22)
rep2K4deltaMCDZA_chrom12 <- DataTrack(rep2K4deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr12",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-1,1),name="H3K9me3 \n(High SAH - Low SAH)",fontsize.legend=22)
motifs_MCDZA<-HighlightTrack(c(rep2K4deltaMCDZA_chrom12,rep1K9deltaMCDZA_chrom12,ap2g,CDS.track),range=gene,genome="51",start=903663,end=904935,inBackground=F,fill="gray50",col="black",alpha=0.3)
chrom12AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr12",name="Pf3D7_12", fill="red")
plotTracks(c(motifs_MCDZA,genomeAxis), from=887000,to=923000,just.group="below",background.panel = "white", background.title = "black",fontsize=14,sizes=c(5,5,0.3,2,1.5))
#plotTracks(c(rep2K4deltaMCDZA_chrom12,rep1K9deltaMCDZA_chrom12,chrom12AnnoT),just.group="below",background.panel = "white", background.title = "black",fontsize=14)


###chromosome13
rep1K9deltaMCDZA_chrom13 <- DataTrack(rep1K9deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr13",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-2.8,0.5)) 
rep2K4deltaMCDZA_chrom13 <- DataTrack(rep2K4deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr13",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-2.8,0.5),size=2)
chrom13AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr13",name="Pf3D7_13", fill="red")
plotTracks(c(rep2K4deltaMCDZA_chrom13,rep1K9deltaMCDZA_chrom13,chrom13AnnoT))

###chromosome14
rep1K9deltaMCDZA_chrom14 <- DataTrack(rep1K9deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr14",type="histogram",col.histogram="red",window=-1, windowSize = 5000,ylim=c(-3,0.5)) 
rep2K4deltaMCDZA_chrom14 <- DataTrack(rep2K4deltaMCDZA,ucscChromosomeNames=FALSE,chromosome = "Chr14",type="histogram",col.histogram="blue",window=-1, windowSize = 5000,ylim=c(-3,0.5),size=2)
chrom14AnnoT<- AnnotationTrack(K9mergedPeaks,chromosome= "Chr14",name="Chr14", fill="red")
plotTracks(c(rep2K4deltaMCDZA_chrom14,rep1K9deltaMCDZA_chrom14,chrom14AnnoT))

