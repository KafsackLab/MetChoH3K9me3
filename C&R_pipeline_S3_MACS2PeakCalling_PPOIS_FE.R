##### Call peaks with macs2: Assemble list of good quality samples#####
# macs2 callpeak -t NF54_plus_MC-DZA_H3K9me3_B_S17.bam -c NF54_plus_MC-DZA_IgG_A_S5.bam -f BAMPE -g 2.3e7 -q 0.05 --keep-dup 1 --broad --nomodel -B --SPMR  --max-gap 500 --outdir ./ -n  NF54_plus_MC-DZA_H3K9me3_B_S17 &
ftable(RunMe$exp,RunMe$treat,RunMe$Ab)
ftable(RunMe[RunMe$qual=="OK",]$exp,RunMe[RunMe$qual=="OK",]$treat,RunMe[RunMe$qual=="OK",]$Ab)
RunMe<-RunMe[RunMe$qual=="OK",]

RunMe$bam[RunMe$ID==RunMe$IgG]
RunMe$bdg<-paste0("./singleSamplePeaks/",RunMe$sample,"_treat_pileup.bdg")

# call significant peaks with MACS2 
SysCommands<-paste(
  "macs2 callpeak",
  "-t",	RunMe$bam[RunMe$Ab!="IgG"], 
  "-c", RunMe$bam[match(RunMe$IgG[RunMe$Ab!="IgG"], RunMe$ID)], #paste(RunMe$bam[RunMe$Ab=="IgG" & RunMe$exp==exp],collapse=" "), 
  "-f BAMPE -B -g 2.3e7 -q 0.05 --nomodel --broad --keep-dup auto --max-gap 500 --outdir ./singleSamplePeaks -n", 
  RunMe$sample[RunMe$Ab!="IgG" ],
  "&")[!file.exists(paste0("./singleSamplePeaks/",RunMe$sample[RunMe$Ab!="IgG"],"_treat_pileup.bdg"))]
 sapply(SysCommands,system,ignore.stdout=T,ignore.stderr=F)
# sapply(SysCommands[grepl("K9",SysCommands,fixed=T)],system,ignore.stdout=T,ignore.stderr=F)

# make tdfs from the MACS2 significant peak tracks
SysCommands<-paste(  
  "/home/kafsacklab/src/igv/igvtools toTDF", 
  paste0("./singleSamplePeaks/",list.files(path = "./singleSamplePeaks/",pattern = "pileup.bdg$")[!grepl("lambda",list.files(path = "./singleSamplePeaks/",pattern = "pileup.bdg$"),fixed=T)]), 
  paste0("./singleSamplePeaks/tdfs/",list.files(path = "./singleSamplePeaks/",pattern = "pileup.bdg$")[!grepl("lambda",list.files(path = "./singleSamplePeaks/",pattern = "pileup.bdg$"),fixed=T)],".tdf"),
  "~/igv/genomes/Pf3D7.genome &")[!file.exists(paste0("./singleSamplePeaks/tdfs/",list.files(path = "./singleSamplePeaks/",pattern = ".bdg$")[!grepl("lambda",list.files(path = "./singleSamplePeaks/",pattern = ".bdg$"),fixed=T)],".tdf"))]
# sapply(SysCommands,system,ignore.stdout=T,ignore.stderr=F)


###Generate fold change peaks using -SPMR setting in MACS2 (MACS will save signal per million  reads  for  fragment  pileup  profiles)
RunMe$SPMRbdg<-paste0("./singleSamplePeaks/SPMR/",RunMe$sample,"_SPMR_treat_pileup.bdg")

SysCommands<-paste(
  "macs2 callpeak",
  "-t",	RunMe$bam[RunMe$Ab!="IgG"], 
  "-c", RunMe$bam[match(RunMe$IgG[RunMe$Ab!="IgG"], RunMe$ID)], #paste(RunMe$bam[RunMe$Ab=="IgG" & RunMe$exp==exp],collapse=" "), 
  "-f BAMPE -B -g 2.3e7 -q 0.05 --nomodel --SPMR --broad --keep-dup auto --max-gap 500 --outdir ./singleSamplePeaks/SPMR -n", 
  paste0(RunMe$sample[RunMe$Ab!="IgG" ],"_SPMR"),
  "&")[!file.exists(paste0("./singleSamplePeaks/SPMR/",RunMe$sample[RunMe$Ab!="IgG"],"_SPMR_treat_pileup.bdg"))]
# sapply(SysCommands,system,ignore.stdout=T,ignore.stderr=F)
sapply(SysCommands[grepl("K9",SysCommands,fixed=T)],system,ignore.stdout=T,ignore.stderr=F)

# make fold change tdfs
SysCommands<-paste(  
  "/home/kafsacklab/src/igv/igvtools toTDF", 
  paste0("./singleSamplePeaks/SPMR/",list.files(path = "./singleSamplePeaks/SPMR/",pattern = "pileup.bdg$")[!grepl("lambda",list.files(path = "./singleSamplePeaks/SPMR/",pattern = "pileup.bdg$"),fixed=T)]), 
  paste0("./singleSamplePeaks/SPMR/tdfs/",list.files(path = "./singleSamplePeaks/SPMR/",pattern = "pileup.bdg$")[!grepl("lambda",list.files(path = "./singleSamplePeaks/SPMR/",pattern = "pileup.bdg$"),fixed=T)],".tdf"),
  "~/igv/genomes/Pf3D7.genome &")[!file.exists(paste0("./singleSamplePeaks/SPMR/tdfs/",list.files(path = "./singleSamplePeaks/",pattern = ".bdg$")[!grepl("lambda",list.files(path = "./singleSamplePeaks/SPMR/",pattern = ".bdg$"),fixed=T)],".tdf"))]
# sapply(SysCommands,system, ignore.stdout=T,ignore.stderr=T, wait=F)
# sapply(SysCommands[grepl("K9",SysCommands,fixed=T)],system,ignore.stdout=T,ignore.stderr=F)


#####  Calculate significant and fold enrichments by comparing each treat BDG to each background ####

# Calculate significance enrichment tracks (Tells you how significant enrichment is compared to control_lambda at each position)
RunMe$cbdg<-paste0("./singleSamplePeaks/",RunMe$sample,"_ppois_comp.bdg")

SysCommands<-paste(
  "macs2 bdgcmp",
  "-t", RunMe$bdg[RunMe$Ab!="IgG"][file.exists(RunMe$bdg[RunMe$Ab!="IgG"])],
  "-c", paste0("./singleSamplePeaks/",RunMe$sample[RunMe$Ab!="IgG"][file.exists(RunMe$bdg[RunMe$Ab!="IgG"])],"_control_lambda.bdg"),
  "-m ppois",
  "-o", paste0("./singleSamplePeaks/",RunMe$sample[RunMe$Ab!="IgG"],"_ppois_comp.bdg"),
  "&")[!file.exists(paste0("./singleSamplePeaks/",RunMe$sample[RunMe$Ab!="IgG"],"_ppois_comp.bdg"))]
# sapply(SysCommands,system, ignore.stdout=T,ignore.stderr=T, wait=F)
# sapply(SysCommands[grepl("K9",SysCommands)],system,ignore.stdout=F,ignore.stderr=T, wait=F)

#convert bdg to TDF
SysCommands<-paste(  
  "/home/kafsacklab/src/igv/igvtools toTDF", 
  paste0("./singleSamplePeaks/",list.files(path = "./singleSamplePeaks/",pattern = ".bdg$")[!grepl("lambda",list.files(path = "./singleSamplePeaks/",pattern = ".bdg$"),fixed=T)]), 
  paste0("./singleSamplePeaks/tdfs/",list.files(path = "./singleSamplePeaks/",pattern = ".bdg$")[!grepl("lambda",list.files(path = "./singleSamplePeaks/",pattern = ".bdg$"),fixed=T)],".tdf"),
  "~/igv/genomes/Pf3D7.genome")[!file.exists(paste0("./singleSamplePeaks/tdfs/",list.files(path = "./singleSamplePeaks/",pattern = ".bdg$")[!grepl("lambda",list.files(path = "./singleSamplePeaks/",pattern = ".bdg$"),fixed=T)],".tdf"))]
# sapply(SysCommands[grepl("comp",SysCommands,fixed=T)],system,ignore.stdout=F,ignore.stderr=T, wait=F)
# sapply(SysCommands[grepl("K9",SysCommands,fixed=T) & grepl("comp",SysCommands,fixed=T)],system,ignore.stdout=F,ignore.stderr=T, wait=F)


# Calculate Fold enrichment tracks (Tells you fold change of enrichment compared to control_lambda at each position)
RunMe$SPMRcbdg<-paste0("./singleSamplePeaks/SPMR/",RunMe$sample,"_SPMR_FE_comp.bdg")

SysCommands<-paste(
  "macs2 bdgcmp",
  "-t", RunMe$SPMRbdg[RunMe$Ab!="IgG"][file.exists(RunMe$SPMRbdg[RunMe$Ab!="IgG"])],
  "-c", paste0("./singleSamplePeaks/SPMR/",RunMe$sample[RunMe$Ab!="IgG"],"_SPMR_control_lambda.bdg"),
  "-m FE",
  "-o", RunMe$SPMRcbdg[RunMe$Ab!="IgG"], #paste0("./singleSamplePeaks/SPMR/",RunMe$sample[RunMe$Ab!="IgG"],"_SPMR_FE_comp.bdg"),
  "&")[!file.exists(paste0("./singleSamplePeaks/SPMR/",RunMe$sample[RunMe$Ab!="IgG"],"_SPMR_FE_comp.bdg"))]
# sapply(SysCommands,system,ignore.stdout=F,ignore.stderr=T, wait=F)
# sapply(SysCommands[grepl("K9",SysCommands,fixed=T)],system,ignore.stdout=F,ignore.stderr=T, wait=F)

# make combined FE tdfs
SysCommands<-paste(
  "/home/kafsacklab/src/igv/igvtools toTDF",
  paste0("./singleSamplePeaks/SPMR/",list.files(path = "./singleSamplePeaks/SPMR/",pattern = "comp.bdg$")[!grepl("lambda",list.files(path = "./singleSamplePeaks/SPMR/",pattern = "comp.bdg$"),fixed=T)]),
  paste0("./singleSamplePeaks/SPMR/tdfs/",list.files(path = "./singleSamplePeaks/SPMR/",pattern = "comp.bdg$")[!grepl("lambda",list.files(path = "./singleSamplePeaks/SPMR/",pattern = "comp.bdg$"),fixed=T)],".tdf"),
  "~/igv/genomes/Pf3D7.genome")[!file.exists(paste0("./singleSamplePeaks/SPMR/tdfs/",list.files(path = "./singleSamplePeaks/SPMR/",pattern = ".bdg$")[!grepl("lambda",list.files(path = "./singleSamplePeaks/SPMR/",pattern = ".bdg$"),fixed=T)],".tdf"))]
# sapply(SysCommands[grepl("comp",SysCommands,fixed=T)],system,ignore.stdout=F,ignore.stderr=T, wait=F)
# sapply(SysCommands[grepl("comp",SysCommands,fixed=T) & grepl("K9",SysCommands,fixed=T)],system,ignore.stdout=F,ignore.stderr=T, wait=F)
