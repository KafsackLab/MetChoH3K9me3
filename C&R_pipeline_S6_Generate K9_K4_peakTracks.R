#########   Merging Peaks 	############

# call K9 peaks from replicated combined bedgraphs
#		1. call peaks on combined replicates
#		2. find regions that are called in all 5 conditions
# 	3. Find peaks in each conditions that overlap with these

# 1. call peaks on combined replicates
file.remove(list.files("./byTreatAbPeaks/","bed",full.names = T))
SysCommands<-paste(
  "macs2 bdgbroadcall",
  "-i", paste0("./byTreatAbPeaks/",names(by.treat)[grepl("K9",names(by.treat))],".ppois.bdg"),
  "-c 10 -C 9 -l 2000 -G 500 --o-prefix",
  paste0("./byTreatAbPeaks/",names(by.treat)[grepl("K9",names(by.treat))],"."),
  "&")
sapply(SysCommands,system)

file.rename(list.files(path = "./byTreatAbPeaks/",pattern = ".bed12",full.names = T),sub("bed12","bed",list.files(path = "./byTreatAbPeaks/",pattern = ".bed12",full.names = T)))

SysCommands<-paste(
  "macs2 bdgbroadcall",
  "-i", paste0("./byTreatAbPeaks/",names(by.treat)[grepl("K4",names(by.treat))],".ppois.bdg"),
  "-c 7 -C 6 --o-prefix",
  paste0("./byTreatAbPeaks/",names(by.treat)[grepl("K4",names(by.treat))],"."),
  "&")
sapply(SysCommands,system)

file.rename(list.files(path = "./byTreatAbPeaks/",pattern = ".bed12",full.names = T),sub("bed12","bed",list.files(path = "./byTreatAbPeaks/",pattern = ".bed12",full.names = T)))

# 2. find regions in all the conditions
beds<-intersect(intersect(list.files(path = "./byTreatAbPeaks",pattern = "broad.bed$", full.names = T),list.files(path = "./byTreatAbPeaks",pattern = "K9", full.names = T)),list.files(path = "./byTreatAbPeaks",pattern = "NF54", full.names = T))
SysCommands<-paste(
  "multiIntersectBed",
  "-i", paste(beds, collapse = " "),
  #'| grep "1,2,3,4,5"',  # region found in all 5 conditions
  "> ./byTreatAbPeaks/K9.intersect2.bed"
)
system(SysCommands)

system("awk '$4 >= 4' ./byTreatAbPeaks/K9.intersect2.bed > ./byTreatAbPeaks/K9.intersect.bed")

# 3. find peaks that overlap the regions shared by all
SysCommands<-paste(
  "bedtools intersect",
  "-wa -a", beds,
  "-b ./byTreatAbPeaks/K9.intersect.bed",
  "> ", paste0("./byTreatAbPeaks/",names(by.treat)[grepl("K9",names(by.treat))],".overlap.bed")
)
sapply(SysCommands,system)

# 4. merge all those peaks
SysCommands<-paste(
  "multiIntersectBed",
  "-wa -wb -i", paste0("./byTreatAbPeaks/",names(by.treat)[grepl("K9",names(by.treat))],".overlap.bed",collapse = " "),
  "| bedtools merge",
  "-i stdin > ./byTreatAbPeaks/K9.mergedPeaks.bed"
)
system(SysCommands)


###K4
beds<-intersect(intersect(list.files(path = "./byTreatAbPeaks",pattern = "broad.bed$", full.names = T),list.files(path = "./byTreatAbPeaks",pattern = "K4", full.names = T)),list.files(path = "./byTreatAbPeaks",pattern = "NF54", full.names = T))
SysCommands<-paste(
  "multiIntersectBed",
  "-i", paste(beds, collapse = " "),
  #'| grep "1,2,3,4,5"',  # region found in all 5 conditions
  "> ./byTreatAbPeaks/K4.intersect2.bed"
)
system(SysCommands)

system("awk '$4 >= 4' ./byTreatAbPeaks/K4.intersect2.bed > ./byTreatAbPeaks/K4.intersect.bed")

SysCommands<-paste(
  "bedtools intersect",
  "-wa -a", beds,
  "-b ./byTreatAbPeaks/K4.intersect.bed",
  "> ", paste0("./byTreatAbPeaks/",names(by.treat)[grepl("K4",names(by.treat))],".overlap.bed")
)
sapply(SysCommands,system)

####Import comp bedgraphs of K9 and K4 peaks and convert to Grange objects


# 4. merge all those peaks
SysCommands<-paste(
  "multiIntersectBed",
  "-wa -wb -i", paste0("./byTreatAbPeaks/",names(by.treat)[grepl("K4",names(by.treat))],".overlap.bed",collapse = " "),
  "| bedtools merge",
  "-i stdin > ./byTreatAbPeaks/K4.mergedPeaks.bed"
)
system(SysCommands)



