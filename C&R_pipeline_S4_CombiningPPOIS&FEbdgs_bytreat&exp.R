#### Combine ppois and FE bdgs replicates within each experiment and treatment type #####

#"MC"    "_"     "LS"    "M"     "MCDZA" "C"
RunMe<-RunMe[RunMe$sample!="X2.S59.X2.NF54.MCDZA.K9.C",]  # wasn't actually treated
#RunMe<-RunMe[RunMe$sample!="X2.S58.X2.NF54.MC.K9.C",]  
RunMe<-RunMe[RunMe$sample!="X3.S14.X3.NF54.MC.K9.B",]
#RunMe<-RunMe[RunMe$sample!="X3.S30.X3.NF54.C.K9.C",]
#RunMe<-RunMe[RunMe$sample!="X3.S46.X2.NF54.MCDZA.K9.B",]

by.exp.treat<-unlist(dlply(RunMe[RunMe$Ab!="IgG",],.(line,exp,Ab,treat),.fun=function(df) paste(df$cbdg,collapse=" ")))
names(by.exp.treat)
SysCommands<-paste(
  "macs2 cmbreps",
  "-i", by.exp.treat,
  "-m fisher",
  "-o", paste0("./byExpTreatAbPeaks/",names(by.exp.treat),".ppois.bdg"),
  "&")[!file.exists(paste0("./byExpTreatAbPeaks/",names(by.exp.treat),".ppois.bdg"))]
# sapply(SysCommands,system,ignore.stdout=T,ignore.stderr=F)
# sapply(SysCommands[grepl("MCDZA",SysCommands,fixed=T)],system,ignore.stdout=F,ignore.stderr=T, wait=F)

by.exp.treat<-unlist(dlply(RunMe[RunMe$Ab!="IgG",],.(line,exp,Ab,treat),.fun=function(df) paste(df$SPMRcbdg,collapse=" ")))
names(by.exp.treat)
SysCommands<-paste(
  "macs2 cmbreps",
  "-i", by.exp.treat,
  "-m mean",
  "-o", paste0("./byExpTreatAbPeaks/",names(by.exp.treat),".FE.mean.bdg"),
  "&")[!file.exists(paste0("./byExpTreatAbPeaks/",names(by.exp.treat),".FE.mean.bdg"))]
# sapply(SysCommands,system,ignore.stdout=T,ignore.stderr=F)
# sapply(SysCommands[grepl("K9",SysCommands,fixed=T)],system,ignore.stdout=F,ignore.stderr=T, wait=F)

#convert bdg to TDF
SysCommands<-paste(
  "/home/kafsacklab/src/igv/igvtools toTDF",
  paste0("./byExpTreatAbPeaks/",list.files(path="./byExpTreatAbPeaks/",pattern=".bdg$")),
  paste0("./byExpTreatAbPeaks/tdfs/",list.files(path="./byExpTreatAbPeaks/",pattern=".bdg$"),".tdf"),
  "~/igv/genomes/Pf3D7.genome &")[!file.exists(paste0("./byExpTreatAbPeaks/tdfs/",list.files(path="./byExpTreatAbPeaks/",pattern=".bdg$"),".tdf"))]
# sapply(SysCommands,system,ignore.stdout=T,ignore.stderr=F)


##### Combine ppois and FE bdgs tracks by treatment type across all experiments and replicates#####
by.treat<-unlist(dlply(RunMe[RunMe$Ab!="IgG",],.(line,Ab,treat),.fun=function(df) paste(df$cbdg,collapse=" ")))
# by.treat<-unlist(dlply(RunMe[RunMe$Ab!="IgG",],.(line,Ab,treat),.fun=function(df) paste(df$bdg,collapse=" ")))
names(by.treat)
SysCommands<-paste(
  "macs2 cmbreps",
  "-i", by.treat,
  "-m fisher",
  "-o", paste0("./byTreatAbPeaks/",names(by.treat),".ppois.bdg"),
  "&")[!file.exists(paste0("./byTreatAbPeaks/",names(by.treat),".ppois.bdg"))]
# sapply(SysCommands,system,ignore.stdout=T,ignore.stderr=F)
# sapply(SysCommands[grepl("K9",SysCommands, fixed=T)],system,ignore.stdout=T,ignore.stderr=F)

by.treat<-unlist(dlply(RunMe[RunMe$Ab!="IgG",],.(line,Ab,treat),.fun=function(df) paste(df$SPMRcbdg,collapse=" ")))
names(by.treat)
SysCommands<-paste(
  "macs2 cmbreps",
  "-i", by.treat,
  "-m mean",
  "-o", paste0("./byTreatAbPeaks/",names(by.treat),".FE.mean.bdg"),
  "&")[!file.exists(paste0("./byTreatAbPeaks/",names(by.treat),".FE.mean.bdg"))]
# sapply(SysCommands,system,ignore.stdout=T,ignore.stderr=F)
# sapply(SysCommands[grepl("K9",SysCommands,fixed=T)],system,ignore.stdout=F,ignore.stderr=T, wait=F)

#convert bdg to TDF
SysCommands<-paste( 
  "/home/kafsacklab/src/igv/igvtools toTDF", 
  paste0("./byTreatAbPeaks/",list.files(path="./byTreatAbPeaks/",pattern="bdg$")), 
  paste0("./byTreatAbPeaks/tdfs/",list.files(path="./byTreatAbPeaks/",pattern="bdg$"),".tdf"), 
  "~/igv/genomes/Pf3D7.genome &")[!file.exists(paste0("./byTreatAbPeaks/tdfs/",list.files(path="./byTreatAbPeaks/",pattern="bdg$"),".tdf"))]
# sapply(SysCommands,system,ignore.stdout=T,ignore.stderr=F)
# sapply(SysCommands[grepl("K9",SysCommands,fixed=T)& grepl("FE",SysCommands,fixed=T)],system,ignore.stdout=F,ignore.stderr=T, wait=F)
