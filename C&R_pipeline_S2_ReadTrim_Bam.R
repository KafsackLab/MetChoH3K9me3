#### quality trim reads #####################################################################################################################
# java -jar ~/jar/trimmomatic-0.38.jar PE H3K9me3-x1_R1.fastq.gz H3K9me3-x1_R2.fastq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:/home/kafsacklab/src/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:30 &

#assembling code (in Sample info file) for adapter and quality read trimming
RunMe$pr1<-paste0("./paired.reads/",RunMe$sample,"_R1.paired.fq.gz")
RunMe$upr1<-paste0("./",RunMe$sample,"_R1.unpaired.fastq.gz")
RunMe$pr2<-paste0("./paired.reads/",RunMe$sample,"_R2.paired.fq.gz")
RunMe$upr2<-paste0("./",RunMe$sample,"_R2.unpaired.fastq.gz")

SysCommands<-paste(
  "java -jar ~/jar/trimmomatic-0.38.jar PE",
  RunMe$r1, 
  RunMe$r2, 
  RunMe$pr1, 
  RunMe$upr1, 
  RunMe$pr2, 
  RunMe$upr2,
  "ILLUMINACLIP:/home/kafsacklab/src/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:30 &"
)[!file.exists(paste0("./paired.reads/",RunMe$sample,"_R1.paired.fq.gz"))]
sapply(SysCommands,system,ignore.stdout=T,ignore.stderr=F)
system("rm *fastq.gz",ignore.stdout=T,ignore.stderr=F) # removes symbolic links and unpaired reads fastq.gz files


#### format genome for BWA ####
SysCommands<-"bwa index ./genome/3D7+K12.fasta &"
# 	sapply(SysCommands,system,ignore.stdout=T,ignore.stderr=F)

##### map paired readsmap with BWA ####
#bwa mem -M -t 5 -R "@RG\tID:kafsacklab\tSM:X1A.NF54.MC.K4\tPL:Illumina\tLB:lib\tPU:unit" ./genome/PlasmoDB-46_Pfalciparum3D7_Genome.fasta X1A.NF54.MC.K4_R1.paired.fq.gz X1A.NF54.MC.K4_R2.paired.fq.gz | samtools sort --threads 5 --output-fmt BAM -o X1A.NF54.MC.K4.bam &
samples2map<-RunMe$sample[!file.exists(RunMe$bam)]

SysCommands<-paste0(
  'bwa mem -M -t 5 -R "@RG\\tID:kafsacklab\\tSM:',
  samples2map,
  '\\tPL:Illumina\\tLB:lib\\tPU:unit"',
  ' ./genome/3D7+K12.fasta ',
  "./paired.reads/",samples2map,'_R1.paired.fq.gz ',
  "./paired.reads/",samples2map,'_R2.paired.fq.gz ',
  '| samtools view -hbf 2 - | samtools sort --threads 5 -n -O BAM -o ',
  paste0("./bams/",samples2map,".bam"),'&\n')

system("rm runme.sh") 
cat(c("#!/bin/bash\nsource ~/.profile\n",SysCommands[1]),file="runme.sh")
system('chmod u+x runme.sh')
system("./runme.sh")  # problem was it didn't have the ENV in my ~./profile

# samtools quickcheck ./bams/*.bam && echo 'all ok' || echo 'fail!'"
SysCommands<-paste(
  "samtools quickcheck ./bams/*.bam && echo 'all ok' || echo 'fail!'"
)

system("rm runme.sh") 
cat(c("#!/bin/bash\nsource ~/.profile\n",SysCommands),file="runme.sh")
system('chmod u+x runme.sh')
# system("./runme.sh") 

RunMe$bamsize<-round(file.size(RunMe$bam)/1024^2,3)
# samtools quickcheck ./bams/*.bam && echo 'all ok' || echo 'fail!'
# remove bad sequences files: Look at Bam Files to asses quality, then edit "qual" column in AllSamples.tab
# RunMe<-RunMe[RunMe$qual=="OK",] #remove bad quality libraries
ftable(RunMe$treat,RunMe$Ab)
ftable(RunMe$exp,RunMe$treat,RunMe$Ab)

###### get reads for K12 + mito + apico for each file ###### 
#sortBam(RunMe$bam, paste0("./bams",RunMe$sample,".sorted.bam"))
RunMe$sbam<-paste0("./bams/",RunMe$sample,".sorted.bam")
samples2sort<-RunMe$sample[!file.exists(RunMe$sbam)]
SysCommands<-paste("samtools sort -O BAM --reference ./genome/Pf3D7+K12.fasta -o", paste0("./bams/",samples2sort,".sorted.bam"),paste0("./bams/",samples2sort,".bam"),"&\n")

system("rm runme.sh") 
cat(c("#!/bin/bash\nsource ~/.profile\n",SysCommands),file="runme.sh")
system('chmod u+x runme.sh')
# system("./runme.sh") 

###### index sorted bams ######
RunMe$bai<-paste0("./bams/",RunMe$sample,".sorted.bam.bai")
SysCommands<-paste("samtools index", RunMe$sbam[file.exists(RunMe$bam) & !file.exists(RunMe$bai)], "&")
# sapply(SysCommands,system, ignore.stdout=T,ignore.stderr=F)

###### Get total reads and reads per ecoli mito and apico
idxstats<-function(sample){
  x<-idxstatsBam(RunMe$sbam[RunMe$sample==sample], index=RunMe$bai[RunMe$sample==sample])
  y<-c(ecR=x$mapped[1],apiR=x$mapped[16],mitoR=x$mapped[17],totalR=sum(x$mapped))
  return(y)
}
RunMe<-cbind(RunMe,ldply(lapply(RunMe$sample,idxstats)))


rm(FQ.files,FQ.files.r1,FQ.files.r2, samples2map, samples2sort)