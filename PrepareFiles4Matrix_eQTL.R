args= commandArgs()
i = 1
while(args[i]!="--args"){
  cat(args[i])
  i=i+1
}

tissue <- args[ i + 1 ]
workdir <- args[ i + 2 ]
# workdir<-"/media/julia/e4a91c72-ba4c-44e5-a062-4bba142849b4/CEDAR_GE_GIGA/"

geno.file <- paste (workdir , "CEDAR_Genotypes.traw" , sep = "")
geno<-read.table(geno.file , header = TRUE , stringsAsFactors = FALSE)
nrow(geno)
# 694769

unique(geno$CHR)
# subset for autosomal only

geno<-geno[geno$CHR %in% c(1:22), ]
nrow(geno)
# 676800

snp.loc<-data.frame(SNP = geno$SNP , chr = paste( "chr" , geno$CHR, sep="")  , pos = geno$POS)
geno<-data.frame(id = geno$SNP , geno[7:ncol(geno)])

samples<-colnames(geno)[2:ncol(geno)]
sampless<- unlist(strsplit(samples , "_"))[seq(1, length(samples) * 2  , by =2)]
colnames(geno)<-c('id' , sampless)
#####################################""
ge.file<-paste(workdir , tissue , "_Log2_RSN_Batch_Not_Corr_4Plink.txt" , sep="")
ge<-read.table(ge.file , header = TRUE , stringsAsFactors = FALSE , check.names =  FALSE)
probes<-colnames(ge)[3:ncol(ge)]
indis<-ge$FID

ge<-t(ge[, 3:ncol(ge)])
colnames(ge) <- indis
ge<-data.frame(id = probes ,ge)

common.indis<-intersect(indis , sampless)
# those are the samples to include to the analysis
# subset ge for common indis and subset genotypes for common indis
ge.cols<-c("id" , common.indis)
ge<-ge[ , ge.cols]

geno<-geno[ , ge.cols]
print("chould be zero and this is :")
sum(colnames(geno) != colnames(ge) )
#OK
###### check SNPs########################

sum(geno$id != snp.loc$SNP)
# OK
########################### PROBES

probes.file<-paste(workdir , "Probes_good_reanno_31137_TSS.txt" , sep="")
probes<-read.table(probes.file , header = TRUE , stringsAsFactors = FALSE)
# subset
probes<-probes[probes$ProbeID %in% ge[,1],  ]
# order
probes<-probes[ order (probes$CHR , probes$From) , ]
ge.ord<-ge[as.character(probes$ProbeID) , ]

sum(rownames(ge.ord) != as.character(probes$ProbeID))
gene.loc<- data.frame(geneid = probes$ProbeID , chr = probes$Chr , s1 = probes$From , s2 = probes$To)
#[1] 0
##################################################
#### write genotypes
geno[1:5 , 1:5 ]
outDir<- paste(workdir , tissue , sep="")

if (!file.exists(outDir))
{
  dir.create(outDir)
  setwd(outDir)
}
print(getwd() )

f2w<-paste( "SNPs4" , tissue , ".txt" , sep="")
write.table(file = f2w, geno , row.names = FALSE, quote = FALSE)

#### write SNP location
snp.loc [1:5 ,]
f2w<-paste( "SNPs_location.txt" , sep="")
write.table(file = f2w, snp.loc , row.names = FALSE, quote = FALSE)
#### write GE 
ge.ord [1:5 , 1:5]
f2w<-paste( "GE_" , tissue , ".txt" , sep="")
write.table(file = f2w, ge.ord , row.names = FALSE, quote = FALSE)
#### write gene location
gene.loc[1:5 , ]
f2w<-paste( "Gene_location_" , tissue , ".txt" , sep="")
write.table(file = f2w, gene.loc , row.names = FALSE, quote = FALSE)

