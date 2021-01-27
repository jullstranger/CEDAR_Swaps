args= commandArgs()
i = 1
while(args[i]!="--args"){
  cat(args[i])
  i=i+1
}

tissue <- args[ i + 1 ]
work.dir <- args[ i + 2 ]

library(MatrixEQTL)

work.dir<-"/media/julia/e4a91c72-ba4c-44e5-a062-4bba142849b4/CEDAR_GE_GIGA/"
workdir<-paste(work.dir , "/" , tissue , sep="")
setwd(workdir)
print( getwd() )

SNP_file_name = paste("SNPs4", tissue , ".txt" , sep="")
SNP_location_file_name = paste( "SNPs_location.txt" , sep="")
expression_file_name = paste( "GE_" , tissue , ".txt" , sep ="")
gene_location_file_name = paste("Gene_location_" , tissue , ".txt", sep="")

pvOutputThreshold = 1e-02

errorCovariance = numeric()
cisDist = 5e5

snps = SlicedData$new()
snps$fileOmitCharacters="NA"
snps$fileDelimiter=" "
snps$fileSliceSize = 2000
snps$LoadFile(SNP_file_name)

## Load gene expression data

gene = SlicedData$new()
gene$fileDelimiter = " "      # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values;
gene$fileSliceSize = 2000      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name)



snpspos = read.table(SNP_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

useModel = modelLINEAR;

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = NULL,
  pvOutputThreshold.cis = pvOutputThreshold,
  pvOutputThreshold     = 0,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

#unlink(output_file_name_tra)
#unlink(output_file_name)

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
#show(me$cis$eqtls)

cat('Detected distant eQTLs:', '\n');
#show(me$trans$eqtls)

plot(me)

# chose only significant cis-eQTLs 
cis_eqtls_sign <- me$cis$eqtls[me$cis$eqtls$FDR< 0.05 , ]
nrow(cis_eqtls_sign)
# [1] 12277
# chose the best eQTL per gene

f2r<-paste ( work.dir , "/Probes_good_reanno_31137_TSS.txt" , sep="")
probes<-read.table(f2r, header = TRUE , stringsAsFactors =  FALSE)

cis_eqtls_sign_anno<-merge(cis_eqtls_sign , probes[ , c("ProbeID" , "Gene")] , by.x = "gene" , by.y = "ProbeID" )

genes<-unique(cis_eqtls_sign_anno$Gene)
length(genes)
#[1] 1501

cis_eqtls_sign_anno_per_gene <-data.frame()
for ( i in 1:length(genes)){
  gene <- genes[i]
  cis_sub<- cis_eqtls_sign_anno[cis_eqtls_sign_anno$Gene == gene , ]
  best<-which.min(cis_sub$FDR)
  if (i == 1 ){
    cis_eqtls_sign_anno_per_gene <- cis_sub[best, ]
  }else{
    cis_eqtls_sign_anno_per_gene <- rbind(cis_eqtls_sign_anno_per_gene , cis_sub[best, ])
  }
}

f2w<-paste(  "Cis_eQTLs_sign_per_gene_", tissue, ".txt" , sep="")
write.table(file = f2w ,cis_eqtls_sign_anno_per_gene , row.names = FALSE , quote = FALSE)

