
#tissue<-"CD15"
#perm<-10
#work.dir<-"/media/julia/e4a91c72-ba4c-44e5-a062-4bba142849b4/CEDAR_GE_GIGA/"

args= commandArgs()
i = 1
while(args[i]!="--args"){
  cat(args[i])
  i=i+1
}

tissue <- args[ i + 1 ]
work.dir <- args[ i + 2 ]
eqtls_prc <-as.numeric(args[ i + 3 ])

workdir<-paste( work.dir , "/" , tissue , sep="")
setwd(workdir)
print(getwd())

SNP_file_name = paste( "SNPs4", tissue , ".txt" , sep="")

expression_file_name = paste("GE_" , tissue , ".txt" , sep ="")

cis_eqtls_file<-paste ( "Cis_eQTLs_sign_per_gene_" , tissue , ".txt" , sep="")

geno<-read.table(SNP_file_name , header = TRUE , stringsAsFactors =  FALSE)
ge<-read.table (expression_file_name , header = TRUE , stringsAsFactors = FALSE)
cis_eqtls<-read.table(cis_eqtls_file , header = TRUE , stringsAsFactors =  FALSE)
# to test 
sum(colnames(geno ) != colnames(ge) )
# [1] 0

# for eqch cis eQTL gene estimate the average gene expression per genotype
# for each individual calculate zet score based on the real gene expression  and predicted gene expression for each gene
# calculate the average scoe per individual

cis_eqtls<-cis_eqtls[  order(cis_eqtls$FDR),  ]

if(eqtls_prc > 0 ){
  n_eqtls<-as.integer(nrow(cis_eqtls)*eqtls_prc/100)
  cis_eqtls <- cis_eqtls[ 1 : n_eqtls , ] # subset cis_eqtls 
}

print(paste("for " , tissue , " we are going to analyse " , nrow(cis_eqtls), " cis_eqtls" , sep=""))

ge_imp<-matrix(nrow=nrow(cis_eqtls) , ncol = (ncol(ge) - 1) )
stds<-matrix(nrow=nrow(cis_eqtls) , ncol = (ncol(ge) - 1) ) # it is not vars, this is std
ge_mat<-matrix(nrow=nrow(cis_eqtls) , ncol = (ncol(ge) - 1) )

indis<-colnames( ge )[2:ncol(ge)]
n_indis <- length(indis)

prc<-0
prev_prc <- 0
to_do <- nrow(cis_eqtls)

for (i in 1 : to_do) {
  snp <- cis_eqtls$snps[i]
  gene <- cis_eqtls$gene[i]
  
  ge_s<-ge[ge$id == gene , 2:ncol(ge)]
  geno_s <- geno[geno$id == snp , 2 : ncol(geno)]
  ge_mat[i , ] <-as.numeric(ge_s)
  
  ge_ss<-rep(NA , n_indis)
  ge_ss[ geno_s == 0 ] <- mean(ge_s[geno_s == 0] , na.rm = TRUE)
  ge_ss[ geno_s == 1] <- mean(ge_s[geno_s == 1] , na.rm = TRUE)
  ge_ss[ geno_s == 2] <- mean(ge_s[geno_s == 2] , na.rm = TRUE)
  
  std_ss<-rep(NA , n_indis)
  std_ss[ geno_s == 0 ] <- sqrt(var(ge_s[geno_s == 0] , na.rm = TRUE) )
  std_ss[ geno_s == 1] <- sqrt(var(ge_s[geno_s == 1] , na.rm = TRUE) )
  std_ss[ geno_s == 2] <- sqrt(var(ge_s[geno_s == 2] , na.rm = TRUE) )
  ge_imp[i , ] <- as.numeric(ge_ss)
  stds[i , ] <- as.numeric(std_ss)

  prc<-as.integer(100 * i / to_do)
  if(prc %% 5 == 0 && prc != prev_prc){
    print (paste(prc , "% done" , sep=""))
    prev_prc <- prc
  }
}

###############
### you have two matrices : imputed gene expression and sigma
### 

ncol(ge_mat)
ncol(ge_imp)
ncol(stds)

z_scores<-matrix(ncol = n_indis , nrow = n_indis)

for( i in 1 : n_indis ) {
  if (i == 1){
    indxs <- 1 : n_indis
  }else{
    indxs<- c(i:n_indis , 1 : (i-1))
  }
  
  z <- abs(ge_mat - ge_imp[ , indxs])/stds[ , indxs]
  z<-colMeans(z , na.rm = TRUE)
  
  for ( j in 1 : length(indxs)){
    z_scores[j , indxs[j] ]<-z[j]
  }
}

rownames(z_scores) <-indis
colnames(z_scores)<-indis

ftp<- paste( "Z_scores_all_vs_all_", tissue , "_" , eqtls_prc , "_prc_eqtls.pdf" , sep="")
pdf(ftp)
image(z_scores)
dev.off()

f2w<- paste( "Z_scores_all_vs_all_", tissue , "_" , eqtls_prc ,  "_prc_eqtls.txt" , sep="")
write.table(file = f2w , z_scores , quote = FALSE)



