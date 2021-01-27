args= commandArgs()
i = 1
while(args[i]!="--args"){
  cat(args[i])
  i=i+1
}

tissue <- args[ i + 1 ]
work.dir <- args[ i + 2 ]
perm<- as.numeric(args[ i + 3 ])

work.dir<-"/media/julia/e4a91c72-ba4c-44e5-a062-4bba142849b4/CEDAR_GE_GIGA/"
workdir<-paste( work.dir , "/" , tissue , sep="")
setwd(workdir)
print(getwd())

SNP_file_name = paste( "SNPs4", tissue , ".txt" , sep="")
#SNP_location_file_name = paste( work.dir , "SNPs_location.txt" , sep="")
expression_file_name = paste("GE_" , tissue , ".txt" , sep ="")
#gene_location_file_name = paste(work.dir ,  "Gene_location.txt", sep="")
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

indis<-colnames(ge)[2:ncol(ge) ]
if (perm > 0 ){
  n_indis<-ncol(geno) -1
  perm_n <- as.integer(perm * n_indis /100)
  s<-sample(2:ncol(geno) , perm_n)
  ss<-sample(s, perm_n)
  geno[ , s]<-geno[, ss ]
  
  colnames(geno)<-colnames(ge)
  sum(colnames(geno) != colnames(ge) )
  #[1] 0
}

ge_4_genotypes <-data.frame()

to_do <-nrow(cis_eqtls)

prc<-0
prev_prc<-0

for (i in 1 : nrow(cis_eqtls) ){
  probe<-cis_eqtls$gene[i]
  snp<-cis_eqtls$snps[i]
  
  ge_val <-unlist(ge[ge$id == probe  , 2:ncol(ge)])
  geno_val <-unlist(geno[geno$id  == snp , 2:ncol(geno) ])
  
  i_0 <- which(geno_val == 0)
  i_1 <- which(geno_val == 1)
  i_2 <- which(geno_val == 2)
  
  
  roww<-data.frame(SNP = snp , 
                   Probe = probe, 
                   AA = mean( ge_val[ i_0 ] , na.rm = TRUE) ,
                   AB = mean( ge_val[ i_1 ] , na.rm = TRUE) ,
                   BB = mean (ge_val[ i_2 ] , na.rm = TRUE) ,
                   sigma_AA = sqrt(var(ge_val[i_0] , na.rm = TRUE) ) ,
                   sigma_AB = sqrt(var(ge_val[i_1] , na.rm = TRUE) ) ,
                   sigma_BB = sqrt(var(ge_val[i_2] , na.rm = TRUE) ) ,
                   N_AA = length(i_0) ,
                   N_AB = length(i_1) ,
                   N_BB = length(i_2) 
                   
  )
  if (i == 1){
    ge_4_genotypes<-roww
  }else{
    ge_4_genotypes<-rbind(ge_4_genotypes , roww)
  }
  
  prc <- as.integer(i*100 / to_do)
  if (prc %% 5 == 0 && prc != prev_prc){
    print(paste(prc , "% done" , sep=""))
    prev_prc = prc
  }
}

###########
# for each individual calculate zeta score for each probe, calculate average for this individual
z_indis<-data.frame()

for ( i in 1: length(indis) ){
  indi<-indis[i]
  z_sum = 0
  z_count = 0
  
  for ( j in 1:nrow(cis_eqtls) ){
    probe<-cis_eqtls$gene[j]
    snp<-cis_eqtls$snps[j]
    
    ge_indi<-ge[ge$id == probe , indi]
    geno_indi<-geno[geno$id == snp , indi]
    
    ge_4genotype<-ge_4_genotypes[ge_4_genotypes$SNP == snp & ge_4_genotypes$Probe == probe , ]
    
    z = NA
    if (! is.na(geno_indi) ) {
      if ( geno_indi == 0 ){
        if (! is.na( ge_4genotype$sigma_AA))
          z = abs(ge_4genotype$AA - ge_indi)/ge_4genotype$sigma_AA
      }
      if (geno_indi == 1){
        if (! is.na( ge_4genotype$sigma_AB))
          z = abs(ge_4genotype$AB - ge_indi)/ge_4genotype$sigma_AB
      }
      if (geno_indi == 2 ){
        if (! is.na( ge_4genotype$sigma_BB))
          z = abs(ge_4genotype$BB - ge_indi)/ge_4genotype$sigma_BB
      }
    }
    
    if (! is.na(z)){
      z_sum = z_sum + z
      z_count = z_count + 1
    }
    
  }
  
  z_avr <- z_sum/z_count
  
  if (i == 1){
    z_indis<-data.frame( ID = indi , z_avr = z_avr , z_count = z_count)
  }else{
    z_indis<-rbind(z_indis , data.frame( ID = indi , z_avr = z_avr , z_count = z_count ))
  }
  print(i)
} 

if (perm == "TRUE"){
  f2r<-paste( "Z_scores_", tissue , ".txt" , sep="") 
  z_scores<- read.table( f2r , header = TRUE , stringsAsFactors =  FALSE)
  
  f2w<-paste( "Z_scores_", tissue , "_perm_", perm, ".txt" , sep="") 
  write.table(file = f2w , z_indis , row.names = FALSE , quote = FALSE)
  
  
  f2p<-paste( "Z_scores_", tissue , "_with_perm_", perm, ".pdf" , sep="") 
  pdf(f2p)
  par(mfrow=c(2,2))
  hist(z_scores$z_avr,breaks = 30 , xlab = "Z" , main = paste( "Z scores hist " , tissue , sep = "" ) )
  plot( sort(z_scores$z_avr) , ylab = "Z_avr")
  
  hist(z_indis$z_avr , breaks = 30 , xlab = "Z" , main = paste( "Z scores hist PERM " , tissue , sep = "" ) )
  plot( sort(z_indis$z_avr) , ylab = "Z_avr")
  
  dev.off()
  
}else{
  f2w<-paste( "Z_scores_", tissue , ".txt" , sep="") 
  write.table(file = f2w , z_indis , row.names = FALSE , quote = FALSE)
  
  f2p<-paste( "Z_scores_", tissue , ".pdf" , sep="") 
  pdf(f2p)
  par(mfrow=c(1,2))
  hist(z_indis$z_avr, breaks = 30 , main = paste( "Z scores hist " , tissue , sep = "" ) )
  plot( sort(z_indis$z_avr) , ylab = "Z_avr")
  dev.off()
  
}



