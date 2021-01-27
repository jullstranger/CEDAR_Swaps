
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
perm<- as.numeric(args[ i + 3 ]) # the prc of individuals to permute

set.seed(137)

n_parts<- 10 # in how many parts eQTLs should be analysed 

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

cis_eqtls<-cis_eqtls[  order(cis_eqtls$FDR),  ]
print(paste("for " , tissue , " we are going to analyse " , nrow(cis_eqtls), " cis_eqtls" , sep=""))
#cis_eqtls <- cis_eqtls [1:105 , ]

indis<-colnames(ge)[2:ncol(ge) ]
if (perm > 0 ){
  genop<-geno[ ,]
  n = as.integer ( perm * length(indis) / 100 )
  s<-sample(2:ncol(geno) , n )
  ss <- sample(s , n )
  genop[, s] <-genop[ , ss]
  
  swapped_indis<-data.frame( ID_1 = colnames( geno) [s] , ID_2 = colnames(geno) [ss])
  
  colnames(genop)<-colnames(ge)
  sum(colnames(genop) != colnames(ge) )
  #[1] 0
}


ge_4_genotypes <-data.frame()

to_do <-nrow(cis_eqtls)

prc<-0
prev_prc<-0

for (i in 1 : nrow(cis_eqtls) ){
# for (i in 1 :  100 ) {
  probe<-cis_eqtls$gene[i]
  snp<-cis_eqtls$snps[i]
  
  ge_val <-unlist(ge[ge$id == probe  , 2:ncol(ge)])
  geno_val <-unlist(geno[geno$id  == snp , 2:ncol(geno) ])
  genop_val <-unlist(genop[genop$id  == snp , 2:ncol(genop) ])
  
  i_0 <- which(geno_val == 0)
  i_1 <- which(geno_val == 1)
  i_2 <- which(geno_val == 2)
  
  i_0_p <- which(genop_val == 0)
  i_1_p <- which(genop_val == 1)
  i_2_p <- which(genop_val == 2)
  

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
                   N_BB = length(i_2) ,
                   
                   AA_p = mean( ge_val[ i_0_p ] , na.rm = TRUE) ,
                   AB_p = mean( ge_val[ i_1_p ] , na.rm = TRUE) ,
                   BB_p = mean (ge_val[ i_2_p ] , na.rm = TRUE) ,
                   sigma_AA_p = sqrt(var(ge_val[i_0_p] , na.rm = TRUE) ) ,
                   sigma_AB_p = sqrt(var(ge_val[i_1_p] , na.rm = TRUE) ) ,
                   sigma_BB_p = sqrt(var(ge_val[i_2_p] , na.rm = TRUE) ) ,
                   N_AA_p = length(i_0_p) ,
                   N_AB_p = length(i_1_p) ,
                   N_BB_p = length(i_2_p) 
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
part<-as.integer(nrow(cis_eqtls)/n_parts )
print( paste("we are going to analyse " , nrow(cis_eqtls) ,  " eQTLs per group of " , part , sep="" ) )
to_do<-length(indis)
prev_prc <- 0

for ( i in  1: to_do ){
  indi<-indis[i]
  z_sum = 0
  z_count = 0
  
  z_sum_p= 0
  z_count_p = 0
  
  roww<-data.frame()
  for ( j in 1:nrow(cis_eqtls) ){
    probe<-cis_eqtls$gene[j]
    snp<-cis_eqtls$snps[j]
    
    ge_indi<-ge[ge$id == probe , indi]
    geno_indi<-geno[geno$id == snp , indi]
    genop_indi<-genop[genop$id == snp , indi]
    
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
    
    z_p = NA
    if (! is.na(genop_indi) ) {
      if ( genop_indi == 0 ){
        if (! is.na( ge_4genotype$sigma_AA_p))
          z_p = abs(ge_4genotype$AA_p - ge_indi)/ge_4genotype$sigma_AA_p
      }
      if (genop_indi == 1){
        if (! is.na( ge_4genotype$sigma_AB_p))
          z_p = abs(ge_4genotype$AB_p - ge_indi)/ge_4genotype$sigma_AB_p
      }
      if (genop_indi == 2 ){
        if (! is.na( ge_4genotype$sigma_BB_p))
          z_p = abs(ge_4genotype$BB_p - ge_indi)/ge_4genotype$sigma_BB_p
      }
    }
    
    if (! is.na(z_p)){
      z_sum_p = z_sum_p + z_p
      z_count_p = z_count_p + 1
    }
    
    k<- as.integer (j / part )
    r <- j %% part
    
    if ( k > 0 & r == 0 | j == nrow(cis_eqtls) ){
      z_avr <- z_sum/z_count
      z_avr_p <- z_sum_p/z_count_p
      
      # print( paste("we have finished " ,  k , " part" , sep="") )
      if ( k == 1 ){
        roww<-data.frame( ID = indi , z_avr = z_avr , z_count = z_count , z_avr_p = z_avr_p ,z_count_p = z_count_p )
        
      }else{
        roww<-cbind(roww,  z_avr = z_avr , z_count = z_count , z_avr_p = z_avr_p ,z_count_p = z_count_p)
       
      }
      if (k == n_parts & r > 0 ){
        k = k + 1
      }
      cols<-paste0( c("z_avr_", "z_count_" , "z_avr_p_" ,"z_count_p_"), rep(k, 4) , sep = "" )
      colnames(roww)[(length(roww) -3) : length(roww)]<-cols
    } 
    
  }
  if (nrow(z_indis) == 0 ){
    z_indis <- roww
  }else{
    z_indis<-rbind(z_indis , roww)
  }
  
  prc<- as.integer(i * 100 / to_do)
  if (prc %% 5 == 0 & prc != prev_prc){
    print(paste(prc, "% done" , sep="") )
    prev_prc <- prc
  }
  
} 

if (perm > 0 ){
  
  f2w<-paste( "Z_scores_", tissue , "_perm_" , perm , ".txt" , sep="") 
  write.table(file = f2w , z_indis , row.names = FALSE , quote = FALSE)
  
  swapped_indis<-merge( swapped_indis , z_indis [ , c(1,grep( "avr_p_" , colnames(z_indis)) )], by.x = "ID_1"  , by.y = "ID")
  
  f2w<-paste( "Z_scores_", tissue , "_perm_" , perm , "_swappd_indis.txt" , sep="") 
  write.table(file = f2w , swapped_indis , row.names = FALSE , quote = FALSE)
  
  # f2p<-paste( "Z_scores_", tissue , "_with_perm_" , perm, ".pdf" , sep="") 
  # pdf(f2p)
  # par(mfcol=c(5,2))
  # hist(z_indis$z_avr_1 ,breaks = 30 , xlab = "Z" , main = paste( "Z scores hist " , tissue , " 10%" , sep = "" ) )
  # hist(z_indis$z_avr_2 ,breaks = 30 , xlab = "Z" , main = paste( "Z scores hist " , tissue , " 20%" , sep = "" ) )
  # hist(z_indis$z_avr_3 ,breaks = 30 , xlab = "Z" , main = paste( "Z scores hist " , tissue , " 30%" , sep = "" ) )
  # hist(z_indis$z_avr_4 ,breaks = 30 , xlab = "Z" , main = paste( "Z scores hist " , tissue , " 40%" , sep = "" ) )
  # hist(z_indis$z_avr_5 ,breaks = 30 , xlab = "Z" , main = paste( "Z scores hist " , tissue , " 50%" , sep = "" ) )
  # hist(z_indis$z_avr_6 ,breaks = 30 , xlab = "Z" , main = paste( "Z scores hist " , tissue , " 60%" , sep = "" ) )
  # hist(z_indis$z_avr_7 ,breaks = 30 , xlab = "Z" , main = paste( "Z scores hist " , tissue , " 70%" , sep = "" ) )
  # hist(z_indis$z_avr_8 ,breaks = 30 , xlab = "Z" , main = paste( "Z scores hist " , tissue , " 80%" , sep = "" ) )
  # hist(z_indis$z_avr_9 ,breaks = 30 , xlab = "Z" , main = paste( "Z scores hist " , tissue , " 90%" , sep = "" ) )
  # hist(z_indis$z_avr_10 ,breaks = 30 , xlab = "Z" , main = paste( "Z scores hist " , tissue , " 100%" , sep = "" ) )
  # 
  # dev.off()
  
  ###############################################
  
  f2r<-paste( "Z_scores_", tissue , "_perm_", perm, ".txt" , sep="")
  z_scores<-read.table(f2r, header = TRUE , stringsAsFactors =  FALSE)
  
  n_indis<-nrow(z_scores)
  n_indis_p<-as.integer( n_indis * perm / 100 )
  
  for ( i in 1 : as.integer(n_parts/2) ) {
    
    part<- ( i * 2 - 1 )
    
    z_col_1<-paste("z_avr_" , part, sep="")
    z_col_p_1<-paste("z_avr_p_" , part , sep="")
    prc_1<- 10 * part
    thresh_1<- sort( z_scores [ , z_col_p_1 ])[ n_indis - n_indis_p ]
    swaps_1<-sum(z_scores[, z_col_1] > thresh_1 )
    
    z_col_2<-paste("z_avr_" , part + 1   , sep="")
    z_col_p_2<-paste("z_avr_p_" , part + 1  , sep="")
    prc_2 <- 10 * (part + 1 )
    thresh_2<- sort( z_scores [ , z_col_p_2 ])[ n_indis - n_indis_p ]
    swaps_2<-sum(z_scores[, z_col_2] > thresh_2 )
    
    
    f2p<-paste("Z_scores_", tissue , "_perm_", perm, "_" , part , "_" , part + 1 , ".pdf" , sep="")
    pdf(f2p)
    
    par(mfcol=c(2,2))
    
    hist(z_scores[, z_col_1] ,breaks = 20 , xlab = "Z" , main = paste( "Z scores hist " , tissue ," ", 
                                                                       prc_1, "% eQTLs \n swaps = " , swaps_1 , sep = "" ) )
    abline(v = thresh_1 , col = "red")
    
    hist(z_scores[, z_col_2 ] ,breaks = 20 , xlab = "Z" , main = paste( "Z scores hist " , tissue , " ",
                                                                        prc_2 ,  "% eQTLs \n swaps = " , swaps_2 , sep = "" ) )
    abline(v = thresh_2 , col = "blue")
    hist(z_scores[, z_col_p_1] ,breaks = 20 , xlab = "Z" , main = paste( "Z scores hist perm ",  perm , "% ", 
                                                                         "\n threshold = " , round(thresh_1 , 3 ) ,
                                                                         "\n expected swaps = " , n_indis_p  , sep = "" ) )
    abline(v = thresh_1 , col = "red")
    
    hist(z_scores[, z_col_p_2 ] ,breaks = 20 , xlab = "Z" , main = paste( "Z scores hist  perm " ,  perm , "% ", 
                                                                          "\n threshold = " , round(thresh_2 , 3 )  , sep = "" ) )
    abline(v = thresh_2 , col = "blue")
    
    dev.off()
    
  }
  
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


###################################################
#z_score<-function(s){
# mu <- mean(s , na.rm =  TRUE)
# sigma <-sqrt (var( s , na.rm =  TRUE) )
# return ( mean(abs( s - mu ) / sigma  , na.rm = TRUE) )
#}

#z_val<-apply(ge[, 2:ncol(ge)], 2 , z_score)

###################################################

###################################################
#c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
#c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
#hist_perm <- hist(z_indis$z_avr , breaks = 30 , main = paste ("Z scores Perm" , tissue , sep = "") )
#hist_real <- hist(z_scores$z_avr , breaks = 30 , main = paste("Z scores " , tissue , sep = "") )
#par(mfrow=c(1,1))
#plot(hist_real , col = c1 )
#plot(hist_perm , col = c2 , add = TRUE)
