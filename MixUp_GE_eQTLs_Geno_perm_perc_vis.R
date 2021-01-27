
args= commandArgs()
i = 1
while(args[i]!="--args"){
  cat(args[i])
  i=i+1
}

tissue <- args[ i + 1 ]
work.dir <- args[ i + 2 ]
perm<- as.numeric(args[ i + 3 ]) # the prc of individuals to permute

perm<-15
n_parts<- 10 # in how many parts eQTLs should be analysed 
work.dir<-"/media/julia/e4a91c72-ba4c-44e5-a062-4bba142849b4/CEDAR_GE_GIGA/"
tissues<-c("CD4" , "CD8" , "CD14" , "CD15" , "CD19" , "PLA" , "IL" , "TR" , "RE")

swaps_tab <-data.frame()
for (tissue in tissues){
  if (tissue == "CD19")
    break
  
  workdir<-paste( work.dir , "/" , tissue , sep="")
  setwd(workdir)
  print(getwd())
  
  f2r <- paste("Z_scores_", tissue , "_perm_", perm , "_swappd_indis.txt" , sep="")
  swapped_indis<-read.table(f2r, header = TRUE , stringsAsFactors =  FALSE)
  same <-which (swapped_indis$ID_1 == swapped_indis$ID_2)
  if (length(same) > 0 ){
    swapped_indis <-swapped_indis[ -same , ]
    print(paste("removed " , same , " individuals from swapped indis" , sep=""))
  }
  f2r<-paste( "Z_scores_", tissue , "_perm_", perm, ".txt" , sep="")
  z_scores<-read.table(f2r, header = TRUE , stringsAsFactors =  FALSE)
  n_indis<-nrow(z_scores)
  n_indis_p<-as.integer( n_indis * perm / 100 )
  swaps<-c()
  for ( i in 1 : as.integer(n_parts/2) ) {
    
    part<- ( i * 2 - 1 )
    
    z_col_1<-paste("z_avr_" , part, sep="")
    z_col_p_1<-paste("z_avr_p_" , part , sep="")
    prc_1<- 10 * part
    thresh_1<- sort( z_scores [ , z_col_p_1 ])[ n_indis - n_indis_p ]
    swaps_1<-sum(z_scores[, z_col_1] > thresh_1 )
    
    if (length(swaps) == 0 ){
      swaps<- z_scores$ID[z_scores[ , z_col_1] > thresh_1 ] 
    }else{
      swaps<-intersect(swaps , z_scores$ID[z_scores[ , z_col_1] > thresh_1 ])
    }
   

    z_col_2<-paste("z_avr_" , part + 1   , sep="")
    z_col_p_2<-paste("z_avr_p_" , part + 1  , sep="")
    prc_2 <- 10 * (part + 1 )
    thresh_2<- sort( z_scores [ , z_col_p_2 ])[ n_indis - n_indis_p ]
    swaps_2<-sum(z_scores[, z_col_2] > thresh_2)
    swaps<-intersect(swaps , z_scores$ID[z_scores[ , z_col_2] > thresh_2 ])
    
    f2p<-paste("Z_scores_", tissue , "_perm_", perm, "_" , part , "_" , part + 1 , "_swapped_indis.pdf" , sep="")
    pdf(f2p)
    
    par(mfcol=c(2,3))
    
    hist(z_scores[, z_col_1] ,breaks = 20 , xlab = "Z" , main = paste( "Z scores hist " , tissue ," ", 
                                                                       prc_1, "% eQTLs \n swaps = " , swaps_1 , sep = "" ) )
    abline(v = thresh_1 , col = "red")
    ###################""
    hist(z_scores[, z_col_2 ] ,breaks = 20 , xlab = "Z" , main = paste( "Z scores hist " , tissue , " ",
                                                                        prc_2 ,  "% eQTLs \n swaps = " , swaps_2 , sep = "" ) )
    abline(v = thresh_2 , col = "blue")
    ##############################
    hist(z_scores[, z_col_p_1] ,breaks = 20 , xlab = "Z" , main = paste( "Z scores hist perm ",  perm , "% ", 
                                                                          "\n threshold = " , round(thresh_1 , 3 ) ,
                                                                         "\n expected swaps = " , n_indis_p  , sep = "" ) )
    abline(v = thresh_1 , col = "red")
    #######################################"
    hist(z_scores[, z_col_p_2 ] ,breaks = 20 , xlab = "Z" , main = paste( "Z scores hist  perm " ,  perm , "% ", 
                                                                            "\n threshold = " , round(thresh_2 , 3 )  , sep = "" ) )
    abline(v = thresh_2 , col = "blue")
    ############################################
    
    hist(swapped_indis [, z_col_p_1 ] , xlab = "Z" , main = "Z scores swapped ")
    hist(swapped_indis[ , z_col_p_2] ,  xlab = "Z" , main = "Z scores swapped " )
    
    dev.off()
    
  }
  
  assign(paste( "swaps_" , tissue , sep="") , sort(swaps) ) 
  roww<-data.frame( Tissue = tissue  ,  Swaps = paste (sort(swaps) , collapse = "," ) ) 
  swaps_tab <-rbind(swaps_tab , roww)
}

f2w<-paste(work.dir , "Swapped_samples_perm_" , perm, ".txt" , sep="")
write.table(file = f2w, swaps_tab , row.names =  FALSE , quote = FALSE)

