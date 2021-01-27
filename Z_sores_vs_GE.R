setwd("/media/julia/e4a91c72-ba4c-44e5-a062-4bba142849b4/CEDAR_GE_GIGA/")

ge<-read.table("CD15/GE_CD15.txt" , header = TRUE , stringsAsFactors = FALSE)
cis_eqtls<-read.table ("CD15/Cis_eQTLs_sign_per_gene_CD15.txt" , header = TRUE , stringsAsFactors =  FALSE)

nrow(ge)
ge<-ge[ge$id %in% cis_eqtls$gene ,  ]
nrow(ge)
mu<-rowMeans(ge [, 2 : ncol(ge)] , na.rm = TRUE)
sigma<-sqrt(apply(ge[, 2 : ncol(ge) ] , 1 , var ) )

z_score<-function(s){
  
  return ( abs(s - mu)/sigma )
}  

d<-data.frame()

# Can you subset ge data for cis_eqtl only, thank you 

for ( i in 2 : ncol(ge) ){
  indi <- colnames(ge)[i]
  gee <- ge[ , i]
  z_scores<- z_score(gee)
  z_scores <- sum(z_scores , na.rm =  TRUE) / sum(!is.na(z_scores))
  roww<-data.frame(ID  = indi , z_score = z_scores )
  if (nrow(d) == 0 ){
    d<-roww
  }else{
    d<-rbind(d , roww)
  }
}

z_scores<-read.table( "CD15/Z_scores_CD15.txt" , header = TRUE , stringsAsFactors = FALSE)


par( mfrow = c(1, 2))
hist(z_scores$z_avr, breaks = 30 , main = "Z-scores vs GeMu ")
abline(v = mean(z_scores$z_avr , na.rm =  TRUE ), col = "red")
hist(d$z_score , breaks = 30 , main = "Z-score vs Mu ")
abline(v = mean(d$z_score , na.rm =  TRUE)  , col = "red")

par( mfrow = c(1 , 1 ) )
pdf("CD15/Z_Score_vs_mu.pdf")
hist(d$z_score , breaks = 30 , main = "Z-score vs Mu ")
abline(v = mean(d$z_score , na.rm =  TRUE)  , col = "red")
dev.off()
