tissue<-"CD15"
eqtls_prc<-10
work.dir<-"/media/julia/e4a91c72-ba4c-44e5-a062-4bba142849b4/CEDAR_GE_GIGA/"

args= commandArgs()
i = 1
while(args[i]!="--args"){
  cat(args[i])
  i=i+1
}
tissue<-args  [ i + 1 ]
work.dir <- args[ i + 2 ]
eqtls_prc <-as.numeric(args[ i + 3 ])
library(plot.matrix)

tissues <-c("CD4" , "CD8" ,"CD14" , "CD15" , "CD19" , "PLA" , "IL" , "TR" , "RE")
total_swapped<- 0
for(tissue in tissues){
  setwd(paste(work.dir , "/" , tissue , sep = "") )
  
  all_vs_all <-read.table( paste("Z_scores_all_vs_all_" , tissue , "_" , eqtls_prc , "_prc_eqtls.txt" , sep = "" ) )
  min_indxs<-apply(all_vs_all , 1, which.min)
  Min_values<-apply(all_vs_all , 1, min)
  
  Diag <-diag(as.matrix(all_vs_all) )
  
  d<-data.frame(Ind = rownames(all_vs_all) , Min_ind = colnames(all_vs_all) [min_indxs] , Diag, Min_values)
  d$Ind<-as.character(d$Ind)
  d$Min_ind <-as.character(d$Min_ind)
  d<-cbind(d , Dist = abs(d$Diag - d$Min_values) )
  indxs<-which(as.character(d$Ind) != as.character(d$Min_ind) )
  total_swapped<- total_swapped + length(indxs)
  print(paste(tissue , " " , length(indxs)))
  bad<-d[indxs, ]
  indis.all<-unique(c(bad$Ind , bad$Min_ind))
  indxs<-which(d$Ind %in% indis.all )
  bad<-d[indxs, ]
  
  f2w<-paste("Swapped_samples_" , tissue , "_v2.txt" , sep="")
  write.table(file  = f2w , bad , row.names = FALSE , quote = FALSE)
  
  indxs<-which(rownames(all_vs_all) %in% indis.all )
  z_scores<-as.matrix(all_vs_all[indxs , indxs])
  ftp<-paste("Z_scores_swapped_", tissue , "_3.pdf" , sep="") 
  pdf(ftp)
  if(tissue != "CD15"){
    plot(z_scores,  main = tissue, 
         fmt.cell='%.2f', text.cell=list(cex=0.5),
         breaks=c(0.5, 1.3),
         axis.col=list(side=1, las=2), axis.row = list(side=2, las=1) , cex.axis = 0.5)
    
  }else{
    plot(z_scores,  main = tissue, 
         border=NA,
         breaks=c(0.5, 1.3),
         axis.col=list(side=1, las=2), axis.row = list(side=2, las=1) , cex.axis = 0.5)
    
  }
  
  dev.off()
  
  
  z_scores<-as.matrix(all_vs_all)
  ftp<-paste("Z_scores_all_", tissue , "_3.pdf" , sep="") 
  pdf(ftp)
  plot(z_scores,  main = tissue, 
       border=NA,
       breaks=c(0.5, 1.3),fmt.key="%.2f", cex = 0.6,
       axis.col=NULL, axis.row=NULL
      )
  dev.off()
  
  print(paste("ready with " , tissue , sep=""))
  
}

print(total_swapped)
work.dir<-"/media/julia/e4a91c72-ba4c-44e5-a062-4bba142849b4/CEDAR_GE_GIGA/"
setwd(work.dir)

reassignment<-read.table("Samples_Reassignment.csv" , header = TRUE , stringsAsFactors = FALSE)
pdf("Samples_Swapped_wrongly_assigned.pdf")
t<-table(reassignment$Tissue)
ttext<-paste("Swaps + Wrongly assigned samples\n" , sum(t) , " total\n", 
             sum(reassignment$Type == "Swap" ) , " swapped", sep = "")
x<-barplot(t, col = "cornflowerblue" , main= ttext , ylim = c(0, 100) )
text(x , y =t + 2 , labels = as.character(t))
dev.off()
