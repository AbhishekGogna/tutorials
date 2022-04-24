c_data <- read.table("chr_to_numeriac_coding_snp_chip_1.txt", header = T)
g_data <- read.table("chr_to_numeriac_coding_snp_chip_2.txt", header = T)

recode_markers<-function(count_data, genodata, n_ind = 19){
  output<-matrix(0 ,dim(count_data)[1], n_ind+1)
  colnames(genodata)[1] <- colnames(count_data)[1]<- "Marker"
  colnames(output)<-colnames(genodata)
  
  system.time({for (i in 1:dim(count_data)[1]){
    marker<-as.matrix(genodata[genodata$Marker==count_data$Marker[i],])
    marker[marker==count_data[i, "major_al"]] <- 0
    marker[marker==count_data[i, "minor_al"]] <- 2
    marker[marker==count_data[i, "het"]] <- 1
    marker[marker=="failed"] <- NA
    output[i,]<-marker
    if (i%%5 == 0){print(paste0("recoding complete for ", i, " marker of ", dim(count_data)[1] ))}
  }})
  
  GNdata<-apply(output[,-1], 2, function(x) as.numeric(x))
  rownames(GNdata)<-output[,1]
  
  return(GNdata)
}

dummy <- recode_markers(c_data, g_data)
