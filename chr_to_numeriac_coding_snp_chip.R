g_data <- read.table("data/chr_to_numeriac_coding_snp_chip_2.txt", header = T)
`%!in%` = Negate(`%in%`)

get_counts <- function(geno_data){
  geno_mat <- as.matrix(geno_data)[, 2:ncol(geno_data)]
  rownames(geno_mat) <- geno_data$Marker
  
  out_counts <- NULL
  for(i in 1:nrow(geno_mat)){
    counts <- table(geno_mat[i, ])
    calls <- names(counts)
    het_calls <- c("R", "Y", "S", "W", "K", "M")
    if(sum(calls %in% het_calls) > 0){
      print(paste0(rownames(geno_mat)[i], " marker has het calls"))
      het <- calls[which(calls %in% het_calls)]
    } else { het <- NA }
    
    counts_no_het <- counts[which(calls %!in% het_calls)]
    calls_no_het <- names(counts_no_het)
    
    if(length(counts_no_het) == 2){
      major_al <- calls_no_het[which(counts_no_het/ncol(geno_mat) > 0.5)]
      minor_al <- calls_no_het[which(counts_no_het/ncol(geno_mat) < 0.5)]
    } else if (length(counts) == 1){
      print(paste0(rownames(geno_mat)[i], " marker is monomorphic"))
      major_al <- calls_no_het[which(counts_no_het/ncol(geno_mat) > 0.5)]
      minor_al <- NA
    }
    
    out_counts <- rbind(out_counts, cbind("Marker" = rownames(geno_mat)[i], major_al, minor_al, het))
  }
  out_counts_df <- as.data.frame(out_counts)
  return(out_counts_df)
}

recode_markers<-function(genodata, n_ind = 19){
  
  count_data <- get_counts(genodata)
  
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

dummy <- recode_markers(g_data)
