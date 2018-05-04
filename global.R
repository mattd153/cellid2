#30/4/18
## heatmap functions
list.of.packages <- c("tidyverse", "reshape")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(ggplot2)
library(tidyverse)
library(reshape)


# Compute the values for the heatmap
computeMarkersOverlap <- function(sample_1, sample_2, denom) {
  sample_1 <- as.data.frame(sample_1)
  sample_2 <- as.data.frame(sample_2)
  sample_1 <- split(sample_1, f = sample_1$cluster)
  sample_2 <- split(sample_2, f = sample_2$cluster)
  sample_1_matrix <- matrix(nrow=length(sample_1), ncol = length(sample_2))
  for (i in 1:length(sample_1)){ #sample_1 -> the rows
    for (j in (1:length(sample_2))){ #sample_2 -> the columns
      ##fill the table using sample 1 as the denominator
      #if table(s1 %in% s2) has 2 dimensions, then there will be some shared genes and some not shared
      if (dim(table((str_trim(unique(sample_1[[i]]$external_gene_name))) %in% (str_trim(unique(sample_2[[j]]$external_gene_name)))))==2){
        #if there are some true and some false, there will be two columns
        shared = table(unique(sample_1[[i]]$external_gene_name) %in% unique(sample_2[[j]]$external_gene_name))[[2]]
        # User specifies denominator as sample_1 or sample_2
        if (denom=="Sample_1"){
          denominator = length(str_trim(unique((sample_1[[i]]$external_gene_name))))
        }
        else if (denom=="Sample_2"){
          denominator = length(str_trim(unique((sample_2[[j]]$external_gene_name))))
        }
        else
          stop("denom argument must equal \"Sample_1\" or \"Sample 2\"")
        
        #this populates the table
        sample_1_matrix[i,j] <- shared/(denominator) ## this uses sample_1 as the denominator
      }
      #if all of the genes of s1 are in s2
      else if (all(str_trim(unique(sample_1[[i]]$external_gene_name)) %in% str_trim(unique(sample_2[[j]]$external_gene_name)))){
        sample_1_matrix[i,j] <- 1
      }
      #if none of the genes of s1 are in s2
      else 
        sample_1_matrix[i,j] <- (0)
    }}
  sample_1_matrix
  rnms <- names(sample_1)
  cnms <- names(sample_2)
  colnames(sample_1_matrix) <- cnms
  rownames(sample_1_matrix) <- rnms
  
  markers_overlap <- reshape::melt(sample_1_matrix)
  
  # Put the clusters in alphanumeric order.
  markers_overlap$X1 <- factor(markers_overlap$X1, as.character(unique(markers_overlap$X1)))
  markers_overlap$X2 <- factor(markers_overlap$X2, as.character(unique(markers_overlap$X2)))
  return(markers_overlap)
}

computeMarkersOverlap_ensid <- function(sample_1, sample_2, denom) {
  sample_1 <- as.data.frame(sample_1)
  sample_2 <- as.data.frame(sample_2)
  sample_1 <- split(sample_1, f = sample_1$cluster)
  sample_2 <- split(sample_2, f = sample_2$cluster)
  sample_1_matrix <- matrix(nrow=length(sample_1), ncol = length(sample_2))
  for (i in 1:length(sample_1)){ #sample_1 -> the rows
    for (j in (1:length(sample_2))){ #sample_2 -> the columns
      #if table(s1 %in% s2) has 2 dimensions, then there will be some shared genes and some not shared
      if (dim(table(str_trim(unique(sample_1[[i]]$ensembl_gene_id)) %in% str_trim(unique(sample_2[[j]]$ensembl_gene_id))))==2){
        #if there are some true and some false, there will be two columns, this take the number in the TRUE column
        shared = table(str_trim(unique(sample_1[[i]]$ensembl_gene_id)) %in% str_trim(unique(sample_2[[j]]$ensembl_gene_id)))[[2]]
        # User specifies denominator as sample_1 or sample_2
        if (denom=="Sample_1"){
          denominator = length(str_trim(unique(sample_1[[i]]$ensembl_gene_id)))
        }
        else if (denom=="Sample_2"){
          denominator = length(str_trim(unique(sample_2[[j]]$ensembl_gene_id)))
        }
        else
          stop("denom argument must equal \"Sample_1\" or \"Sample 2\"")
        sample_1_matrix[i,j] <- shared/(denominator) 
      }
      #if all of the genes of s1 are in s2
      else if (all(str_trim(unique(sample_1[[i]]$ensembl_gene_id)) %in% str_trim(unique(sample_2[[j]]$ensembl_gene_id)))){
        sample_1_matrix[i,j] <- 1
      }
      #if none of the genes of s1 are in s2
      else 
        sample_1_matrix[i,j] <- (0)
    }}
  sample_1_matrix
  rnms <- names(sample_1)
  cnms <- names(sample_2)
  colnames(sample_1_matrix) <- cnms
  rownames(sample_1_matrix) <- rnms
  
  markers_overlap <- reshape::melt(sample_1_matrix)
  
  # Put the clusters in alphanumeric order.
  markers_overlap$X1 <- factor(markers_overlap$X1, as.character(unique(markers_overlap$X1)))
  markers_overlap$X2 <- factor(markers_overlap$X2, as.character(unique(markers_overlap$X2)))
  return(markers_overlap)
}


# Create the heatmap
heatmapMarkersOverlap <- function(markers_overlap, sample_1, sample_2) {
  title= paste0("Portion of marker genes shared in each cluster", "\n")
  #ylab=paste0(((substitute(sample_1))), " cluster" )
  #xlab=paste0(((substitute(sample_2)))," cluster") 
  p <- ggplot(markers_overlap, aes(X2, X1)) + #X1 is sample_1, X2 is sample_2
    geom_tile(aes(fill = value)) + 
    geom_text(aes(label = round(value, 2))) +
    scale_fill_gradient(low = "white", high = "red", limits = c(0,1)) +
    labs(x=xlab, y =ylab)  +
    ggtitle(title) +
    theme(plot.margin=unit(c(0.5,1,1.5,1.2),"cm"), 
          plot.title = element_text(hjust = 0.5, size = 25), ###
          axis.text.x=element_text(angle=90, hjust=1, size = 12),
          axis.text.y = element_text(hjust=1,size=12, vjust = 0.5),
          panel.background = element_blank(), 
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          axis.line = element_line(colour = "black")) ###
  return(p)
}


## Function that returns a list of genes
marker_return <- function(sample_1, sample_2, cluster_s1, cluster_s2, genelist_type){
  sample_1 <- as.data.frame(sample_1)
  sample_2 <- as.data.frame(sample_2)
  #get sections of tables
  gtable_sc1 <- sample_1[sample_1$cluster==cluster_s1,]
  gtable_sc2 <- sample_2[sample_2$cluster==cluster_s2,]
  #index for intersects
  if (genelist_type=="Gene_Symbols"){
    common_genes <- intersect(str_trim(unique(gtable_sc1$external_gene_name)), str_trim(unique(gtable_sc2$external_gene_name)))
    other_genes_cs1 <- setdiff(str_trim(unique(gtable_sc1$external_gene_name)), str_trim(unique(gtable_sc2$external_gene_name)))
    other_genes_cs2 <- setdiff(str_trim(unique(gtable_sc2$external_gene_name)), str_trim(unique(gtable_sc1$external_gene_name)))
  }
  else if (genelist_type=="Ensemble_IDs"){
    common_genes <- intersect(str_trim(unique(gtable_sc1$ensembl_gene_id)), str_trim(unique(gtable_sc2$ensembl_gene_id)))
    other_genes_cs1 <- setdiff(str_trim(unique(gtable_sc1$ensembl_gene_id)), str_trim(unique(gtable_sc2$ensembl_gene_id)))
    other_genes_cs2 <- setdiff(str_trim(unique(gtable_sc2$ensembl_gene_id)), str_trim(unique(gtable_sc1$ensembl_gene_id)))
  }
  shared <- paste0(common_genes, collapse = ", ")
  not_shared_cs1 <- paste0(other_genes_cs1, collapse = ", ")          
  not_shared_cs2 <- paste0(other_genes_cs1, collapse = ", ")
  
  out1 <- paste0("Shared: ", shared)
  out2 <- paste0("Other genes in Cluster 1: ", not_shared_cs1)
  out3 <- paste0("Other genes in Cluster 2: ", not_shared_cs2)
  all_out <- list(out1, out2, out3)
  return(all_out)
}



# Barplot functions

### Now the individual cluster barplots
cluster_barplot <- function(sample_1,sample_2, denom, clust_select){
  ## Clust select is an individual cluster from the cluster sepcified
  ## by "denom =" 
  sample_1 <- as.data.frame(sample_1)
  sample_2 <- as.data.frame(sample_2)
  sample_1 <- split(sample_1, f = sample_1$cluster)
  sample_2 <- split(sample_2, f = sample_2$cluster)
  
  # The three columns of bartable_s2
  ct_vector1 <- c() ## The name of each cell type, repeated twice
  sns_vector1 <- c() ## The words "Shared" and "out of"
  ng_vector1 <- c() ## The number of genes (shared and not shared for each cell type)
  # And bartable_s2
  ct_vector2 <- c()
  sns_vector2 <- c()
  ng_vector2 <- c()
  
  #create the blank data frames
  bar_table_s1 <- matrix(ncol = 3, nrow = (length(names(sample_1)))*2)
  colnames(bar_table_s1) <- c("cell_type","numerator/denominator", "number_of_genes")
  bar_table_s1 <- as.data.frame(bar_table_s1)
  bar_table_s2 <- matrix(ncol = 3, nrow = (length(names(sample_2)))*2)
  colnames(bar_table_s2) <- c("cell_type","numerator/denominator", "number_of_genes")
  bar_table_s2 <- as.data.frame(bar_table_s2)
  
  celltype_s1 <- names(sample_1)
  celltype_s2 <- names(sample_2)
  
  # Fill the data frames
  # sample 1 data
  for (i in 1:length(sample_1)){ #sample_1
    #shared, out of
    sns1 <- c("shared","not_shared")
    sns_vector1 <- c(sns_vector1, sns1)
    #celltype
    ct1 <- rep(celltype_s1[i], 2)
    ct_vector1 <- append(ct_vector1, ct1)
    #number of genes
    if (dim(table(str_trim(unique(sample_1[[i]]$external_gene_name)) %in% str_trim(unique(sample_2[[clust_select]]$external_gene_name))))==2){
      shared1 = table(str_trim(unique(sample_1[[i]]$external_gene_name)) %in% str_trim(unique(sample_2[[clust_select]]$external_gene_name)))[[2]]
      out_of1 <- length(str_trim(unique(sample_2[[clust_select]]$external_gene_name)))-shared1
    }
    else if (all((sample_1[[i]]$external_gene_name) %in% (sample_2[[clust_select]]$external_gene_name))){
      shared1 = length(str_trim(unique(sample_2[[clust_select]]$external_gene_name)))
      out_of1 <- 0
    }
    else {
      shared1 = length(str_trim(unique(sample_2[[clust_select]]$external_gene_name)))
      out_of1 <- 0
    }
    ng1 <- c(shared1, out_of1)
    ng_vector1 <- c(ng_vector1, ng1)  
  }
  bar_table_s1$cell_type <- ct_vector1
  bar_table_s1$`numerator/denominator` <- sns_vector1
  bar_table_s1$number_of_genes <- ng_vector1
  
  
  #sample_2 table
  for (j in (1:length(sample_2))){ #sample_2 
    #shared, out of
    sns2 <- c("shared","not_shared")
    sns_vector2 <- c(sns_vector2, sns2)
    #celltype
    ct2 <- rep(celltype_s2[j], 2)
    ct_vector2 <- append(ct_vector2, ct2)
    #number of genes
    if (dim(table((unique(sample_2[[j]]$external_gene_name)) %in% (unique(sample_1[[clust_select]]$external_gene_name))))==2){
      shared2 <- table((unique(sample_2[[j]]$external_gene_name)) %in% (unique(sample_1[[clust_select]]$external_gene_name)))[[2]] # returns number of TRUE cases
      out_of2 <- length(str_trim(unique(sample_1[[clust_select]]$external_gene_name))) -shared2
    } 
    else if (all((sample_1[[clust_select]]$external_gene_name) %in% (sample_2[[j]]$external_gene_name))){
      shared2 <- length(str_trim(unique(sample_1[[clust_select]]$external_gene_name)))
      out_of2 <- 0
    }
    else { 
      shared2 <- 0
      out_of2 <- length(str_trim(unique(sample_1[[clust_select]]$external_gene_name)))
    }
    ng2 <- c(shared2, out_of2)
    ng_vector2 <- c(ng_vector2, ng2)
  }
  bar_table_s2$cell_type <- ct_vector2
  bar_table_s2$`numerator/denominator` <- sns_vector2
  bar_table_s2$number_of_genes <- ng_vector2
  
  if (denom=="Sample_2"){
    title <- paste0("Cluster: ", clust_select, " of ", denom)
    p <- ggplot(data = bar_table_s1, aes(x=cell_type,y=number_of_genes, fill=`numerator/denominator`, label = number_of_genes))+ geom_bar(stat="identity") 
    final_p <- p + coord_flip() + guides(fill = guide_legend(reverse = TRUE)) + ggtitle(title)+ 
      theme(plot.title = element_text(hjust = 0.5, size = 25),
            axis.text.y = element_text(hjust=1,size=12, vjust = 0.5),
            axis.title.x = element_text(size=20),
            axis.title.y = element_text(size=20),
            panel.background = element_blank()
            )+geom_text(size = 4)
    show(final_p)
  }
  
  if (denom=="Sample_1"){
    title <- paste0("Cluster: ", clust_select, " of ", denom)
    p <- ggplot(data = bar_table_s2, aes(x=cell_type,y=number_of_genes, fill=`numerator/denominator`, label = number_of_genes))+ geom_bar(stat="identity") 
    final_p <- p + coord_flip() + guides(fill = guide_legend(reverse = TRUE)) + ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, size = 25),
          axis.text.y = element_text(hjust=1,size=12, vjust = 0.5),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          panel.background = element_blank()
          )+geom_text(size = 4)
    show(final_p)    
  }
}






