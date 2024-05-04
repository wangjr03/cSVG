fn_cluster_genes<-function(SE_genes,list1,DEC_genes){

  final1=cbind(SE_genes,list1)
  final=final1[!final1$SE_genes %in% DEC_genes,]
  items=strsplit(final$V1,",")
  unique_genes=c()
  mat_new=matrix(0,nrow=1,ncol=2)
  for(g in 1:nrow(final)){
    gene_list=c(final[g,2],unlist(items[g]))
    gene_list=SE_genes[as.numeric(gene_list)]
    if(length(gene_list)==1){
      unique_genes=c(unique_genes,gene_list)
    }else{
      mat=t(combn(gene_list, 2))
      mat_new=rbind(mat_new,mat)
    }
  }
  el=mat_new[-1,]
  library(igraph)
  gr=graph_from_edgelist(el, directed = FALSE)
  com=cluster_leiden(gr,
                   objective_function = "modularity",
                   weights = NULL,
                   resolution_parameter = 1,
                   beta = 0.01,
                   initial_membership = NULL,
                   n_iterations = 5,
                   vertex_weights = NULL
  )
  out1=membership(com) 
  other_genes=DEC_genes[!DEC_genes %in% names(out1)]
  if(length(other_genes)>1){
    clusters=(max(out1)+1):(max(out1)+length(other_genes))
    names(clusters)=other_genes
    out1=c(out1,clusters)
  }

return(out1)
}