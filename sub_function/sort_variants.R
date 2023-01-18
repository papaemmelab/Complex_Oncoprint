sort_variants  <- function(vars, genes.order= NULL, group.label=NULL, variants.class= NULL){
  
  vars$nrow <- NULL
  # vars.freq <- ddply(vars, c("GENE"), "nrow",.drop = TRUE)
  # vars <- merge(vars, vars.freq, by="GENE")
  
  gene.list <- NULL
  
  vars <- vars %>% group_by(GENE) %>% mutate(nrow=n())

  if (!is.null(genes.order)){
    gene.list$GENE <- genes.order
    
  } else {
    vars <- vars[order(-vars$nrow),]
    }
  
  new.genes <- data.frame(GENE= c(genes.order, setdiff(vars$GENE, genes.order)) , LAB= group.label) %>% mutate(CLASS= variants.class)

  # gene.list <- rbind(gene.list, new.genes)

  return(list(vars= vars, gene.list= new.genes))
    
  }
