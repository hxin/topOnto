library('ograph')
ograph::initOGraph('HDO')
g<-new("ontGraph",ontology='HDO')
mg<-subGraphByNodes(g@graph,nodes=c('DOID:10652','DOID:14330','DOID:5419'))
mg<-mapGene2Graph(mg,file=system.file("extdata/annotation","human_gene2HDO", package ="topOnto"),rollup=TRUE)
GO2geneID=.get.genesInNodes(mg,V(mg)$name)
e2s<-entrez2symbol()
GO2geneID.symbol<-lapply(GO2geneID,function(x){
  unname(e2s[x][!is.na(e2s[x])])
})
terms<-Term(ONTTERM)

result<-elim(mg,GO2geneID.symbol)

# temp <- readLines('/home/xin/Downloads/GSEA-P-R/GeneSetDatabases/C2.gmt')
# 
# l=list()
# l<-lapply(temp,function(row){
#   row=unlist(strsplit(row, "\t"))
#   #name=row[1]
#   row<-row[-1]
#   row<-row[-1]
#   row
# })
# 
# names(l)<-sapply(temp,function(row){
#   row=unlist(strsplit(row, "\t"))
#   row[1]
# })



##update gene set
for(x in names(GO2geneID.symbol)){
  if(exists(x,result$elimGenes.LookUP, mode = "character"))
    GO2geneID.symbol[[x]]=setdiff(GO2geneID.symbol[[x]],get(x,result$elimGenes.LookUP, mode = 'character'))
}

gs.db=NULL
terms<-Term(ONTTERM)
for(i in names(GO2geneID.symbol)){
  s=paste(i,terms[[i]],paste(GO2geneID.symbol[[i]],collapse="\t"),sep="\t")
  gs.db=c(gs.db,s)
  s=''
}

r=GSEA(                                                                    # Input/Output Files :-------------------------------------------
                                                                          input.ds =  "/home/xin/Downloads/GSEA-P-R/Datasets/Lung_Stanford.gct",           # Input gene expression Affy dataset file in RES or GCT format
                                                                          input.cls = "/home/xin/Downloads/GSEA-P-R/Datasets/Lung_Stanford.cls",           # Input class vector (phenotype) file in CLS format
                                                                          #gs.db =     "/home/xin/Downloads/GSEA-P-R/GeneSetDatabases/C2.gmt",         # Gene set database in GMT format
                                                                          gs.db =     gs.db,         # Gene set database in GMT format
                                                                          output.directory      = "/home/xin/Downloads/GSEA-P-R/TEST/",        # Directory where to store output and results (default: "")
                                                                          #  Program parameters :-------------------------------------------------------------------------------------------------------------------------
                                                                          doc.string            = "TeSt",   # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
                                                                          non.interactive.run   = F,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
                                                                          reshuffling.type      = "sample.labels", # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels" 
                                                                          nperm                 = 1000,              # Number of random permutations (default: 1000)
                                                                          weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
                                                                          nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
                                                                          fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
                                                                          fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
                                                                          topgs                 = 20,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
                                                                          adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
                                                                          gs.size.threshold.min = 5,              # Minimum size (in genes) for database gene sets to be considered (default: 25)
                                                                          gs.size.threshold.max = 800,             # Maximum size (in genes) for database gene sets to be considered (default: 500)
                                                                          reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
                                                                          preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
                                                                          random.seed           = 3338,            # Random number generator seed. (default: 123456)
                                                                          perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
                                                                          fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
                                                                          replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
                                                                          save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
                                                                          OLD.GSEA              = F,               # Use original (old) version of GSEA (default: F)
                                                                          use.fast.enrichment.routine = T          # Use faster routine to compute enrichment for random permutations (default: T)
)


#plot sig
plotSiginificant(mg,r$report1,file='/home/xin/Desktop/Untitled Folder/R_noElim/sig1')


plotSiginificant<-function(graph,report,number_of_node=10,file='',sub=TRUE){
  ids=as.vector(report$GS)
  ps=as.numeric(as.vector(report$'NOM p-val'))
  names(ps)=ids
  x=sort(ps)
  
  if(length(x)>number_of_node)
    x=x[1:number_of_node]
  
  log.x=log10(x)
  color <- round(log.x - range(log.x)[1] + 1)
  colorMap <- heat.colors(max(color))
  
  color<-sapply(names(color),function(x){
    colorMap[color[x]]
  })
  
  if(sub){
    s<-subGraphByNodes(graph,nodes=ids)
    graph.result=ograph::set.node.attribute(s,attr_name='color',attr_value=color,nodes=names(color))
  }else{
    graph.result=ograph::set.node.attribute(graph,attr_name='color',attr_value=color,nodes=names(color))
  }
  
  if(file==''){
    treeplot(graph.result,label=0,vertex.size=5,vertex.label.cex=1)
  }else{
    plot2file(filename=file)
    treeplot(graph.result,label=0,vertex.size=5,vertex.label.cex=1)
    dev.off()
  }
}