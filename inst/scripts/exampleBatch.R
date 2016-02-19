#source("http://bioconductor.org/biocLite.R")
cat("##############################################################################\n")
cat("# The script is an example run with HDO & HPO (takes about 10s to run).\n")
cat("# The annotation data is from omim with a test list of 105 genes.\n")
cat("##############################################################################\n")

library(topOnto) 
initWHAT()
res<-run.batch(ontologys = c('HDO','HPO','RECTOMEPATHWAY','GOBP','GOCC','GOMF'),
          annotation.file = c(system.file("extdata/annotation","human_gene2HDO", package ="topOnto"),
                              system.file("extdata/annotation","human_gene2HPO", package ="topOnto"),
                              system.file("extdata/annotation","human_gene2reactome", package ="topOnto"),
                              system.file("extdata/annotation","human_gene2GO", package ="topOnto"),
                              system.file("extdata/annotation","human_gene2GO", package ="topOnto"),
                              system.file("extdata/annotation","human_gene2GO", package ="topOnto")
                              ),
          gene.file = system.file("extdata/genelist","ARC", package ="topOnto"),
          algorithm = c('classic','elim','weight','weight01','parentchild'),statistic = c('fisher','fisher','fisher','fisher','fisher'),
          topNodes = 9999,useLevels=TRUE,cutoff=1,clip = NULL,orderBy=2,ranksOf=1
)

for(i in names(res)){
  print(res[[i]]$tableView)
  printGraph(res[[i]]$ONTdata, firstSigNodes = 10, res[[i]]$result$elimfisher,res[[i]]$result$classicfisher,fn.prefix=i,useInfo = "def")
}
# 
# head(res[[3]]$tableView,50)
# 
# showSigOfNodes(res$HDO$ONTdata, score(res$HDO$result$classicfisher), firstSigNodes = 3, useInfo = 'def',)
# showSigOfNodes(res$RECTOMEPATHWAY$ONTdata, score(res$RECTOMEPATHWAY$result$classicfisher), firstSigNodes = 15, useInfo = 'def')
# printGraph(res$HDO$ONTdata, res$HDO$result$classicfisher, firstSigNodes = 3, res$HDO$result$elimfisher, fn.prefix = "tGO", useInfo = "def")
# printGraph(res$RECTOMEPATHWAY$ONTdata, res$RECTOMEPATHWAY$result$classicfisher, firstSigNodes = 5, res$RECTOMEPATHWAY$result$elimfisher, fn.prefix = "tGO", useInfo = "def")

# def=Term(ONTTERM)
# x=score(resultElimFis)
# y=-log(x)
# freq=y/sum(y)
# wordcloud(words=def[names(y)],freq=freq,scale=c(3,0.1),random.order=FALSE, max.words=30,rot.per=0.35, use.r.layout=FALSE, colors=rev(heat.colors(30, alpha = 1)))
# #brewer.pal(8, 'Dark2')
# termStat(GOdata, names(score(resultElimFis)))[1:10,]
# showSigOfNodes(GOdata, score(resultElimFis), firstSigNodes = 99, useInfo = 'all')
# 

cat("Loading HDO objects from db...\n")
topOnto::initONT('HDO')

a<-'/home/xin/Workspace/DisEnt/disent/DisEntR/topOnto/inst/extdata/annotation/human_gene2GO'
g<-system.file("extdata/genelist","ARC", package ="topOnto")

geneID2TERM <- readMappings(file = a)
geneNames=names(geneID2TERM)
myInterestingGenes=(read.csv(header = FALSE, file = g))$V1
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames

##clip?
# terms<-c('DOID:10652')
# term2geneID<-filter.ontology.annotation(terms,term2genes=revmap(geneID2TERM))
# geneID2TERM<-revmap(term2geneID)

ONTdata <- new("topONTdata", ontology = "GOCC", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2TERM)
resultFis <- runTest(ONTdata, algorithm = "classic", statistic = "fisher")
resultElimFis<- runTest(ONTdata, algorithm = "elim", statistic = "fisher")
allRes <- GenTable(ONTdata,  fisher = resultFis,elimfisher = resultElimFis,topNodes = 50,useLevels=TRUE,cutoff=1,orderBy=2,ranksOf=1)
GenTable(ONTdata,  fisher = resultFis,elimfisher = resultElimFis,topNodes = 50,useLevels=TRUE,cutoff=1,orderBy=1,ranksOf=1)
print(allRes)
cat("Demo done..Seems working!\n")

