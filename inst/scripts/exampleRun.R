#source("http://bioconductor.org/biocLite.R")
cat("##############################################################################\n")
cat("# The script is an example run with HDO(takes about 10s to run).\n")
cat("# The annotation data is from omim with a test list of 105 genes.\n")
cat("##############################################################################\n")

library(topOnto) 

cat("Loading HDO objects from db...\n")
topOnto::initONT('HDO')

a<-system.file("extdata/annotation","human_gene2HDO_o", package ="topOnto")
g<-system.file("extdata/genelist","age", package ="topOnto")

geneID2TERM <- readMappings(file = a)
geneNames=names(geneID2TERM)
myInterestingGenes=(read.csv(header = FALSE, file = g))$V1
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames

GOdata <- new("topGOdata", ontology = "HDO", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2TERM)
resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultElimFis<- runTest(GOdata, algorithm = "elim", statistic = "fisher")
allRes <- GenTable(GOdata,  elimfisher = resultElimFis,fisher = resultFis,topNodes = 15,useLevels=TRUE,cutoff=0.05)
print(allRes)
cat("Demo done..Seems working!\n")


# def=Term(ONTTERM)
# x=score(resultElimFis)
# y=-log(x)
# freq=y/sum(y)
# wordcloud(words=def[names(y)],freq=freq,scale=c(3,0.1),random.order=FALSE, max.words=30,rot.per=0.35, use.r.layout=FALSE, colors=rev(heat.colors(30, alpha = 1)))
# #brewer.pal(8, 'Dark2')
# termStat(GOdata, names(score(resultElimFis)))[1:10,]
# showSigOfNodes(GOdata, score(resultElimFis), firstSigNodes = 99, useInfo = 'all')
# 
