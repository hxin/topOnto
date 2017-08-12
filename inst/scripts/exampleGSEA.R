setwd('R/')
source.files<-list.files()
sapply(source.files[c(1,3,11,5,6,7,8,9,10)],source)
require(Rgraphviz)


require(topOnto)

##init ontology
topOnto::initONT('HDO')

##read in the annotation
geneID2TERM <- readMappings(file = system.file("extdata/annotation","human_gene2HDO_o", package ="topOnto"))
##or with score
geneID2TERM <- readRDS(system.file("extdata/annotation","human_gene2HDO_weighted_g2t.Rds", package ="topOnto"))
geneID2TERM <- revmap(read.gmt(file = '/Users/xinhe/Documents/Rworkspace/topOnto/inst/extdata/ontology/GSEA/GeneSetDatabases/HDO_o.gmt'))
##clip?
terms<-c("DOID:0050686","DOID:0060083","DOID:1036","DOID:1037","DOID:11868","DOID:1240","DOID:12603","DOID:12965","DOID:14566","DOID:162","DOID:2531","DOID:3264","DOID:4","DOID:6004","DOID:7757","DOID:8692","DOID:8761","DOID:9254","all")
term2geneID<-filter.ontology.annotation(terms,term2genes=revmapWithScore(geneID2TERM))
geneID2TERM<-revmapWithScore(term2geneID)
##


##read exp data
exp<-read.gct('/Users/xinhe/Documents/Rworkspace/topOnto/inst/extdata/ontology/GSEA/Datasets/Leukemia.gct')
pty<-read.cls('/Users/xinhe/Documents/Rworkspace/topOnto/inst/extdata/ontology/GSEA/Datasets/Leukemia.cls')
class.labels<-pty$class.v


##in this example we map the annotation from entrez_id to symbol
e2s<-ograph::entrez2symbol()
names(geneID2TERM)<-e2s[names(geneID2TERM)]
geneID2TERM<-geneID2TERM[!is.na(names(geneID2TERM))]
geneID2TERM<-geneID2TERM[1:50]


ONTdata <- new("topONTdata",ontology = "HDO", annot = annFUN.gene2GO.Score, gene2GO = geneID2TERM,useScore=T,exp=exp,pty=pty)


##perpare premutation
O<-GSEA.GeneRanking(A=as.matrix(exp),class.labels = pty$class.v,nperm = 10,reshuffling.type = 'sample.labels',random.seed= 3338)

##exp type can be c(0,1,2)
##0 equal step
##1 step base on correl
##2 step base on correl and annotation weight
resultClassic <- runTest(ONTdata, algorithm = "classicgsea", statistic = "ks.csw",geneRanking=O,exp.type=1,output.directory='/home/xin/Desktop/tmp/gseaelim/classic/',topgs=10,cutOff=0.05,min.size=1,max.size=500,doc.string='classic')

##elim.type can be c('score','simple')
##score means the annotation weight will be remove from parents once found significant 
##simple means the gene will be remove from parents once found significant 
##elim.gene.type can be c('all','core') means either removing all the genes from parents when found significant, or only remove the core genes
system.time({resultElim <- runTest(ONTdata, algorithm = "elimgsea", statistic = "ks.csw",geneRanking=O,exp.type=2,elim.type='score',
                      elim.gene.type='core',output.directory='~/Desktop/tmp/gseaelim/elim/',topgs=10,cutOff=0.05,min.size=15,max.size=500,doc.string='elimscore')})

resultElim <- runTest(ONTdata, algorithm = "elimgsea", statistic = "ks.csw",geneRanking=O,exp.type=2,elim.type='simple',
                      elim.gene.type='core',output.directory='~/Desktop/tmp/gseaelim/',topgs=10,cutOff=0.99,min.size=2,max.size=2000,doc.string='elimsimple')


##see result
merge(resultClassic@global.report$ALL,resultElim@global.report$ALL,by='GS',all=T,suffixes=c('.classic','.elim'))
merge(resultClassic@global.report$AML,resultElim@global.report$AML,by='GS',all=T,suffixes=c('.classic','.elim'))

###see detail report for a gene set
resultClassic@gs.report$`DOID:12603`

###see detail plot
openPDF(normalizePath(resultClassic@plots$`DOID:1441`))



system.time(O<-GSEA.GeneRanking(A=as.matrix(exp),class.labels = pty$class.v,nperm = 10,reshuffling.type = 'sample.labels',random.seed= 3338))
system.time(resultClassic <- runTest(ONTdata, algorithm = "classicgsea", statistic = "ks.csw",geneRanking=O,exp.type=1,
                                     output.directory='/home/xin/Desktop/tmp/gseaelim/',topgs=10,
                                     cutOff=0.5,min.size=15,max.size=500,doc.string='classic'))




