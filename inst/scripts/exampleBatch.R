#source("http://bioconductor.org/biocLite.R")
cat("##############################################################################\n")
cat("# The script is an example run with HDO & HPO (takes about 10s to run).\n")
cat("# The annotation data is from omim with a test list of 105 genes.\n")
cat("##############################################################################\n")

library(topOnto) 

res<-run.batch(ontologys = c('HDO','HPO'),
          annotation.file = c(system.file("extdata/annotation","human_gene2HDO_o", package ="topOnto"),system.file("extdata/annotation","human_gene2HPO_o", package ="topOnto")),
          gene.file = system.file("extdata/genelist","age", package ="topOnto"),
          algorithm = c('classic','elim','classic'),statistic = c('fisher','fisher','ks'),
          topNodes = 15,useLevels=TRUE,cutoff=1#,clip = NULL
          )

for(i in names(res)){
  print(res[[i]]$tableView)
}

showSigOfNodes(res$HDO$ONTdata, score(res$HDO$result$classicfisher), firstSigNodes = 3, useInfo = 'def')
#printGraph(res$HDO$ONTdata, res$HDO$result$classicfisher, firstSigNodes = 3, res$HDO$result$elimfisher, fn.prefix = "tGO", useInfo = "def")

# def=Term(ONTTERM)
# x=score(resultElimFis)
# y=-log(x)
# freq=y/sum(y)
# wordcloud(words=def[names(y)],freq=freq,scale=c(3,0.1),random.order=FALSE, max.words=30,rot.per=0.35, use.r.layout=FALSE, colors=rev(heat.colors(30, alpha = 1)))
# #brewer.pal(8, 'Dark2')
# termStat(GOdata, names(score(resultElimFis)))[1:10,]
# showSigOfNodes(GOdata, score(resultElimFis), firstSigNodes = 99, useInfo = 'all')
# 
