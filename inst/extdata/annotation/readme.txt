###reactomepathway
annotation is downloaded from http://www.reactome.org/pages/download-data/

Reactome Pathways Gene Set:
http://www.reactome.org/download/current/ReactomePathways.gmt.zip
Note:This is already rolled up!!So parent pathway has all its offsprings' annotation


The REACTOMEPATHWAY ontology is create using with only human pathways:
http://www.reactome.org/download/current/ReactomePathwaysRelation.txt
http://www.reactome.org/download/current/ReactomePathways.txt


####################################
def<-read.table(file = '/home/xin/Desktop/thesis/data/ReactomePathways.def',sep = '\t',header = F,stringsAsFactors = F,quote = "")
unique(def$V3)
def<-def[def$V3=='Homo sapiens',]
write.table(def,file = '/home/xin/Desktop/thesis/data/ReactomePathways.def.human',quote = F,row.names = F,sep = '\t',col.names = F)


def<-read.table(file = '/home/xin/Desktop/thesis/data/ReactomePathways.def.human',sep = '\t',header = F,stringsAsFactors = F,quote = "")

id2term<-split(def$V2,f = def$V1)
relation<-read.table(file = '/home/xin/Desktop/thesis/data/ReactomePathwaysRelation.human',sep = '\t',header = F,stringsAsFactors = F,quote = "")
parent2children<-split(relation$V2,relation$V1)
children2parent<-split(relation$V1,relation$V2)

sink('/home/xin/Workspace/DisEnt/disent/data/reactomePathway.obo')
for(id in names(id2term)){
  cat('[Term]\n')
  cat('id: ',id,"\n",sep='')
  cat('name: ',id2term[[id]],"\n",sep='')
  cat('def: ',id2term[[id]],"\n",sep='')
  if(length(children2parent[[id]])>0){
    for (x in children2parent[[id]]) {
      cat('is_a: ',x," ! ", id2term[[x]],"\n",sep="")
    }
  }
  cat('\n')
}
sink()

annotation<-read.csv('/home/xin/Desktop/thesis/data/ReactomePathways.gmt',header = F,sep = '\t',quote = "",stringsAsFactors = F)
p2g<-split(annotation$V2,annotation$V1)
tmp=unlist(id2term)[match(names(p2g),unlist(id2term))]
p2g<-p2g[which(!is.na(tmp))]
names(p2g)<-names(tmp[which(!is.na(tmp))])
head(p2g,1)
p2g<-sapply(p2g,strsplit,split = ',')
g2p<-revmap(p2g)
names(g2p)[!names(g2p) %in% entrez2symbol()]

require(org.Hs.eg.db)
ls('package:org.Hs.eg.db')
library(AnnotationFuncs)

mapping<-translate(names(g2p), org.Hs.egSYMBOL2EG)
tmp<-list()
for(i in names(g2p)){
  if(length(mapping[[i]])>0){
    for(j in mapping[[i]]){
      tmp[[j]]<-g2p[[i]]
    }
  }
}
p2g.entrez<-revmap(tmp)

sink('/home/xin/Workspace/DisEnt/disent/DisEntR/topOnto/inst/extdata/annotation/reactomepathway')
for(i in names(p2g.entrez)){
  cat(i,'\t',paste(p2g.entrez[[i]],collapse = ','),sep='')
  cat('\n')
}
sink()

sink('/home/xin/Workspace/DisEnt/disent/DisEntR/topOnto/inst/extdata/annotation/reactomepathway_g2p')
for(i in names(tmp)){
  cat(i,'\t',paste(tmp[[i]],collapse = ','),sep='')
  cat('\n')
}
sink()




################################################################################





####GO
annotation is from 



