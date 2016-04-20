read.gct<-function(filename=NULL){
  content <- readLines(filename)
  content <- content[-1]
  content <- content[-1]
  col.names <- noquote(unlist(strsplit(content[1], "\t")))
  col.names <- col.names[c(-1, -2)]
  num.cols <- length(col.names)
  content <- content[-1]
  content <- content[which(content!="")]
  num.lines <- length(content)
  
  row.nam <- vector(length=num.lines, mode="character")
  row.des <- vector(length=num.lines, mode="character")
  m <- matrix(0, nrow=num.lines, ncol=num.cols)
  
  for (i in 1:num.lines) {
    line.list <- noquote(unlist(strsplit(content[i], "\t")))
    row.nam[i] <- noquote(line.list[1])
    row.des[i] <- noquote(line.list[2])
    line.list <- line.list[c(-1, -2)]
    for (j in 1:length(line.list)) {
      m[i, j] <- as.numeric(line.list[j])
    }
  }
  ds <- data.frame(m)
  names(ds) <- col.names
  row.names(ds) <- row.nam
  return(ds)
}

read.cls<-function(file=NULL){
  cls.cont <- readLines(file)
  num.lines <- length(cls.cont)
  class.list <- unlist(strsplit(cls.cont[[3]], " "))
  s <- length(class.list)
  t <- table(class.list)
  l <- length(t)
  phen <- vector(length=l, mode="character")
  phen.label <- vector(length=l, mode="numeric")
  class.v <- vector(length=s, mode="numeric")
  for (i in 1:l) {
    phen[i] <- noquote(names(t)[i])
    phen.label[i] <- i - 1
  }
  for (i in 1:s) {
    for (j in 1:l) {
      if (class.list[i] == phen[j]) {
        class.v[i] <- phen.label[j]
      }
    }
  }
  return(list(phen = phen, class.v = class.v))
}

GSEA.GeneRanking <- function(A, class.labels, nperm=10, reshuffling.type='sample.labels', permutation.type = 0, sigma.correction = "GeneCluster", fraction=1.0, replace=F, reverse.sign= F,random.seed=123456) { 
  
  # This function ranks the genes according to the signal to noise ratio for the actual phenotype and also random permutations and bootstrap  
  # subsamples of both the observed and random phenotypes. It uses matrix operations to implement the signal to noise calculation 
  # in stages and achieves fast execution speed. It supports two types of permutations: random (unbalanced) and balanced. 
  # It also supports subsampling and bootstrap by using masking and multiple-count variables.  When "fraction" is set to 1 (default)
  # the there is no subsampling or boostrapping and the matrix of observed signal to noise ratios will have the same value for 
  # all permutations. This is wasteful but allows to support all the multiple options with the same code. Notice that the second 
  # matrix for the null distribution will still have the values for the random permutations 
  # (null distribution). This mode (fraction = 1.0) is the defaults, the recommended one and the one used in the examples.
  # It is also the one that has be tested more thoroughly. The resampling and boostrapping options are intersting to obtain 
  # smooth estimates of the observed distribution but its is left for the expert user who may want to perform some sanity 
  # checks before trusting the code.
  #
  # Inputs:
  #   A: Matrix of gene expression values (rows are genes, columns are samples) 
  #   class.labels: Phenotype of class disticntion of interest. A vector of binary labels having first the 1's and then the 0's 
  #   gene.labels: gene labels. Vector of probe ids or accession numbers for the rows of the expression matrix 
  #   nperm: Number of random permutations/bootstraps to perform 
  #   permutation.type: Permutation type: 0 = unbalanced, 1 = balanced. For experts only (default: 0) 
  #   sigma.correction: Correction to the signal to noise ratio (Default = GeneCluster, a choice to support the way it was handled in a previous package) 
  #   fraction: Subsampling fraction. Set to 1.0 (no resampling). For experts only (default: 1.0) 
  #   replace: Resampling mode (replacement or not replacement). For experts only (default: F) 
  #   reverse.sign: Reverse direction of gene list (default = F)
  #
  # Outputs:
  #   s2n.matrix: Matrix with random permuted or bootstraps signal to noise ratios (rows are genes, columns are permutations or bootstrap subsamplings
  #   obs.s2n.matrix: Matrix with observed signal to noise ratios (rows are genes, columns are boostraps subsamplings. If fraction is set to 1.0 then all the columns have the same values
  #   order.matrix: Matrix with the orderings that will sort the columns of the obs.s2n.matrix in decreasing s2n order
  #   obs.order.matrix: Matrix with the orderings that will sort the columns of the s2n.matrix in decreasing s2n order
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  set.seed(seed=random.seed, kind = NULL)
  A <- A + 0.00000001
  
  N <- nrow(A)
  Ns <- ncol(A)
  
  subset.mask <- matrix(0, nrow=Ns, ncol=nperm)
  reshuffled.class.labels1 <- matrix(0, nrow=Ns, ncol=nperm)
  reshuffled.class.labels2 <- matrix(0, nrow=Ns, ncol=nperm)
  class.labels1 <- matrix(0, nrow=Ns, ncol=nperm)
  class.labels2 <- matrix(0, nrow=Ns, ncol=nperm)
  
  order.matrix <- matrix(0, nrow = N, ncol = nperm)
  obs.order.matrix <- matrix(0, nrow = N, ncol = nperm)
  s2n.matrix <- matrix(0, nrow = N, ncol = nperm)
  obs.s2n.matrix <- matrix(0, nrow = N, ncol = nperm)
  
  obs.gene.labels <- vector(length = N, mode="character")
  obs.gene.descs <- vector(length = N, mode="character")
  obs.gene.symbols <- vector(length = N, mode="character")
  
  M1 <- matrix(0, nrow = N, ncol = nperm)
  M2 <- matrix(0, nrow = N, ncol = nperm)
  S1 <- matrix(0, nrow = N, ncol = nperm)
  S2 <- matrix(0, nrow = N, ncol = nperm)
  
  C <- split(class.labels, class.labels)
  class1.size <- length(C[[1]])
  class2.size <- length(C[[2]])
  class1.index <- seq(1, class1.size, 1)
  class2.index <- seq(class1.size + 1, class1.size + class2.size, 1)
  
  for (r in 1:nperm) {
    class1.subset <- sample(class1.index, size = ceiling(class1.size*fraction), replace = replace)
    class2.subset <- sample(class2.index, size = ceiling(class2.size*fraction), replace = replace)
    class1.subset.size <- length(class1.subset)
    class2.subset.size <- length(class2.subset)
    subset.class1 <- rep(0, class1.size)
    for (i in 1:class1.size) {
      if (is.element(class1.index[i], class1.subset)) {
        subset.class1[i] <- 1
      }
    }
    subset.class2 <- rep(0, class2.size)
    for (i in 1:class2.size) {
      if (is.element(class2.index[i], class2.subset)) {
        subset.class2[i] <- 1
      }
    }
    subset.mask[, r] <- as.numeric(c(subset.class1, subset.class2))
    fraction.class1 <- class1.size/Ns
    fraction.class2 <- class2.size/Ns
    
    if (permutation.type == 0) { # random (unbalanced) permutation
      full.subset <- c(class1.subset, class2.subset)
      label1.subset <- sample(full.subset, size = Ns * fraction.class1)
      reshuffled.class.labels1[, r] <- rep(0, Ns)
      reshuffled.class.labels2[, r] <- rep(0, Ns)
      class.labels1[, r] <- rep(0, Ns)
      class.labels2[, r] <- rep(0, Ns)
      for (i in 1:Ns) {
        m1 <- sum(!is.na(match(label1.subset, i)))
        m2 <- sum(!is.na(match(full.subset, i)))
        reshuffled.class.labels1[i, r] <- m1
        reshuffled.class.labels2[i, r] <- m2 - m1
        if (i <= class1.size) {
          class.labels1[i, r] <- m2
          class.labels2[i, r] <- 0
        } else {
          class.labels1[i, r] <- 0
          class.labels2[i, r] <- m2
        }
      }
      
    } else if (permutation.type == 1) { # proportional (balanced) permutation
      
      class1.label1.subset <- sample(class1.subset, size = ceiling(class1.subset.size*fraction.class1))
      class2.label1.subset <- sample(class2.subset, size = floor(class2.subset.size*fraction.class1))
      reshuffled.class.labels1[, r] <- rep(0, Ns)
      reshuffled.class.labels2[, r] <- rep(0, Ns)
      class.labels1[, r] <- rep(0, Ns)
      class.labels2[, r] <- rep(0, Ns)
      for (i in 1:Ns) {
        if (i <= class1.size) {
          m1 <- sum(!is.na(match(class1.label1.subset, i)))
          m2 <- sum(!is.na(match(class1.subset, i)))
          reshuffled.class.labels1[i, r] <- m1
          reshuffled.class.labels2[i, r] <- m2 - m1
          class.labels1[i, r] <- m2
          class.labels2[i, r] <- 0
        } else {
          m1 <- sum(!is.na(match(class2.label1.subset, i)))
          m2 <- sum(!is.na(match(class2.subset, i)))
          reshuffled.class.labels1[i, r] <- m1
          reshuffled.class.labels2[i, r] <- m2 - m1
          class.labels1[i, r] <- 0
          class.labels2[i, r] <- m2
        }
      }
    }
  }
  
  # compute S2N for the random permutation matrix
  ##from each, randomly take X samples from the two class lables and get the mean as the null value of that gene's expression
  P <- reshuffled.class.labels1 * subset.mask
  n1 <- sum(P[,1])         
  M1 <- A %*% P
  M1 <- M1/n1      
  gc()
  A2 <- A*A        
  S1 <- A2 %*% P   
  S1 <- S1/n1 - M1*M1    
  S1 <- sqrt(abs((n1/(n1-1)) * S1))   
  gc()
  P <- reshuffled.class.labels2 * subset.mask
  n2 <- sum(P[,1])           
  M2 <- A %*% P           
  M2 <- M2/n2          
  gc()
  A2 <- A*A           
  S2 <- A2 %*% P      
  S2 <- S2/n2 - M2*M2 
  S2 <- sqrt(abs((n2/(n2-1)) * S2))
  rm(P)
  rm(A2)
  gc()
  
  if (sigma.correction == "GeneCluster") {  # small sigma "fix" as used in GeneCluster
    S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
    S2 <- ifelse(S2 == 0, 0.2, S2)
    S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
    S1 <- ifelse(S1 == 0, 0.2, S1)
    gc()
  }
  
  M1 <- M1 - M2
  rm(M2)
  gc()
  S1 <- S1 + S2
  rm(S2)
  gc()
  
  s2n.matrix <- M1/S1
  
  if (reverse.sign == T) {
    s2n.matrix <- - s2n.matrix
  }
  gc()
  
  for (r in 1:nperm) {
    order.matrix[, r] <- order(s2n.matrix[, r], decreasing=T)            
  }
  
  # compute S2N for the "observed" permutation matrix
  
  P <- class.labels1 * subset.mask
  n1 <- sum(P[,1])         
  M1 <- A %*% P
  M1 <- M1/n1      
  gc()
  A2 <- A*A        
  S1 <- A2 %*% P   
  S1 <- S1/n1 - M1*M1    
  S1 <- sqrt(abs((n1/(n1-1)) * S1))   
  gc()
  P <- class.labels2 * subset.mask
  n2 <- sum(P[,1])           
  M2 <- A %*% P           
  M2 <- M2/n2          
  gc()
  A2 <- A*A           
  S2 <- A2 %*% P      
  S2 <- S2/n2 - M2*M2 
  S2 <- sqrt(abs((n2/(n2-1)) * S2))
  rm(P)
  rm(A2)
  gc()
  
  if (sigma.correction == "GeneCluster") {  # small sigma "fix" as used in GeneCluster
    S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
    S2 <- ifelse(S2 == 0, 0.2, S2)
    S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
    S1 <- ifelse(S1 == 0, 0.2, S1)
    gc()
  } 
  
  M1 <- M1 - M2
  rm(M2)
  gc()
  S1 <- S1 + S2
  rm(S2)
  gc()
  
  obs.s2n.matrix <- M1/S1
  gc()
  
  if (reverse.sign == T) {
    obs.s2n.matrix <- - obs.s2n.matrix
  }
  
  for (r in 1:nperm) {
    obs.order.matrix[,r] <- order(obs.s2n.matrix[,r], decreasing=T)            
  }
  
  #rescale correl to (-1,1]
  for(k in 1:ncol(s2n.matrix)){
    x<-s2n.matrix[,k]
    s2n.matrix[,k]<-x/max(abs(x))
  }
  
  for(k in 1:ncol(obs.s2n.matrix)){
    x<-obs.s2n.matrix[,k]    
    obs.s2n.matrix[,k]<-x/max(abs(x))
  }
  
  if(reshuffling.type=='gene.labels'){
    for (r in 1:nperm) {
      order.matrix[,r]<- sample(1:nrow(obs.s2n.matrix))
      s2n.matrix[,r] <- obs.s2n.matrix[,1]
    }
  }
  
  o.m=cbind(obs.order.matrix[,1],order.matrix)
  s2n.m = cbind(obs.s2n.matrix[,1],s2n.matrix)
  s2n.m.sort.abs<-matrix(NA, nrow = nrow(s2n.m), ncol = ncol(s2n.m))
  for(i in 1:ncol(s2n.m)){
    s2n.m.sort.abs[,i]<-abs(s2n.m[,i][o.m[,i]])
  }
  
  
  
  
  return(list(#s2n.matrix = s2n.matrix, 
              #obs.s2n.matrix = obs.s2n.matrix, 
              #order.matrix = order.matrix,
              #obs.order.matrix = obs.order.matrix,
              o.m=o.m,
              s2n.m = s2n.m,
              s2n.m.sort.abs = s2n.m.sort.abs,
              reshuffled.class.labels1 = reshuffled.class.labels1))
}

GSEA.EnrichmentScore.weighted <- function(gene.list, gene.set, correl.vector = NULL,gene.labels,exp.type=1,scores=NULL) {  
  #
  # Computes the weighted GSEA score of gene.set in gene.list. 
  # The weighted score type is the exponent of the correlation 
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
  # necessary to input the correlation vector with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
  #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
  #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1) 
  #   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain") 
  #   RES: Numerical vector containing the running enrichment score for all locations in the gene list 
  #   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list 
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  #type = 0 equ step  1 use s2n 2 use weight
  
  #browser()
  gene.labels<-correl.vector
  correl.vector<-sort(correl.vector,decreasing = T)
  if(exp.type == 0){
    weights=1/scores
    t.i=c(0)
  }else if(exp.type == 1){
    weights=1/scores
    t.i=c(1)
  }else if(exp.type == 2){
    weights<-1/scores
    t.i<-match(gene.list, gene.set, nomatch=0)
    t.i[t.i>0]<-weights[t.i[t.i>0]]
  }
  names(weights)<-names(gene.labels[gene.set])
  
  tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
  no.tag.indicator <- 1 - tag.indicator 
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 

  correl.vector <- abs(correl.vector)
  
  sum.correl.tag    <- sum(tag.indicator * correl.vector^t.i)
  ###nomalize on gs.db
  norm.tag    <- 1.0/sum.correl.tag
  norm.no.tag <- 1.0/Nm
  
  
  RES <- cumsum(tag.indicator * correl.vector^t.i * norm.tag - no.tag.indicator * norm.no.tag)      
  
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (max.ES > - min.ES) {
    #      ES <- max.ES
    ES <- signif(max.ES, digits = 5)
    arg.ES <- which.max(RES)
  } else {
    #      ES <- min.ES
    ES <- signif(min.ES, digits=5)
    arg.ES <- which.min(RES)
  }
  #browser()
  ##debug  
  ##position_in_rank_list,gene_name,score_norm,corl,corl^score_norm
  a<-names(sort(gene.labels[gene.set],decreasing = T))
  names(a)<-match(a,names(sort(gene.labels,T)))
  #position in rank list / gene name / correl to class label / comfident score / norm step length(corl^p*norm)
  dbg<-paste(names(a),a,gene.labels[a],weights[a],correl.vector[a]^weights*norm.tag,sep='/')
  #names(dbg)<-names(a)
  #plot(c(0,RES,0),type='b')
  ##debug
  gc()
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator,dbg=dbg))    
}


GSEA.EnrichmentScore.weighted.batch <- function(gene.list, gene.set, correl.vector = NULL,correl.vector.sorted=NULL,gene.labels,exp.type=1,scores=NULL) {  
  #browser()

  correl.vector.sorted<-correl.vector.sorted
  # for(i in 1:ncol(correl.vector)){
  #   correl.vector.sorted[,i]<-abs(correl.vector[,i][gene.list[,i]])
  # }

  
  if(exp.type == 0){
    weights=1/scores
    t.i=matrix(rep(0,length(gene.list)),ncol = ncol(gene.list))
  }else if(exp.type == 1){
    weights=1/scores
    t.i=matrix(rep(1,length(gene.list)),ncol = ncol(gene.list))
  }else if(exp.type == 2){
    weights<-1/scores
    t.i<-apply(gene.list,2,match,table=gene.set,nomatch=0)
    t.i<-apply(t.i,2,function(x){x[x>0]<-weights[x[x>0]];x})
  }
  
  tag.indicator<-sign(apply(gene.list,2,match,table=gene.set,nomatch=0))
  no.tag.indicator <- 1 - tag.indicator
  N <- nrow(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 

  
  step<-tag.indicator * correl.vector.sorted^t.i
  sum.correl.tag <- colSums(step)
  
  ###nomalize on gs.db
  norm.tag    <- 1.0/sum.correl.tag
  norm.no.tag <- 1.0/Nm
  
  step<-sweep(step, 2, sum.correl.tag, FUN="/")
  
  RES <- vector(mode = "list",length(norm.tag))
  for(i in 1:length(norm.tag)){
    RES[[i]]<-cumsum(step[,i] - no.tag.indicator[,i] * norm.no.tag)   
    #names(RES[[i]])<-names(correl.vector[,i][gene.list[,i]])
  }
  
  max.ES <- sapply(RES,max)
  min.ES <- sapply(RES,min)
  
  ES=vector(mode = "numeric", length = length(max.ES))
  arg.ES=vector(mode = "numeric", length = length(max.ES))
  for(i in 1:length(max.ES)){
    if (max.ES[i] > - min.ES[i]) {
      ES[i] <- signif(max.ES[i], digits = 5)
      arg.ES[i] <- which.max(RES[[i]])
    } else {
      ES[i] <- signif(min.ES[i], digits=5)
      arg.ES[i] <- which.min(RES[[i]])
    }
  }
  
  
  #browser()
  ##debug  
  ##position_in_rank_list,gene_name,score_norm,corl,corl^score_norm

  # dbg<-vector(mode = "list", length = ncol(gene.list))
  # gene.labels<-rownames(correl.vector)
  # for(i in 1:ncol(correl.vector)){
  #   #gene.labels<-correl.vector[,i]
  #   #names(RES[[i]])<-names(gene.labels[gene.list[,i]])
  #   hit<-tag.indicator[,i]==1
  #   gname<-gene.labels[gene.list[,i][hit]]
  #   position<-which(hit)
  #   #a<-sort(gene.labels[gene.set],decreasing = T)
  #   #position<-match(names(a),names(gene.labels[gene.list[,i]]))
  #   dbg[[i]]<-paste(position,gname,correl.vector[,i][gene.list[,i][hit]],weights[gname],step[,i][position],sep='/')
  #   #dbg[[i]]<-paste(position,names(a),a,weights[names(a)],step[,i][position],sep='/')
  #   #names(dbg[[i]])=rep(c("position_in_rank_list/gene_name/correl_to_class/comfident_score/norm_step_length(corl^p*norm)"),length(dbg[[i]]))
  # }
  # #position in rank list / gene name / correl to class label / comfident score / norm step length(corl^p*norm)
  
  gene.labels<-rownames(correl.vector)
  hit<-tag.indicator[,1]==1
  gname<-gene.labels[gene.list[,1][hit]]
  position<-which(hit)
  #position in rank list / gene name / correl to class label / comfident score / norm step length(corl^p*norm)
  dbg<-paste(position,gname,correl.vector[,i][gene.list[,i][hit]],weights[gname],step[,i][position],sep='/')

  return(list(ES = ES, arg.ES = arg.ES, RES = RES[[1]], indicator = tag.indicator[,1],dbg=dbg))    
}

plot.result<-function(A=A,O=O,rl=rl,pty=pty,nom.p.val.threshold=0.05,
                      output.directory='~/Desktop/tmp/2/',doc.string=NULL,topgs=10,adjust.param=0.5){
  
  
  
  if(is.null(doc.string)){
    doc.string='gsea'
  }
  
  #browser()
  pdf(file=paste(output.directory, doc.string, ".global.plots", sep="", collapse=""), height = 10, width = 10)
  nf <- layout(matrix(c(1,1,2,3), 2, 2, byrow=T), c(1,1), c(1,1), TRUE)
  
  Ng<-length(rl)
  Ns<-ncol(A)
  
  obs.s2n=sort(O$s2n.m[,1],T)
  N=length(obs.s2n)
  location <- 1:N
  max.corr <- max(obs.s2n)
  min.corr <- min(obs.s2n)
  phen=pty$phen
  
  x <- plot(location, obs.s2n, ylab = "Signal to Noise Ratio (S2N)", xlab = "Gene List Location", main = "Gene List Correlation (S2N) Profile", type = "l", lwd = 2, cex = 0.9, col = 1)            
  
  for (i in seq(1, N, 20)) {
    lines(c(i, i), c(0, obs.s2n[i]), lwd = 3, cex = 0.9, col = colors()[12]) # shading of correlation plot
  }
  x <- points(location, obs.s2n, type = "l", lwd = 2, cex = 0.9, col = 1)            
  lines(c(1, N), c(0, 0), lwd = 2, lty = 1, cex = 0.9, col = 1) # zero correlation horizontal line
  temp <- order(abs(obs.s2n), decreasing=T)
  arg.correl <- temp[N]
  lines(c(arg.correl, arg.correl), c(min.corr, 0.7*max.corr), lwd = 2, lty = 3, cex = 0.9, col = 1) # zero correlation vertical line
  
  area.bias <- signif(100*(sum(obs.s2n[1:arg.correl]) + sum(obs.s2n[arg.correl:N]))/sum(abs(obs.s2n[1:N])), digits=3)
  area.phen <- ifelse(area.bias >= 0, phen[1], phen[2])
  delta.string <- paste("Corr. Area Bias to \"", area.phen, "\" =", abs(area.bias), "%", sep="", collapse="")
  zero.crossing.string <- paste("Zero Crossing at location ", arg.correl, " (",  signif(100*arg.correl/N, digits=3), " %)")
  leg.txt <- c(delta.string, zero.crossing.string)
  legend(x=N/10, y=max.corr, bty="n", bg = "white", legend=leg.txt, cex = 0.9)
  
  leg.txt <- paste("\"", phen[1], "\" ", sep="", collapse="")
  text(x=1, y=-0.05*max.corr, adj = c(0, 1), labels=leg.txt, cex = 0.9)
  
  leg.txt <- paste("\"", phen[2], "\" ", sep="", collapse="")
  text(x=N, y=0.05*max.corr, adj = c(1, 0), labels=leg.txt, cex = 0.9)
  
  
  ###########
  ES.norm=t(sapply(rl,function(x){
    x$ES.premut.norm
    }))
  obs.ES.norm<-sapply(rl,function(x){
    x$ES.obs.norm
  })
  obs.ES<-sapply(rl,function(x){x$ES.obs})
  nperm<-ncol(O$reshuffled.class.labels1)
  if(Ng>0){
    phi.densities.pos <- matrix(0, nrow=512, ncol=nperm)
    phi.densities.neg <- matrix(0, nrow=512, ncol=nperm)
    obs.phi.densities.pos <- matrix(0, nrow=512, ncol=nperm)
    obs.phi.densities.neg <- matrix(0, nrow=512, ncol=nperm)
    phi.density.mean.pos <- vector(length=512, mode = "numeric")
    phi.density.mean.neg <- vector(length=512, mode = "numeric")
    obs.phi.density.mean.pos <- vector(length=512, mode = "numeric")
    obs.phi.density.mean.neg <- vector(length=512, mode = "numeric")
    phi.density.median.pos <- vector(length=512, mode = "numeric")
    phi.density.median.neg <- vector(length=512, mode = "numeric")
    obs.phi.density.median.pos <- vector(length=512, mode = "numeric")
    obs.phi.density.median.neg <- vector(length=512, mode = "numeric")
    x.coor.pos <-  vector(length=512, mode = "numeric")
    x.coor.neg <-  vector(length=512, mode = "numeric")
    #browser()
    for (i in 1:nperm) {
      pos.phi <- ES.norm[ES.norm[, i] >= 0, i]
      if (length(pos.phi) > 2) {
        temp <- density(pos.phi, adjust=adjust.param, n = 512, from=0, to=3.5)
      } else {
        temp <- list(x = 3.5*(seq(1, 512) - 1)/512, y = rep(0.001, 512))
      }
      phi.densities.pos[, i] <- temp$y
      norm.factor <- sum(phi.densities.pos[, i])
      phi.densities.pos[, i] <- phi.densities.pos[, i]/norm.factor
      if (i == 1) {
        x.coor.pos <- temp$x
      }
      
      neg.phi <- ES.norm[ES.norm[, i] < 0, i]
      if (length(neg.phi) > 2) {
        temp <- density(neg.phi, adjust=adjust.param, n = 512, from=-3.5, to=0)
      } else {
        temp <- list(x = 3.5*(seq(1, 512) - 1)/512, y = rep(0.001, 512))
      }
      phi.densities.neg[, i] <- temp$y
      norm.factor <- sum(phi.densities.neg[, i])
      phi.densities.neg[, i] <- phi.densities.neg[, i]/norm.factor
      if (i == 1) {
        x.coor.neg <- temp$x
      }
      
      
      pos.phi <- obs.ES.norm[obs.ES.norm >= 0]
      if (length(pos.phi) > 2) {
        temp <- density(pos.phi, adjust=adjust.param, n = 512, from=0, to=3.5)
      } else {
        temp <- list(x = 3.5*(seq(1, 512) - 1)/512, y = rep(0.001, 512))
      }
      obs.phi.densities.pos[, i] <- temp$y
      norm.factor <- sum(obs.phi.densities.pos[, i])
      obs.phi.densities.pos[, i] <- obs.phi.densities.pos[, i]/norm.factor
      
      neg.phi <- obs.ES.norm[obs.ES.norm < 0]
      if (length(neg.phi)> 2) {  
        temp <- density(neg.phi, adjust=adjust.param, n = 512, from=-3.5, to=0)
      } else {
        temp <- list(x = 3.5*(seq(1, 512) - 1)/512, y = rep(0.001, 512))
      }
      obs.phi.densities.neg[, i] <- temp$y
      norm.factor <- sum(obs.phi.densities.neg[, i])
      obs.phi.densities.neg[, i] <- obs.phi.densities.neg[, i]/norm.factor
    }
    
    phi.density.mean.pos <- apply(phi.densities.pos, 1, mean)
    phi.density.mean.neg <- apply(phi.densities.neg, 1, mean)
    
    obs.phi.density.mean.pos <- apply(obs.phi.densities.pos, 1, mean)
    obs.phi.density.mean.neg <- apply(obs.phi.densities.neg, 1, mean)
    
    phi.density.median.pos <- apply(phi.densities.pos, 1, median)
    phi.density.median.neg <- apply(phi.densities.neg, 1, median)
    
    obs.phi.density.median.pos <- apply(obs.phi.densities.pos, 1, median)
    obs.phi.density.median.neg <- apply(obs.phi.densities.neg, 1, median)
    
    x <- c(x.coor.neg, x.coor.pos)
    x.plot.range <- range(x)
    y1 <- c(phi.density.mean.neg, phi.density.mean.pos)
    y2 <- c(obs.phi.density.mean.neg, obs.phi.density.mean.pos)
    y.plot.range <- c(-0.3*max(c(y1, y2)),  max(c(y1, y2)))
    
    #print(c(y.plot.range, max(c(y1, y2)), max(y1), max(y2)))
    
    plot(x, y1, xlim = x.plot.range, ylim = 1.5*y.plot.range, type = "l", lwd = 2, col = 2, xlab = "NES", ylab = "P(NES)", main = "Global Observed and Null Densities (Area Normalized)")
    
    y1.point <- y1[seq(1, length(x), 2)]
    y2.point <- y2[seq(2, length(x), 2)]
    x1.point <- x[seq(1, length(x), 2)]
    x2.point <- x[seq(2, length(x), 2)]
    points(x, y1, type = "l", lwd = 2, col = colors()[555])
    points(x, y2, type = "l", lwd = 2, col = colors()[29])
    
    for (i in 1:Ng) {
      col <- ifelse(obs.ES.norm[i] > 0, 2, 3) 
      lines(c(obs.ES.norm[i], obs.ES.norm[i]), c(-0.2*max(c(y1, y2)), 0), lwd = 1, lty = 1, col = 1)
    }
    leg.txt <- paste("Neg. ES: \"", phen[2], " \" ", sep="", collapse="")
    text(x=x.plot.range[1], y=-0.25*max(c(y1, y2)), adj = c(0, 1), labels=leg.txt, cex = 0.9)
    leg.txt <- paste(" Pos. ES: \"", phen[1], "\" ", sep="", collapse="")
    text(x=x.plot.range[2], y=-0.25*max(c(y1, y2)), adj = c(1, 1), labels=leg.txt, cex = 0.9)
    
    leg.txt <- c("Null Density", "Observed Density", "Observed NES values")
    c.vec <- c(colors()[555], colors()[29], 1)
    lty.vec <- c(1, 1, 1)
    lwd.vec <- c(2, 2, 2)
    legend(x=0, y=1.5*y.plot.range[2], bty="n", bg = "white", legend=leg.txt, lty = lty.vec, lwd = lwd.vec, col = c.vec, cex = 0.9)
    
    
    B <- A[order(obs.s2n,decreasing = T),]
    if (N > 300) {
      C <- rbind(B[1:100,], rep(0, Ns), rep(0, Ns), B[(floor(N/2) - 50 + 1):(floor(N/2) + 50),], rep(0, Ns), rep(0, Ns), B[(N - 100 + 1):N,])
    } 
    rm(B)
    GSEA.HeatMapPlot(V = C, col.labels = pty$class.v, col.classes = pty$phen, main = "Heat Map for Genes in Dataset")
    dev.off()
  }
  
  
  
  p.vals=sapply(rl,function(x){x$p})
  
  tag.frac<-sapply(rl,function(x){
    if(x$ES.obs>0)
      sum(x$GSEA.results$indicator[1:x$GSEA.results$arg.ES[1]])/sum(x$GSEA.results$indicator)
    else
      sum(x$GSEA.results$indicator[x$GSEA.results$arg.ES[1]:length(x$GSEA.results$indicator)])/sum(x$GSEA.results$indicator)
  })
  
  gene.frac<-sapply(rl,function(x){
    if(x$ES.obs>0)
      x$GSEA.results$arg.ES[1]/length(x$GSEA.results$indicator)
    else
      1-x$GSEA.results$arg.ES[1]/length(x$GSEA.results$indicator)
  })
  
  
  #browser()
  report1 <- data.frame(cbind(names(rl), Term(ONTTERM)[names(rl)],sapply(rl,function(x){sum(x$GSEA.results$indicator)}),  round(obs.ES,3), 
                              round(obs.ES.norm,3), round(p.vals,3), round(tag.frac,3), round(gene.frac,3)),row.names = NULL,stringsAsFactors = F)
  names(report1) <- c("GS","DEF","SIZE", "ES", "NES", "p",  "Tag %", "Gene %")
  #       print(report)
  report2 <- report1
  report.index2 <- order(obs.ES.norm, decreasing=T)
  for (i in 1:Ng) {
    report2[i,] <- report1[report.index2[i],]
  }   
  report3 <- report1
  report.index3 <- order(obs.ES.norm, decreasing=F)
  for (i in 1:Ng) {
    report3[i,] <- report1[report.index3[i],]
  }   
  phen1.rows <- length(obs.ES.norm[obs.ES.norm >= 0])
  phen2.rows <- length(obs.ES.norm[obs.ES.norm < 0])
  report.phen1 <- report2[1:phen1.rows,]
  report.phen2 <- report3[1:phen2.rows,]
  
  if (output.directory != "")  {
    if (phen1.rows > 0) {
      filename <- paste(output.directory, doc.string, ".SUMMARY.RESULTS.REPORT.", pty$phen[1],".txt", sep="", collapse="")
      write.table(report.phen1, file = filename, quote=F, row.names=F, sep = "\t")
    }
    if (phen2.rows > 0) {
      filename <- paste(output.directory, doc.string, ".SUMMARY.RESULTS.REPORT.", pty$phen[2],".txt", sep="", collapse="")
      write.table(report.phen2, file = filename, quote=F, row.names=F, sep = "\t")
    }
  }
  
  ####plot for each gene set
  if (topgs > floor(Ng/2)) {
    topgs <- floor(Ng/2)
  }
  
  #browser()
  result<-list()
  result[['report']][[pty$phen[1]]]<-report.phen1
  result[['report']][[pty$phen[2]]]<-report.phen2
  #browser()
  for (i in 1:Ng) {
    # result[[names(rl)[i]]][['p']]<-p.vals[i]
    # result[[names(rl)[i]]][['ES']]<-obs.ES[i]
    # result[[names(rl)[i]]][['ES.norm']]<-obs.ES.norm[i]
    # if(names(rl)[i]=='DOID:0014667'){
    #   browser()
    #   1
    # }
    if (p.vals[i] <= nom.p.val.threshold ||
        (is.element(i, c(order(obs.ES.norm,decreasing = T)[1:topgs], order(obs.ES.norm,decreasing = T)[(Ng - topgs + 1): Ng]))))
      {
      #  produce report per gene set
      obs.index<-which(rl[[i]]$GSEA.results$indicator==1)
      gene.s2n<-obs.s2n[obs.index]
      gene.names<-names(gene.s2n)
      gene.RES<-rl[[i]]$GSEA.results$RES[obs.index]
      loc <- match(i,order(obs.ES.norm,decreasing = T))
      if (rl[[i]]$'ES.obs' >= 0) {
        core.enrichment<-obs.index <= rl[[i]]$GSEA.results$arg.ES[1]
        loc <- match(i,order(obs.ES.norm,decreasing = T))
        phen.tag <- phen[1]
      } else {
        core.enrichment<-obs.index >= rl[[i]]$GSEA.results$arg.ES[1]
        loc <- Ng - match(i,order(obs.ES.norm,decreasing = T)) + 1
        phen.tag <- phen[2]
      }
      
     #browser()
      gene.report <- data.frame(cbind(c(1:length(gene.names)), gene.names, obs.index, gene.s2n, 
                                      round(1/as.numeric(sapply(rl[[i]]$GSEA.results$dbg,function(x){strsplit(x,split = '/',fixed = T)[[1]][4]}))),
                                      sapply(rl[[i]]$GSEA.results$dbg,function(x){strsplit(x,split = '/',fixed = T)[[1]][4]}),
                                      sapply(rl[[i]]$GSEA.results$dbg,function(x){strsplit(x,split = '/',fixed = T)[[1]][5]}),
                                      gene.RES, core.enrichment),row.names = NULL,stringsAsFactors = F)
      names(gene.report) <- c("gene.number", "GENE","LIST LOC", "S2N","CS","expont p","STEP","RES","CORE_ENRICHMENT")
      result[[phen.tag]][[names(rl)[i]]][['report']]<-gene.report
      result[[phen.tag]][[names(rl)[i]]][['p']]<-p.vals[[i]]
      result[[phen.tag]][[names(rl)[i]]][['ES']]<-obs.ES[[i]]
      result[[phen.tag]][[names(rl)[i]]][['ES.norm']]<-obs.ES.norm[[i]]
      #       print(gene.report)
      
      if (output.directory != "")  {
        
        filename <- paste(output.directory, doc.string, ".", phen.tag, ".", loc,".", names(rl)[i], ".report.",  ".txt", sep="", collapse="")
        write.table(gene.report, file = filename, quote=F, row.names=F, sep = "\t")
        
        gs.filename <- paste(output.directory, doc.string, ".", phen.tag, ".", loc,".", names(rl)[i], ".plot", ".pdf", sep="", collapse="")
        pdf(file=gs.filename, height = 6, width = 14)
        result[[phen.tag]][[names(rl)[i]]][['plot']]<-gs.filename
      }
      
      nf <- layout(matrix(c(1,2,3), 1, 3, byrow=T), 1, c(1, 1, 1), TRUE)
      
      ind <- 1:N
      min.RES <- min(rl[[i]]$GSEA.results$RES)
      max.RES <- max(rl[[i]]$GSEA.results$RES)
      if (max.RES < 0.3) max.RES <- 0.3
      if (min.RES > -0.3) min.RES <- -0.3
      delta <- (max.RES - min.RES)*0.50
      min.plot <- min.RES - 2*delta
      max.plot <- max.RES
      max.corr <- max(obs.s2n)
      min.corr <- min(obs.s2n)
      Obs.correl.vector.norm <- (obs.s2n - min.corr)/(max.corr - min.corr)*1.25*delta + min.plot
      zero.corr.line <- (- min.corr/(max.corr - min.corr))*1.25*delta + min.plot
      col <- ifelse(obs.ES[i] > 0, 2, 4)
      
      # Running enrichment plot
      
      sub.string <- paste("Number of genes: ", N, " (in list), ", length(gene.names), " (in gene set)", sep = "", collapse="")
      
      main.string <- paste("Gene Set ", i, ":", names(rl)[i])
      plot(ind, rl[[i]]$GSEA.results$RES, main = main.string, sub = sub.string, xlab = "Gene List Index", ylab = "Running Enrichment Score (RES)", xlim=c(1, N), ylim=c(min.plot, max.plot), type = "l", lwd = 2, cex = 1, col = col)
      for (j in seq(1, N, 20)) {
        lines(c(j, j), c(zero.corr.line, Obs.correl.vector.norm[j]), lwd = 1, cex = 1, col = colors()[12]) # shading of correlation plot
      }
      lines(c(1, N), c(0, 0), lwd = 1, lty = 2, cex = 1, col = 1) # zero RES line
      lines(c(rl[[i]]$GSEA.results$arg.ES[1], rl[[i]]$GSEA.results$arg.ES[1]), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = col) # max enrichment vertical line
      for(j in obs.index){
        lines(c(j, j), c(min.plot + 1.25*delta, min.plot + 1.75*delta), lwd = 1, lty = 1, cex = 1, col = 1)  # enrichment tags
      }

      lines(ind, Obs.correl.vector.norm, type = "l", lwd = 1, cex = 1, col = 1)
      lines(c(1, N), c(zero.corr.line, zero.corr.line), lwd = 1, lty = 1, cex = 1, col = 1) # zero correlation horizontal line
      temp <- order(abs(obs.s2n), decreasing=T)
      arg.correl <- temp[N]
      lines(c(arg.correl, arg.correl), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = 3) # zero crossing correlation vertical line
      
      leg.txt <- paste("\"", phen[1], "\" ", sep="", collapse="")
      text(x=1, y=min.plot, adj = c(0, 0), labels=leg.txt, cex = 1.0)
      
      leg.txt <- paste("\"", phen[2], "\" ", sep="", collapse="")
      text(x=N, y=min.plot, adj = c(1, 0), labels=leg.txt, cex = 1.0)
      
      adjx <- ifelse(obs.ES[i] > 0, 0, 1)
      
      leg.txt <- paste("Peak at ", rl[[i]]$GSEA.results$arg.ES[1], sep="", collapse="")
      text(x=rl[[i]]$GSEA.results$arg.ES[1], y=min.plot + 1.8*delta, adj = c(adjx, 0), labels=leg.txt, cex = 1.0)
      
      leg.txt <- paste("Zero crossing at ", arg.correl, sep="", collapse="")
      text(x=arg.correl, y=min.plot + 1.95*delta, adj = c(adjx, 0), labels=leg.txt, cex = 1.0)
      
      # nominal p-val histogram
      
      sub.string <- paste("ES =", signif(obs.ES[i], digits = 3), " NES =", signif(obs.ES.norm[i], digits=3), " Nom. p-val=", signif(rl[[i]]$p, digits = 3))
      temp <- density(rl[[i]]$ES.premut, adjust=adjust.param)
      x.plot.range <- range(temp$x)
      y.plot.range <- c(-0.125*max(temp$y), 1.5*max(temp$y))
      plot(temp$x, temp$y, type = "l", sub = sub.string, xlim = x.plot.range, ylim = y.plot.range, lwd = 2, col = 2, main = "Gene Set Null Distribution", xlab = "ES", ylab="P(ES)")
      x.loc <- which.min(abs(temp$x - obs.ES[i]))
      lines(c(obs.ES[i], obs.ES[i]), c(0, temp$y[x.loc]), lwd = 2, lty = 1, cex = 1, col = 1)
      lines(x.plot.range, c(0, 0), lwd = 1, lty = 1, cex = 1, col = 1)
      
      leg.txt <- c("Gene Set Null Density", "Observed Gene Set ES value")
      c.vec <- c(2, 1)
      lty.vec <- c(1, 1)
      lwd.vec <- c(2, 2)
      legend(x=-0.2, y=y.plot.range[2], bty="n", bg = "white", legend=leg.txt, lty = lty.vec, lwd = lwd.vec, col = c.vec, cex = 1.0)
      
      leg.txt <- paste("Neg. ES \"", phen[2], "\" ", sep="", collapse="")
      text(x=x.plot.range[1], y=-0.1*max(temp$y), adj = c(0, 0), labels=leg.txt, cex = 1.0)
      leg.txt <- paste(" Pos. ES: \"", phen[1], "\" ", sep="", collapse="")
      text(x=x.plot.range[2], y=-0.1*max(temp$y), adj = c(1, 0), labels=leg.txt, cex = 1.0)
      
      # create pinkogram for each gene set
      pinko<-A[names(obs.s2n[obs.index]),]
      ##in case only one gene in the set!
      if(is.null(rownames(pinko))){
        pinko<-t(as.matrix(pinko))
        rownames(pinko)<-names(obs.s2n[obs.index])
      }
      ##
      #browser()
      GSEA.HeatMapPlot(V = pinko, row.names = rownames(pinko), col.labels = pty$class.v, col.classes = pty$phen, col.names = colnames(pinko), main =" Heat Map for Genes in Gene Set", xlab=" ", ylab=" ")
      dev.off()

    } # if p.vals thres
    
  } # loop over gene sets
  return(result)
}



GSEA.HeatMapPlot <- function(V, row.names = F, col.labels, col.classes, col.names = F, main = " ", xlab=" ", ylab=" ") {
  #
  # Plots a heatmap "pinkogram" of a gene expression matrix including phenotype vector and gene, sample and phenotype labels
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  n.rows <- length(V[,1])
  n.cols <- length(V[1,])
  row.mean <- apply(V, MARGIN=1, FUN=mean)
  row.sd <- apply(V, MARGIN=1, FUN=sd)
  row.n <- length(V[,1])
  for (i in 1:n.rows) {
    if (row.sd[i] == 0) {
      V[i,] <- 0
    } else {
      V[i,] <- (V[i,] - row.mean[i])/(0.5 * row.sd[i])
    }
    V[i,] <- ifelse(V[i,] < -6, -6, V[i,])
    V[i,] <- ifelse(V[i,] > 6, 6, V[i,])
  }
  
  mycol <- c("#0000FF", "#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D", "#FF0000") # blue-pinkogram colors. The first and last are the colors to indicate the class vector (phenotype). This is the 1998-vintage, pre-gene cluster, original pinkogram color map
  
  mid.range.V <- mean(range(V)) - 0.1
  heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
  heatm[1:n.rows,] <- V[seq(n.rows, 1, -1),]
  heatm[n.rows + 1,] <- ifelse(col.labels == 0, 7, -7)
  image(1:n.cols, 1:(n.rows + 1), t(heatm), col=mycol, axes=FALSE, main=main, xlab= xlab, ylab=ylab)
  
  if (length(row.names) >= 1) {
    numC <- nchar(row.names)
    size.row.char <- 35/(n.rows + 5)
    size.col.char <- 25/(n.cols + 5)
    maxl <- floor(n.rows/1.6)
    for (i in 1:n.rows) {
      row.names[i] <- substr(row.names[i], 1, maxl)
    }
    row.names <- c(row.names[seq(n.rows, 1, -1)], "Class")
    axis(2, at=1:(n.rows + 1), labels=row.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.row.char, font.axis=2, line=-1)
  }
  
  if (length(col.names) > 1) {
    axis(1, at=1:n.cols, labels=col.names, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=2, line=-1)
  }
  
  C <- split(col.labels, col.labels)
  class1.size <- length(C[[1]])
  class2.size <- length(C[[2]])
  axis(3, at=c(floor(class1.size/2),class1.size + floor(class2.size/2)), labels=col.classes, tick=FALSE, las = 1, cex.axis=1.25, font.axis=2, line=-1)
  
  return()
}

read.gmt<-function(file='/home/xin/Downloads/GSEA-P-R/GeneSetDatabases/HDO_o.gmt'){
  l<-list()
  gs.db<-readLines(file)
  for (i in 1:length(gs.db)) {
    temp<-strsplit(gs.db[[i]], "\t")[[1]]
    l[[temp[1]]]<-temp[-c(1,2)]
  }
  l
}

write.gmt<-function(term2geneID,file){
  def<-Term(ONTTERM)
  s=''
  for(i in 1:length(term2geneID)){
    n<-names(term2geneID)[i]
    s<-paste(s,paste(n,def[n],paste(term2geneID[[i]],collapse = '\t'),sep = '\t'),'\n',sep='')
  }
  sink(file)
  cat(s)
  sink()
}

