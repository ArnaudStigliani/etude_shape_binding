rm(list=ls())
library(Biostrings)
nRegion <- 600
library(stringr)
source("PNfonctions.r")                 # fonctions auxiliaires



#-------------------------------------read matrices ------------------------------------------


pfm_ARF2<- read.table("m_ARF2.txt",header=TRUE,sep="\t",skip=1)
pfm_ARF2 <- round((t(as.matrix(pfm_ARF2)))*nRegion)+1 ;pfm_ARF2
maxi_ARF2 <- apply(pfm_ARF2,FUN=max, 2)
maxi_ARF2 <- matrix(nrow=4, rep(maxi_ARF2,4),byrow=TRUE)
pwm_ARF2 <- log(pfm_ARF2/maxi_ARF2)
pwm_ARF2_rev <- pwm_ARF2 - minScore(pwm_ARF2)/dim(pwm_ARF2)[2] 

pwm_ARF2 <-  reverseComplement(pwm_ARF2_rev) ; pwm_ARF2

#-------------------------------------read fasta files-----------------------------------------


ARF2_pos <- readDNAStringSet('ARF2.fas')#[(1:1000)]
width_pos <- width(ARF2_pos)
seq_pos <- as.character(ARF2_pos)



#-------------------------------------Compute Scores-----------------------------------------
th <- maxScore(pwm_ARF2) - 9

#pos

scores_ARF2_pos<- mapply(seq_pos,FUN=PWMscoreStartingAt,SIMPLIFY=FALSE,  starting.at=mapply(seq,1,width_pos-dim(pwm_ARF2)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF2))

scores_ARF2_rev_pos <- mapply(seq_pos,FUN=PWMscoreStartingAt,SIMPLIFY=FALSE,  starting.at=mapply(seq,1,width_pos-dim(pwm_ARF2_rev)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF2_rev))

density_pos <- (sapply(FUN=sum,lapply(FUN=">",scores_ARF2_pos,th))+sapply(FUN=sum,lapply(FUN=">",scores_ARF2_rev_pos,th)))/(width_pos*2)

