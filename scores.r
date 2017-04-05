rm(list=ls())
library(Biostrings)
nRegion <- 600
library(stringr)




#-------------------------------------read matrices ------------------------------------------


pfm_ARF<- read.table("m_ARF2.txt",header=TRUE,sep="\t",skip=1)
pfm_ARF <- round((t(as.matrix(pfm_ARF)))*nRegion)+1 ;pfm_ARF
maxi_ARF <- apply(pfm_ARF,FUN=max, 2)
maxi_ARF <- matrix(nrow=4, rep(maxi_ARF,4),byrow=TRUE)
pwm_ARF <- log(pfm_ARF/maxi_ARF)

pwm_ARF_rev <- pwm_ARF - minScore(pwm_ARF)/dim(pwm_ARF)[2] 

pwm_ARF <-  reverseComplement(pwm_ARF_rev) ; pwm_ARF

#-------------------------------------read fasta files-----------------------------------------


ARF_pos <- readDNAStringSet('ARF2.fas')#[(1:1000)]
width_pos <- width(ARF_pos)
seq_pos <- as.character(ARF_pos)
seq_rev_pos <- as.character(reverseComplement(ARF_pos))



#-------------------------------------Compute Scores-----------------------------------------
th <- -3
#
scores_ARF_pos<- mapply(seq_pos,FUN=PWMscoreStartingAt,SIMPLIFY=TRUE,  starting.at=mapply(seq,1,width_pos-dim(pwm_ARF)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF)) - maxScore(pwm_ARF)
scores_ARF_rev_pos<- mapply(seq_rev_pos,FUN=PWMscoreStartingAt,SIMPLIFY=TRUE,  starting.at=mapply(seq,1,width_pos-dim(pwm_ARF)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF)) - maxScore(pwm_ARF)
#
pos_max <- apply(FUN=which.max,scores_ARF_pos,2)
pos_max_rev<- apply(FUN=which.max,scores_ARF_rev_pos,2)
#
seq_plus <- seq_pos[pos_max==93]
seq_rev<- seq_rev_pos[pos_max_rev==93]
#
scores_plus <- scores_ARF_pos[,pos_max==93]
scores_rev <- scores_ARF_rev_pos[,pos_max_rev==93]
#
seq <- c(seq_plus,seq_rev)
scores <- c(scores_plus,scores_rev)
#
seq_sel <- apply(FUN=max,scores_rev,2) > th & apply(FUN=max,scores_rev,2) < th + 1
