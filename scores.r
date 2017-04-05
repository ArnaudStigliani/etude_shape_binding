rm(list=ls())
library(Biostrings)
nRegion <- 600
library(stringr)




#-------------------------------------read matrices ------------------------------------------


pfm_ARF<- read.table("m_ARF5.txt",header=TRUE,sep="\t",skip=1)
pfm_ARF <- round((t(as.matrix(pfm_ARF)))*nRegion)+1 ;pfm_ARF
pfm_ARF <- cbind(rep(151,4),rep(151,4),pfm_ARF)
maxi_ARF <- apply(pfm_ARF,FUN=max, 2)
maxi_ARF <- matrix(nrow=4, rep(maxi_ARF,4),byrow=TRUE)
pwm_ARF <- log(pfm_ARF/maxi_ARF)
pwm_ARF_rev <- pwm_ARF - minScore(pwm_ARF)/dim(pwm_ARF)[2] 
pwm_ARF <-  reverseComplement(pwm_ARF_rev) ; pwm_ARF

#-------------------------------------read fasta files-----------------------------------------


ARF_pos <- readDNAStringSet('ARF5.fas')#[(1:1000)]
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
seq_plus <- seq_pos[pos_max==96]
seq_rev<- seq_rev_pos[pos_max_rev==96]
#
scores_plus <- scores_ARF_pos[,pos_max==96]
scores_rev <- scores_ARF_rev_pos[,pos_max_rev==96]
#
seq <- c(seq_plus,seq_rev)
scores <- cbind(scores_plus,scores_rev)
#
seq_ind <- apply(FUN=max,scores,2) > th & apply(FUN=max,scores,2) < th + 1
seq_sel <- seq[seq_ind]


###########################C ompute Dinucleotide ###############################

pwm <- matrix(c(1,5,9,13,0,1,2,3),byrow=FALSE,nrow=4)
rownames(pwm) <- c("A","C","G","T")

#----------------------- read table ------------------------------------

dinuc <- read.table("table.txt",header=TRUE,sep="\t")

#----------------------- read fasta ------------------------------------

reg_size <- mean(width(ARF))

#------------------------- compute  Dinucleotide  -------------------------

seq_sel<- as.character(ARF)


scores <- as.data.frame(lapply(seq,FUN=PWMscoreStartingAt,pwm=pwm,starting.at=(1:(reg_size-1))))
dim_score<- dim(scores)
scores <- matrix(as.vector(unlist(scores)),nrow=dim_score[1])

tableau <- matrix(0,125,reg_size-1)
rownames(tableau) <- dinuc[,2]
for (i in 1:(dim_score[2]))
{
    for (j in 1:(dim_score[1]))
    {
        tableau[,j] <- tableau[,j] + dinuc[,scores[j,i]+2]
    }
}
tableau <- tableau/dim_score[2]


names <- rownames(tableau)
tableau <- tableau[!sapply(names,FUN=str_detect,pattern="RNA"),]


