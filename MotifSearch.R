library(Biostrings)
library(seqinr)

rm(list = ls())
setwd('C:/Users/hanic/Desktop/PRG/exercise_07')


############################## MOTIF SEARCH ####################################
##################### The Brute-Force Motif Search #############################

################################################################################
############ TASK 1 - Counting consensus score for consensus string ############
Score <- function(start_indexes, sequences, motif_length){
  num_seqs <- length(sequences)
  score <- 0
  
  for (pos in 0:(motif_length - 1)){
    column_bases <- character(num_seqs)
    
    for (seq in 1:num_seqs){
      start_pos <- start_indexes[seq]
      base <- subseq(sequences[[seq]], start = start_pos + pos, width = 1)
      column_bases[seq] <- as.character(base)
    }
    freq_table <- table(column_bases)
    score <- score + max(freq_table)
  }
  return(score)
}  


ScoreModified <- function(start_indexes, i, sequences, motif_length){
  num_seqs <- length(sequences)
  score <- 0
  
  for (pos in 0:(motif_length - 1)){
    column_bases <- character(num_seqs)
    
    for (s in 1:i){
      start_pos <- start_indexes[s]
      base <- subseq(sequences[[s]], start = start_pos + pos, width = 1)
      column_bases[s] <- as.character(base)
    }
    freq_table <- table(column_bases)
    score <- score + max(freq_table)
  }
  return(score)
}  
  
  
# seqs <- readDNAStringSet('seq_score.fasta', format = 'fasta')
# idxs <- c(1,1,1,1,1)
# motif_length <- 5
# Score(idxs, seqs, motif_length)


################################################################################
############## TASK 2 - creating array of starting indexes #####################
NextLeaf <- function(s, t, k){
  for (i in t:1){
    if (s[i] < k){
      s[i] <- s[i] + 1
      return(s)
    }
    s[i] <- 1
  }
  return(s)
}


# s = c(1,1,1,1,1) # array of starting indexes
# t = length(seqs)   # number of sequences
# n = width(seqs)[1]  # length of sequences
# l = 6   # motif_length
# k = n - l + 1
# NextLeaf(s, t, k)


####################################################################################
# TASK 3 - searching for array of starting positions with best score for consensus #
BFMotifSearch <- function(DNA, t, n, l){
  s <- c(1, 1, 1, 1)
  bestScore <- Score(s, DNA, l)
  while (TRUE){
    s <- NextLeaf(s, t, (n - l + 1))
    if (Score(s, DNA, l) > bestScore){
      bestScore <- Score(s, DNA, l)
      bestMotif <- c(s)
    }
    if (all( s == rep(1, t))){
      return(bestMotif)
    }
  }
}

#sequences <- readDNAStringSet('seq_score.fasta', format = 'fasta')
DNA <- readDNAStringSet('seq_motif.fasta', format = 'fasta')
motif_length <- 5
k <- width(DNA)[1] + 1 - motif_length
s = c(1,1,1,1)
t <- length(DNA) # number of sequences
n <- width(DNA)[1] # length of each sequence
l <- 5 # motif length
BFMotifSearch(DNA, t, n, l)


################################################################################
##################### The Branch-and-Bound Motif Search ########################

################################################################################
###################### TASK 4 - Next Vertex in a tree ##########################
NextVertex <- function(s, i, t, k){
  if (i < t){
    s[i+1] <- 1
    ret <- c(s, (i + 1))
    return(ret)
  }
  else{
    for (j in t:1){
      if (s[j] < k){
        s[j] <- s[j] + 1
        ret <- c(s,j)
        return(ret)
      }
    }
  }
  ret <- c(s,0)
  return(ret)
}


# sequences <- readDNAStringSet('seq_motif.fasta', format = 'fasta')
# s <- c(1,1,1,1)
# i <- 2
# t <- length(sequences) # number of sequences
# n <- width(sequences)[1] # length of each sequence
# l <- 5 # motif length
# k <- n - l + 1
# 
# NextVertex(s, i, t, k)


################################################################################
################# TASK 5 - Next leaf after a skip of a subtree #################
ByPass <- function(s, i, t, k){
  for (j in i:1){
    if (s[j] < k){
      s[j] <- s[j] + 1
      return(c(s,j))
    }
  }
  return(c(s,0))
}

# sequences <- readDNAStringSet('seq_motif.fasta', format = 'fasta')
# s <- c(1,1,1,1)
# i <- 2
# t <- length(sequences) # number of sequences
# n <- width(sequences)[1] # length of each sequence
# l <- 5 # motif length
# k <- n - l + 1
# ByPass(s,i,t,k)


####################################################################################
# TASK 6 - searching for array of starting positions with best score for consensus #
BBMotifSearch <- function(DNA, t, n, l){
  s <- rep(1,t)
  bestScore <- 0
  i <- 1
  while (i > 0){
    if (i < t){
      optimisticScore <- (ScoreModified(s, i, DNA, l) + (t - 1) * l)
      if (optimisticScore < bestScore){
        temp <- ByPass(s, i, t, n - l + 1)
      }
      else{
        temp <- NextVertex(s, i, t, n - l + 1)
      }
    }
    else{
      if (ScoreModified(s, t, DNA, l) > bestScore){
        bestScore <- ScoreModified(s, t, DNA, l)
        bestMotif <- c(s)
      }
      temp <- NextVertex(s, i, t, n - l + 1)
    }
    s <- temp[1:t]
    i <- temp[t+1]
  }
  return(bestMotif)
}

DNA <- readDNAStringSet('seq_motif.fasta', format = 'fasta')
t <- length(DNA) # number of sequences
n <- width(DNA)[1] # length of each sequence
l <- 5 # motif length
BBMotifSearch(DNA, t, n, l)

