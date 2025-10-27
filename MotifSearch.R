library(Biostrings)
library(seqinr)

rm(list = ls())
setwd('C:/Users/hanic/Desktop/PRG/exercise_07')

## TASK 1 - Counting consensus score
Score <- function(start_indexes, sequences, motif_length){
#   seq_matrix <- as.matrix(seqs)
#   seq_consensus <- consensus(seq_matrix)
#   for (i in 1:ncol(seq_matrix)){
#     column <- seq_matrix[ ,1]
#   }
#   return(score)
# }
  
  num_seqs <- length(sequences)
  score <- 0
  
  for (pos in 0:(motif_length - 1)){
    column_bases <- character(num_seqs)
    
    for (s in 1:num_seqs){
      start_pos <- start_indexes[s]
      base <- subseq(sequences[[s]], start = start_pos + pos, width = 1)
      column_bases[s] <- as.character(base)
    }
    freq_table <- table(column_bases)
    score <- score + max(freq_table)
  }
  return(score)
}  
  
  
#seqs <- readDNAStringSet('seq_score.fasta', format = 'fasta')
#idxs <- c(1,1,1,1,1)
#length <- 5
#Score(idxs, seqs, length)


## TASK 2
NextLeaf <- function(s, t, k){
  for (i in t:1){
    if (s[i] < k){
      s[i] <- s[i] + k
      return(s)
    }
    s[i] <- 1
  }
  return(s)
}


s = c(1,1,1,1,1) # array of starting indexes
#t = 5   # number of sequences
#n = 70  # length of sequences  
#l = 6   # motif_length
#k = n - l + 1
#NextLeaf(s, t, k)


## TASK 3
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

sequences <- readDNAStringSet('seq_score.fasta', format = 'fasta')
DNA <- readDNAStringSet('seq_motif.fasta', format = 'fasta')
motif_length <- 5
k <- width(sequences)[1] + 1 - motif_length
s = c(1,1,1,1)
t <- length(DNA) # number of sequences
n <- width(DNA)[1] # length of each sequence
l <- 5 # motif length

BFMotifSearch(DNA, t, n, l)


## TASK 4
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


sequences <- readDNAStringSet('seq_motif.fasta', format = 'fasta')
s <- c(1,1,1,1)
i <- 2
t <- length(sequences) # number of sequences
n <- width(sequences)[1] # length of each sequence
l <- 5 # motif length
k <- n - l + 1

NextVertex(s, i, t, k)

## TASK 5
ByPass <- function(s, i, t, k){
  for (j in i:1){
    if (s[j] < k){
      s[j] <- s[j] + 1
      return(c(s,j))
    }
  }
  return(c(s,0))
}

sequences <- readDNAStringSet('seq_motif.fasta', format = 'fasta')
s <- c(1,1,1,1)
i <- 2
t <- length(sequences) # number of sequences
n <- width(sequences)[1] # length of each sequence
l <- 5 # motif length
k <- n - l + 1

ByPass(s,i,t,k)


## TASK 6
BBMotifSearch <- function(DNA, t, n, l){
  s <- c(1, 1, 1, 1)
  bestScore <- 0
  i <- i
  while (i > 0){
    if (i < t){
      optimisticScore <- (Score(s, i, DNA, l) + (t - 1) * l)
      if (optimisticScore < bestScore){
        c(s, i) <- ByPass(s, i, t, n - l + 1)
      }
      else{
        c(s, i) <- NextVertex(s, i, t, n - l + 1)
      }
    }
    else{
      if (Score(s, t, DNA, l) > bestScore){
        bestScore <- Score(s, t, DNA, l)
        bestMotif <- c(s)
      }
      c(s,i) <- NextVertex(s, i, t, n - l + 1)
    }
  }
  return(bestMotif)
}

DNA <- readDNAStringSet('seq_motif.fasta', format = 'fasta')
t <- length(sequences) # number of sequences
n <- width(sequences)[1] # length of each sequence
l <- 5 # motif length
BBMotifSearch(DNA, t, n, l)

