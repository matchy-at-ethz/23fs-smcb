##########################################################################
###########             Project 2:  PROFILE HMM                 ##########
##########################################################################

#
# Parse an alignment from a file
#
parseAlignment <- function(alignmentFile) {
  alignment <- scan(file=alignmentFile, what=character(0))
  # Split rows and convert to matrix
  alignment.mat <- matrix(nrow=length(alignment), ncol=nchar(alignment[[1]]))
  for(i in 1:length(alignment)) {
    alignment.mat[i,] <- unlist(strsplit(alignment[i],""))
  }
  return(alignment.mat)
}


#
# Learn HMM from data
#
learnHMM <- function(alignment, alphabet=NULL, counts.only=FALSE) {
  # Lists with all amino acids as symbols and the state transitions
  transitions <- c("MM","MD","MI","IM","ID","II","DM","DD","DI")
  if (is.null(alphabet)) {
    alphabet <- c("A","C","D","E","F","G","H","I","K","L","M",
                  "N","P","Q","R","S","T","V","W","Y")
  }
  alRows <- nrow(alignment)
  alCols <- ncol(alignment)
  L <- 0	# Length of the model
  stateNum <- array(0, dim=alCols) # Keeping track of state number
  assignedPos <- array(FALSE, dim=alCols)
  # Getting the model states and length
  for(i in 1:alCols) {
    hyphenfreq <- sum(alignment[,i] == '-') # Hyphen frequencies
    if(hyphenfreq < alRows / 2){
      L = L +1
      assignedPos[i] <- TRUE
    }
    stateNum[i] <- L
  }
  # Initialise matrices for: (index shift because of begin state)
  # match emissions (mE)
  # insertion emssions (iE)
  # transitions (T)
  mE <- matrix(0, nrow=length(alphabet), ncol=L+1, dimnames=list(alphabet, paste(0:L, sep="")))
  iE <- matrix(0, nrow=length(alphabet), ncol=L+1, dimnames=list(alphabet, NULL))
  T  <- matrix(0, nrow=length(transitions), ncol=L+1, dimnames=list(transitions, NULL))

  mE[, 1] <- NA

  for(i in 1:alRows) {
    prev <- 'M' 		# Assume starting with a match
    prevStateNum <- 0 	# Assume a begin state
    for(j in 1:alCols) {
      aa <- alignment[i, j]
      if(assignedPos[j] == TRUE) {
        if(aa == "-")
          newState <- 'D'
        else {
          newState <- 'M'
          mE[aa, stateNum[j] + 1] <- mE[aa, stateNum[j] + 1] + 1
        }
      } else {
        if(aa == "-")
          next
        else {
          iE[aa,stateNum[j] + 1] <- iE[aa, stateNum[j] + 1] + 1
          newState <- 'I'
        }
      }
      # Update transition matrix
      # j -> j+1 save at entry j
      transition <- paste(prev, newState, sep="")
      T[transition, prevStateNum +1] <- T[transition, prevStateNum +1] +1
      prevStateNum <- stateNum[j]
      prev <- newState
    }

    # Last state is always match
    transition<- paste(prev, 'M', sep="")
    T[transition, prevStateNum +1] <- T[transition, prevStateNum +1] +1
  }

  # If !counts.only, add a pseudocount of 1
  if (!counts.only) {
    # Transition matrix
    T <- T+1
    for(i in 1:ncol(T)) {
      T[1:3,i] <- T[1:3,i] / sum(T[1:3,i])
      T[4:6,i] <- T[4:6,i] / sum(T[4:6,i])
      T[7:9,i] <- T[7:9,i] / sum(T[7:9,i])
    }
    # Match emission matrix
    mE <- mE + 1
    mE <- apply(mE, 2, function(x) x / sum(x))
    # Insertion emission matrix
    iE <- iE + 1
    iE <- apply(iE, 2, function(x) x / sum(x))
  }
  return(list(T=T, mE=mE, iE=iE, L=L, alphabet=alphabet, transitions=transitions,
              stateNum=stateNum))
}

#
# parse file with one protein per line
#
parseProteins <-function(proteinsFile) {
  proteins <- scan(file= proteinsFile, what = character(0))
  proteinList <- as.list(strsplit(proteins,""))
  return(proteinList)
}


#
# The Forward algorithm
#
forward <- function(HMM, seq) {
  L <- HMM$L
  mE <- HMM$mE
  iE <- HMM$iE
  T <- HMM$T
  len <- length(seq)

  # Random model: same emission probability for all letters in the alphabet
  qr <- 1/nrow(mE)
  # Initialize score matrices for match, insert and delete
  Fm<- Fi<- Fd<- matrix(-Inf, nrow=len+1, ncol=L+1)
  #First column for the begin state (M0) (indices get shifted by 1)
  Fm[1,1]<-0

  #We then allow transitions to I0 and D1 (start with insertion or deletion)
  for(i in 2:len+1){
    Fi[i,1] <- log(iE[seq[i-1],1]/qr) +
      log(T['MI',1]*exp(Fm[i-1,1]) +
            T['II',1]*exp(Fi[i-1,1]))
  }

  # i is the pos in the sequence, j is the profile state
  for(i in 2:(len+1)) {
    for(j in 2:(L+1)) {
      Fm[i,j] <- log(mE[seq[i-1],j]/qr) +
        log(T['MM',j-1]*exp(Fm[i-1,j-1]) +
              T['IM',j-1]*exp(Fi[i-1,j-1]) +
              T['DM',j-1]*exp(Fd[i-1,j-1]))
      Fi[i,j] <- log(iE[seq[i-1],j]/qr) +
        log(T['MI',j]*exp(Fm[i-1,j]) +
              T['II',j]*exp(Fi[i-1,j]) +
              T['DI',j]*exp(Fd[i-1,j]))
      Fd[i,j] <- log( T['MD',j-1]*exp(Fm[i,j-1]) +
                        T['ID',j-1]*exp(Fi[i,j-1]) +
                        T['DD',j-1]*exp(Fd[i,j-1]))
    }
  }

  # Termination step: Fe = Fm[len,L+2] = log(P(x|M)/P(x|R))
  Fe <- log(T['MM',L+1]*exp(Fm[len+1,L+1]) +
              T['IM',L+1]*exp(Fi[len+1,L+1]) +
              T['DM',L+1]*exp(Fd[len+1,L+1]))
  return(Fe)
}
