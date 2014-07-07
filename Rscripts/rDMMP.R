# -----------------------------------------
# rDMMP = R Drug Markov Mean Properties
# -----------------------------------------
# by Cristian R Munteanu, muntisa@gmail.com
# BiGCaT, Maastricht University
# -----------------------------------------
library(ChemmineR)
library(base)
library(expm)

start.time <- Sys.time()

#-------------------------------------------------------------------------------
# PARAMETERS
#-------------------------------------------------------------------------------
kPower <- 5                          # power for Markov chains (0 - kPower)
SFile <- "SMILES.txt"                # SMILES file
WFile <- "AtomProperties.txt"        # weigths file:  AtomE EM PKJ vdWArea AC2P (tab-separated)
#-------------------------------------------------------------------------------
# Read SMILES formulas
#-------------------------------------------------------------------------------
smiset <- read.SMIset(SFile)         # read a file with SMILES codes as object
smiles <- as.character(smiset)       # convert the object with SMILES into a text list
nSmiles <- length(smiles)            # number of SMILES formulas
#-------------------------------------------------------------------------------
# Read WEIGHTs from file
#-------------------------------------------------------------------------------
dfW=read.table(WFile,header=T)       # read weigths file
Headers <- names(dfW)                # list variable names into the dataset
wNoRows=dim(dfW)[1]                  # number of rows = atoms with properties
wNoCols=dim(dfW)[2]                  # number of cols = atomic element, properties 
#-------------------------------------------------------------------------------
# Initialize the header of RESULTS string
sResults <- c("Molecule")
for(h in 2:wNoCols){ 
  sResults <- c(sResults,sprintf("MMP_%s",Headers[h]))
}
sResults <- c(sResults,"\n")
#-------------------------------------------------------------------------------
# PROCESS EACH SMILES
# - calculate MMPs for each SMILES, each pythical-chemical property and atom type
#   averaged for all powers (0-kPower)
#-------------------------------------------------------------------------------
for(s in 1:nSmiles){                 # process each SMILES
  smi <- smiles[[s]]                 # SMILES formula
  molName <- names(smiles)[s]        # molecule label
  print(c(s,molName))
  sdf <- smiles2sdf(smi)             # convert one smiles to sdf format
  BM <- conMA(sdf,exclude=c("H"))    # bond matrix (complex list!)
  
  # Connectivity matrix CM
  CM <- BM[[1]]                      # get only the matrix  
  CM[CM > 0] <- 1                    # convert bond matrix into connectivity matrix/adjacency matrix CM

  # Degrees
  deg <- rowSums(CM)                 # atom degree (no of chemical bonds)
  atomNames <- (rownames(CM))        # atom labels from the bond table (atom_index)
  nAtoms <- length(atomNames)        # number of atoms
  
  # Get list with atoms and positions
  Atoms <- list()                    # inicialize the list of atoms 
  AtIndexes <- list()                # inicialize the list of atom indixes
  for(a in 1:nAtoms){                # process each atom in bond table
    Atom <- atomNames[a]                      # pick one atom
    Elem_Ind <- strsplit(Atom,'_')            # split atom labels (atom_index)
    AtElem <- Elem_Ind[[1]][1]                # get atomic element
    AtIndex <- Elem_Ind[[1]][2]               # get index of atom element
    Atoms[a] <- c(AtElem)                     # add atom element to a list
    AtIndexes[a] <- as.numeric(c(AtIndex))    # add index of atom elements to a list
  }
  
  # Weights data frame (for all atom properties)
  # -----------------------------------------------------------------------
  # Atoms data frame
  dfAtoms <- data.frame(Pos=as.numeric(t(t(AtIndexes))),AtomE=t(t(Atoms)))
  # Weights data frame
  dfAtomsW <- merge(dfW,dfAtoms,by="AtomE") # merge 2 data frame using values of AtomE
  dfAtomsW <- dfAtomsW[order(dfAtomsW$Pos),1:wNoCols] # order data frame and remove Pos
  rownames(dfAtomsW) <- seq(1:nAtoms)
  # NEED CORRECTION FOR ATOMS that are not in the properties file
  
  # For each atom property
  # -----------------------------------------
  vMMP <- c()                        # final MMP descriptors for one molecule
  for(prop in 2:wNoCols){                    # for each property
    w <- t(data.matrix(dfAtomsW[prop]))[1,]    # weigths VECTOR
    W <- t(CM * w)                             # weigthed MATRIX
    p0j <- w/sum(w)                            # absolute initial probabilities vector
    
    # Probability Matrix P
    # ----------------------
    degW <- rowSums(W) # degree vector
    P <- W * ifelse(degW == 0, 0, 1/degW)      # transition probability matrix (corrected for zero division)
    
    # Average all descriptors by power (0-k)
    # ------------------------------------------------------------
    vMMPk <- c()                               # MMPs for k values
    for(k in 0:kPower){                        # power between 0 and kPower
      Pk <- (P %^% k)                            # Pk = k-powered transition probability matrix
      MMPk = t(p0j) %*% Pk %*% t(t(w))           # MMPk for all atoms for one k = vMv type product
      vMMPk <- c(vMMPk, MMPk)                    # vector with all MMPs for all ks and atom properties
    }
    vMMP <- c(vMMP,mean(vMMPk))                # average value for all k values ( +1 because of k starting with 0)
  }
  # print final MMPs for one molecule, for each property (averaged by all k values), for all types of atoms
  # Molecule_Name MMP_EM MMP_PKJ MMP_vdWA MMP_AC2P
  sResults <- c(sResults,molName,vMMP,"\n")
  # additional calculations will be added for each type of atom!!!!
}

# print final output
cat(sResults, sep = "\t")


end.time <- Sys.time()
time.taken <- start.time - end.time
time.taken

# 500 SMILES in 13.5 s (on i5,8G RAM,SSD)
