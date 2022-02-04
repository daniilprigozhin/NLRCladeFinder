#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: Minimal_Chimera_Entropy_Script.R
##
## Purpose of script: Output chimera entropy attribute file for a protein in a given alignment as a command line script.
##
## Author: Daniil Prigozhin
##
## Date Created: 2021-02-23
##
## Copyright (c) Daniil Prigozhin, 2021
## Email: daniilprigozhin@lbl.gov
##
## ---------------------------
##
## Notes: script will attempt to install several packages if not already available
##   
##
## ---------------------------


#install.packages("optparse")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("msa")

## load packages --------------------
package_list<-c("optparse","entropy","dplyr","msa","tidyverse")

load_pack <- function(x){
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      install.packages( i , dependencies = TRUE )
      require( i , character.only = TRUE , quietly = T)
    }
  }
}
load_pack(package_list)
#library(odseq)
cat("=======================================\nLoaded packages\n=======================================\n")


## get input file and options--------------------
option_list = list(
  make_option( c("-f", "--file") , type = "character" , default=NULL, 
              help="dataset file name", metavar="character"),
  make_option( c("-n", "--name") , type = "character" , default=NULL, 
              help="protein name, has to be in the alignment name space", metavar="character"),
  make_option(c("-c", "--common"), type="character", default=NULL, 
              help="protein common name (optional)", metavar="character"),
  make_option(c("-d", "--directory"), type="character", default="Test", 
              help="directory name (optional)", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("Provide protein alignment with -f option", call.=FALSE)
}
if (is.null(opt$name)){
  print_help(opt_parser)
  stop("Provide reference protein name with -n option", call.=FALSE)
}
gene <- opt$n
CN <- opt$c
file <- opt$f
maa <- readAAMultipleAlignment(file)

### Check that the protein name matches one alignment key ------------------
lm <- length(grep(pattern = gene, x=rownames(maa)))
if (lm == 1) {
  cat ("Found Reference Sequence\n")}else 
    if (lm ==0) {
      stop("Protein name not found in alignment", call.=FALSE)}else
        if (lm > 1) {
          stop("More than one protein matched the name provided", call.=FALSE)}else{stop("Error at protein name", call.=FALSE)}
            

## Masking columns by reference gene ----------------

Alph_21 <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","-")
Alph_20 <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")

RefGene <- gene
RefSeq <- as(maa, "AAStringSet")[grep(pattern = gene, x=rownames(maa))]
GapMask<- NULL
for (i in 1:width(RefSeq)){
  c<-as.vector(RefSeq[[1]][i]) %in% c("-")
  GapMask<-append(GapMask,c,length(GapMask))
}
colmask(maa) <- IRanges(GapMask)

## Masking reference genes ---------------------------
# RefIDs <- grep(pattern = "Athaliana", x=names(unmasked(maa)))
# rowmask(maa) <- IRanges(start = RefIDs, end = RefIDs)

#Retrieving the non-masked subset ------------------
RefAli <- as(maa, "AAStringSet")
RefLen <- width(RefAli[1])

## Calculating Consensus Matrix -------------------
Tidy_CM<-as_tibble(t(consensusMatrix(RefAli, baseOnly = T)))

## Compensating for consensus matrix not keeping full alphabet in output
for (a in setdiff(Alph_21,colnames(Tidy_CM))){
  vec <- as_tibble(0*(1:nrow(Tidy_CM)))
  colnames(vec) <- paste(a)
  Tidy_CM <- as_tibble(cbind(Tidy_CM,vec))
} 
##Selecting relevant columns
Tidy_CM_NoGaps <- select(Tidy_CM,all_of(Alph_20))

##Entropy Calculation Ignoring Gaps ----------------------
entNG <- apply(Tidy_CM_NoGaps, 1, entropy,unit="log2") %>% as_tibble()
colnames(entNG)<-paste0("EntropyNoGaps_",gene)

## Prepare output directory ------------------------------
if (!dir.exists(opt$d)){dir.create(opt$d)}
OutputDirectory <- paste0(opt$d,"/",gene,"_",CN,"/")
if (!dir.exists(OutputDirectory)){dir.create(OutputDirectory)}

#Output entropy results to file --------------------------
sink(file = paste0(OutputDirectory,gene,"_",CN,".ChimeraEntropy.txt"),append = F)
cat("attribute: shannonEntropy\n")
cat("match mode: 1-to-1\n")
cat("recipient: residues\n")
for (ii in seq_along(entNG[[1]])){
  cat("\t")
  cat(paste0(":",ii))
  cat("\t")
  cat(sprintf("%.5f", entNG[[1]][ii]))
  cat("\n")
}
sink()   

Ent <- as_tibble(cbind(1:nrow(entNG),entNG))
colnames(Ent)<-c("Position","Entropy")

ggplot(Ent, aes(x = Position))+
  geom_line(aes(y = Entropy), color = "red")+
  ylim(0,3)+
  ggtitle(paste(gene)) +
  xlab("Position") + 
  ylab("Shannon Entropy")+
  theme_classic()
ggsave(paste0(OutputDirectory,gene,"_","Entropy_MaskedNG",".pdf"))


## Report result location --------------------------------
cat("=======================================\nFinished analysis of alignment: ",opt$file,
    "\nResults in:",OutputDirectory,"\n=======================================\n")
