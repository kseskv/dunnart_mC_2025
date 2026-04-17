.libPaths("/g/data/vn68/aa3618/tools/R/4.4.2/libs/")

library(data.table)
library(GenomicRanges)
library(plyr)
library(Biostrings)

setwd("/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/bismark/meth_bedGraph/")
bedgraphs <- list.files(pattern = "\\.bismark\\.cov$")
ids <- gsub("\\.bismark\\.cov$", "", bedgraphs)

####--- import chrnames for dunnart list
dunnartchrs <- read.table("/g/data/vn68/aa3618/genomes/Sminthopsis_crassicaudata/T2T_ONT_genome/lambda/dunnart_male_T2T_21082025.lambda.sizes", sep="\t")$V1

samples <- list()

for (i in bedgraphs){
  tab <- fread(i)
  tab <-  GRanges(tab[[1]], IRanges(tab[[2]], width=1), C=tab[[5]], cov=tab[[5]]+tab[[6]])
  samples <- c(samples, tab)
  print(i)
}

samples <- lapply(samples, function (x){seqlevels(x) <- dunnartchrs; x})

####--- import reference
dunnartgenome <- readDNAStringSet("/g/data/vn68/aa3618/genomes/Sminthopsis_crassicaudata/T2T_ONT_genome/lambda/dunnart_male_T2T_21082025.lambda.fasta")
chrs <- names(dunnartgenome)
matches <- vmatchPattern("CG", dunnartgenome)
matches <- matches[unlist(lapply(matches, length)) > 0]
chrs <- chrs[chrs %in% names(matches)]
CpGs <- GRanges(unlist(sapply(1:length(chrs), function (i) paste0(chrs[i], ":", matches[[i]]))))

for(i in 1:length(samples)){
  m <- findOverlaps(CpGs, samples[[i]])
  c <- rep(0, length(CpGs))
  c[queryHits(m)] <- samples[[i]]$C[subjectHits(m)]
  elementMetadata(CpGs)[paste0(ids[i], ".C")] <- c
  cov <- rep(0, length(CpGs))
  cov[queryHits(m)] <- samples[[i]]$cov[subjectHits(m)]
  elementMetadata(CpGs)[paste0(ids[i], ".cov")] <- cov
  message(i)
}

colnames(mcols(CpGs)) <- gsub("Dunart-d9-", "", colnames(mcols(CpGs)))

bigtable <- CpGs
meta <- data.frame(ID=gsub("Dunart-d9-", "", ids))
meta$C <- paste0(meta$ID, ".C")
meta$cov <- paste0(meta$ID, ".cov")
meth.ratios <- as.matrix(values(bigtable[,meta$C]))/as.matrix(values(bigtable[,meta$cov]))
colnames(meth.ratios) <- meta$ID
rownames(meth.ratios) <- paste(CpGs)
meth.ratios <- meth.ratios[apply(meth.ratios, 1, function (x) any(!is.na(x))),]

save(bigtable, meth.ratios, file="/scratch/vn68/aa3618/dunnart_reprogramming/scBSseq/Sminthopsis_crassicaudata_E9_scBSseq_250828/bismark/meth_bedGraph/bigtable/dunnart_E9_scBSseq_bigtable.RData")
