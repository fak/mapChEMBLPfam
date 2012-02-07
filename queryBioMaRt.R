library(biomaRt)


ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ensemblGenes <-  getBM(attributes = c("uniprot_swissprot_accession"),
filters = c('biotype','with_ox_uniprotswissprot'), 
value = list(biotype = 'protein_coding', with_uniprotswissprot = TRUE),
mart = ensembl)

seqs <- sapply(ensemblGenes, function(z) getSequence(id = z, 
	type = "uniprot_swissprot_accession", seqType = "peptide", 
	mart = ensembl))

seqFrame <- as.data.frame(cbind(seqs[[1]], seqs[[2]]), stringsAsFactors=FALSE)
write.table(seqFrame, file = 'data/proteinCoding.tab', quote=F, sep = '\t', 
col.names = TRUE, row.names = F)



