

library(biomaRt)


ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ensemblGenes <-  getBM(attributes = c("uniprot_swissprot_accession"),
filters = c('biotype','with_ox_uniprotswissprot'), 
value = list(biotype = 'protein_coding', with_uniprotswissprot = TRUE),
mart = ensembl)

write.table(ensemblGenes, file = 'data/proteinCoding.tab', quote=F, sep = '\t', 
col.names = TRUE, row.names = F)



