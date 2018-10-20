### Biofinormatics treatment for metabarcoding datasets under Mothur  ###
#   Tragin et Vaulot (2018)
####


set.dir(input=/projet/sbr/osd/pathtoinputfiles, output=/projet/sbr/osd/pathtooutputfiles) 


# database cleaning

unique.seqs(fasta=database.fasta)
list.seqs(fasta=database.unique.fasta)
get.seqs(accnos=database.unique.accnos,taxonomy=database.taxo)


# ref alignment SILVA
summary.seqs(fasta=/projet/sbr/osd/database/silva.seed_v119.Euk.V2.align)


# dataset
unique.seqs(fasta=dataset.fasta)
summary.seqs(fasta=dataset.unique.fasta)

count.seqs(name=dataset.names, group=dataset.groups, processors=1)
split.abund(count=dataset.count_table, fasta=dataset.unique.fasta, cutoff=1, accnos=T)  
screen.seqs(fasta=datasetV4.unique.abund.fasta, minlength=300, maxambig=0, count=dataset.abund.count_table)             


# align to SILVA and filter
align.seqs(fasta=dataset.unique.abund.good.fasta, reference=/projet/sbr/osd/database/silva.seed_v119.Euk.pcr.V2.align, flip=T, processors=1)
filter.seqs(fasta=dataset.unique.abund.good.align)
unique.seqs(fasta=dataset.unique.abund.good.filter.fasta, count=dataset.abund.good.count_table)

# Chimeras checking

chimera.uchime(fasta=dataset.unique.abund.good.filter.unique.fasta, count=dataset.unique.abund.good.filter.unique.count_table, processors=1)
remove.seqs(fasta=dataset.unique.abund.good.filter.unique.fasta, accnos=dataset.unique.abund.good.filter.unique.uchime.accnos, count=dataset.unique.abund.good.filter.unique.count_table)



# classification with PR²
split.abund(count=dataset.unique.abund.good.filter.unique.pick.count_table, fasta=dataset.unique.abund.good.filter.unique.pick.fasta, cutoff=1, accnos=T)
classify.seqs(fasta=dataset.unique.abund.good.filter.unique.pick.abund.fasta, count=dataset.unique.abund.good.filter.unique.pick.abund.count_table, template=database.unique.fasta, taxonomy=database.pick.taxo, processors=1, probs=T) 

summary.tax(taxonomy=dataset.unique.abund.good.filter.unique.pick.abund.pick.wang.taxonomy, count=dataset.unique.abund.good.filter.unique.pick.abund.count_table)


# Sequence Extraction (example for Micromonas sequences)
get.lineage(taxon=Eukaryota;Archaeplastida;Chlorophyta;Mamiellophyceae;Mamiellales;Mamiellaceae;Micromonas; ,taxonomy=dataset.unique.abund.good.filter.unique.pick.abund.pick.wang.taxonomy, fasta=dataset.unique.abund.good.filter.unique.pick.abund.fasta, count=dataset.unique.abund.good.filter.unique.pick.abund.unique.good.count_table)

