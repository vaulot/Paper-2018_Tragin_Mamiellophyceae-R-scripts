###--      Bioinformatics treatment for metabarcoding datasets under Mothur          --###
# Novel diversity within marine Mamiellophyceae (Chlorophyta) unveiled by metabarcoding
# Tragin et Vaulot (2018)
###

#-- Input files --
# 1. The "database.fasta" contains the reference sequences used to assign the environmental sequences.
#    Each sequence from the "database.fasta" is associated to a taxonomic assignation in the "database.taxo" file.
#
# 2. The "silva.seed_v119.Euk.V2.align" is used as reference alignment to align the environmental sequences.
#
# 3. The "dataset.fasta" contains the raw environmental sequences.
#    Each sequences from the "dataset.fasta" is associated with the station, where it was sampled in the "dataset.group" file.
#--------

# Set up the path to imput files and the directory for outputs

set.dir(input=/projet/sbr/osd/pathtoinputfiles, output=/projet/sbr/osd/pathtooutputfiles) 


# Clean databases

# Keep a single copy of every sequences and create a "count_table" file
unique.seqs(fasta=database.fasta)

# Create the a list of the sequences name "accnos" file
list.seqs(fasta=database.unique.fasta)

# Extract the taxonomy of uniques reference sequences
get.seqs(accnos=database.unique.accnos,taxonomy=database.taxo)


# ref alignment SILVA - Optional step
summary.seqs(fasta=/projet/sbr/osd/database/silva.seed_v119.Euk.V2.align)


# dataset
unique.seqs(fasta=dataset.fasta)
summary.seqs(fasta=dataset.unique.fasta)

# Create a "count_table" for the dataset.
count.seqs(name=dataset.names, group=dataset.groups, processors=1)

# Delete singletons
split.abund(count=dataset.count_table, fasta=dataset.unique.fasta, cutoff=1, accnos=T)

# Filter sequences on length and presence of ambiguities     
screen.seqs(fasta=dataset.unique.abund.fasta, minlength=300, maxambig=0, count=dataset.abund.count_table)


# Align to SILVA and filter

# Align the dataset to the reference alignment
align.seqs(fasta=dataset.unique.abund.good.fasta, reference=/projet/sbr/osd/database/silva.seed_v119.Euk.pcr.V2.align, flip=T, processors=1)

# Delete positions with only gaps in the alignment
filter.seqs(fasta=dataset.unique.abund.good.align)
unique.seqs(fasta=dataset.unique.abund.good.filter.fasta, count=dataset.abund.good.count_table)

# Check and remove chimeras checking from the dataset
chimera.uchime(fasta=dataset.unique.abund.good.filter.unique.precluster.fasta, count=dataset.unique.abund.good.filter.unique.precluster.count_table, processors=1)
remove.seqs(fasta=dataset.unique.abund.good.filter.unique.precluster.fasta, accnos=dataset.unique.abund.good.filter.unique.precluster.uchime.accnos, count=dataset.unique.abund.good.filter.unique.precluster.count_table)

# Remove singletons
split.abund(count=dataset.unique.abund.good.filter.unique.precluster.pick.count_table, fasta=dataset.unique.abund.good.filter.unique.precluster.pick.fasta, cutoff=1, accnos=T)

# Assign taxonomy using PR2 as reference   
classify.seqs(fasta=dataset.unique.abund.good.filter.unique.precluster.pick.abund.fasta, count=dataset.unique.abund.good.filter.unique.precluster.pick.abund.count_table, template=database.unique.fasta, taxonomy=database.pick.taxo, processors=1, probs=T) 

#Summarize taxonomy       
summary.tax(taxonomy=dataset.unique.abund.good.filter.unique.precluster.pick.abund.pick.wang.taxonomy, count=dataset.unique.abund.good.filter.unique.precluster.pick.abund.count_table)

# Export sequence from a given lineage (example for Micromonas sequences)
get.lineage(taxon=Eukaryota;Archaeplastida;Chlorophyta;Mamiellophyceae;Mamiellales;Mamiellaceae;Micromonas; ,taxonomy=dataset.unique.abund.good.filter.unique.precluster.pick.abund.pick.wang.taxonomy, fasta=dataset.unique.abund.good.filter.unique.precluster.pick.abund.fasta, count=dataset.unique.abund.good.filter.unique.precluster.pick.abund.unique.good.count_table)

