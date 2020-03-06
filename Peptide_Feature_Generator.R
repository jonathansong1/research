library("Peptides", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
# setwd("/Users/jonathansong/Research_Folder")

Mutant_peptide_name = 'Mutant_peptide'
f = '/Users/jonathansong/Research_Folder/clean_data.tsv'
data <- read.delim(f)
data$Peptide.length <- lengthpep(seq=data$Mutant_peptide)
data$Molecular.weight <- mw(seq=data$Mutant_peptide)

# what ph and pKscale should be used?
data$Net.charge <-  charge(seq = data$Mutant_peptide, pH = 7, pKscale = "EMBOSS")
data$Isoelectric.point <- pI(seq = data$Mutant_peptide, pKscale = "EMBOSS")

# what scale should be used
data$Hydrophobicity <- hydrophobicity(seq = data$Mutant_peptide, scale = "Eisenberg")

# membrane position
# data$Membrane.position <- as.data.frame(membpos(seq = data$Mutant_peptide_name, angle = 100))$MembPos

data$Aliphatic.index <- aIndex(seq = data$Mutant_peptide)
data$Instability.index <- instaIndex(seq = data$Mutant_peptide)
data$Boman.index <- boman(seq = data$Mutant_peptide)
write.table(data, file = '/Users/jonathansong/Research_Folder/clean_data.tsv', sep="\t")

