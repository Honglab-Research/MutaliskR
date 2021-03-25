# This Rscript prepares the PCAWG (SignatureAnalyzer) pentanucleotide
# SBS mutational signature reference matrix for identification with MutaliskR.


library(dplyr)

# Step 1. Load data.
reference.file <- "<PATH>/MutaliskR/inst/extdata/01_Reference_Signatures/01_PCAWG_Signatures/SignatureAnalyzer_COMPOSITE_SBS_W1536.signature.031918.txt"
df <- read.csv(reference.file, sep = "\t")

# Step 2. Remove the last column ("feature").
column.names <- colnames(df)
column.names <- column.names[column.names != "feature"]

# Step 3. Relabel the context in the following format: NN[Ref>Alt]NN
mutation.types <- c()
for (curr.feature in df[,"feature"]) {
  curr.ref <- substr(curr.feature, start = 1, stop = 1)
  curr.alt <- substr(curr.feature, start = 3, stop = 3)
  curr.upstream <- substr(curr.feature, start = 8, stop = 9)
  curr.downstream <- substr(curr.feature, start = 11, stop = 12)
  mutation.type <- paste0(curr.upstream, "[", curr.ref, ">", curr.alt, "]", curr.downstream)
  mutation.types <- c(mutation.types, mutation.type)
}
df <- df[,column.names]
df$Mutation_Type <- mutation.types

# Step 4. Make Mutation_Type first column
df <- df %>%
  select(Mutation_Type, everything())

# Step 5. Write to file
write.csv(x = df,
          file = "<PATH>/MutaliskR/inst/extdata/01_Reference_Signatures/01_PCAWG_Signatures/SignatureAnalyzer_COMPOSITE_SBS_Pentanucleotide_signatures_031918.csv",
          row.names = FALSE, quote = FALSE)

