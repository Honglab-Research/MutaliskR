# Functions to prepare indel data for signature identification.
#
# Author: Andy Jinseok Lee


#' @title Preprocess indels
#'
#' @description This function preprocesses a data.frame of indels.
#'
#' @param df A data.frame with the following columns: Chr, Pos, Ref, Alt.
#' @param bsg BSgenome object.
#' @param padding.len Number of bases to fetch upstream and downstream of each indel.
#' @param analyze.variants A boolean value that indicates whether variant-level signature analysis should be performed.
#'
#' @return A data.frame with the following columns:
#' \item{Chr}{Chromosome name.}
#' \item{Pos}{Genomic position.}
#' \item{Ref}{Reference nucleotide sequence.}
#' \item{Alt}{Alternate nucleotide sequence.}
#' \item{Type}{Either 'INS' (insertion) or 'DEL' (deletion).}
#' \item{Len}{Length (number of bases) of this indel.}
#' \item{Mutated_Seq}{Nucleotide sequence of this indel.}
#' \item{Upstream_Seq}{Upstream nucleotide sequence of this indel.}
#' \item{Upstream_Start}{Upstream nucleotide sequence start position.}
#' \item{Upstream_End}{Upstream nucleotide sequence end position.}
#' \item{Downstream_Seq}{Downstream nucleotide sequence of this indel.}
#' \item{Downstream_Start}{Downstream nucleotide sequence start position..}
#' \item{Downstream_End}{Downstream nucleotide sequence end position.}
#'
#' @export
#' @import dplyr
#' @import BSgenome
PreprocessIndels <- function(df, bsg, padding.len, analyze.variants) {
  # Filter out variants that are not indels based on lengths of REF and ALT
  df <- df[nchar(df$Ref) != nchar(df$Alt),]

  if (nrow(df) == 0) {
    df.temp <- as.data.frame(matrix(, ncol = 11, nrow = 0))
    colnames(df.temp) <- c("Chr", "Pos", "Ref", "Alt",
                           "Type", "Len", "Mutated_Seq",
                           "Upstream_Start", "Upstream_End",
                           "Downstream_Start", "Downstream_End")
    return(df.temp)
  }

  PreprocessIndelsHelper <- function(x) {
    curr.chr <- as.character(x[["Chr"]])
    curr.pos <- as.integer(x[["Pos"]])
    curr.ref <- toupper(as.character(x[["Ref"]]))
    curr.alt <- toupper(as.character(x[["Alt"]]))

    # Determine whether the current variant is a small insertion or a small deletion
    if (nchar(curr.ref) < nchar(curr.alt)) {
      curr.indel.type <- "INS"
      curr.mutated.seq <- substring(text = curr.alt, first = 2, last = nchar(curr.alt))
    } else {
      curr.indel.type <- "DEL"
      curr.mutated.seq <- substring(text = curr.ref, first = 2, last = nchar(curr.ref))
    }
    curr.indel.len <- nchar(curr.mutated.seq)

    # Get 'padding.len' number of padding sequences
    if (curr.indel.type == "INS") {
      upstream.start <- curr.pos - padding.len
      upstream.end <- curr.pos - 1
      downstream.start <- curr.pos + 1
      downstream.end <- curr.pos + padding.len
    }
    if (curr.indel.type == "DEL") {
      upstream.start <- curr.pos - padding.len
      upstream.end <- curr.pos - 1
      downstream.start <- (curr.pos + nchar(curr.ref) - 1) + 1
      downstream.end <- (curr.pos + nchar(curr.ref) - 1) + padding.len
    }
    df.temp <- data.frame(Chr = curr.chr,
                          Pos = curr.pos,
                          Ref = curr.ref,
                          Alt = curr.alt,
                          Type = curr.indel.type, # 'INS' or 'DEL'
                          Len = curr.indel.len, # indel length
                          Mutated_Seq = curr.mutated.seq, # the inserted or deleted sequence
                          Upstream_Start = upstream.start,
                          Upstream_End = upstream.end,
                          Downstream_Start = downstream.start,
                          Downstream_End = downstream.end,
                          stringsAsFactors = FALSE)

    if (analyze.variants == TRUE) {
      df.temp$Variant_Gene <- as.character(x[["Variant_Gene"]])
      df.temp$Variant_Group <- as.character(x[["Variant_Group"]])
    }

    return(df.temp)
  }

  df.indels <- apply(df, MARGIN = 1, FUN = PreprocessIndelsHelper)
  df.indels <- bind_rows(df.indels)

  # Get upstream and downstream sequences
  df.indels$Upstream_Seq <- as.character(BSgenome::getSeq(bsg,
                                                          names = df.indels$Chr,
                                                          start = df.indels$Upstream_Start,
                                                          end = df.indels$Upstream_End))
  df.indels$Downstream_Seq <- as.character(BSgenome::getSeq(bsg,
                                                            names = df.indels$Chr,
                                                            start = df.indels$Downstream_Start,
                                                            end = df.indels$Downstream_End))
  return(df.indels)
}

#' @title Returns indel repeat size
#'
#' @description This function computes and returns the repeat size of an indel variant.
#'
#' @param downstream.seq Downstream nucleotide sequence.
#' @param mutated.seq Indel nucleotide sequence.
#' @param type Indel type. Either 'INS' or 'DEL'.
#'
#' @return Repeat size (integer).
#'
#' @export
#' @import stringr
GetIndelRepeatSize <- function(downstream.seq, mutated.seq, type) {
  # Identify all occurrences of 'mutated.seq' in 'downstream.seq'
  df.locations <- as.data.frame(str_locate_all(string = downstream.seq, pattern = mutated.seq))

  # Enumerate the number of consecutive repeats
  repeat.size <- 0
  if (nrow(df.locations) == 0) { # 0 repeat found
    return(repeat.size)
  }
  if (head(df.locations, 1)$start == 1) { # first repeat is found at the start of the downstream.seq
    repeat.size <- repeat.size + 1
    for (i in 1:nrow(df.locations)) {
      if (i > 1) {
        if (as.integer(df.locations[i, 'start']) == prev.end + 1) { # consecutive repeat
          repeat.size <- repeat.size + 1
        } else {
          break
        }
      }
      prev.end <- as.integer(df.locations[i,'end'])
    }
  }

  return(repeat.size)
}

#' @title Returns microhomology size of a deletion
#'
#' @description This function computes and returns the microhomology size of a small deletion.
#'
#' @param upstream.seq Upstream nucleotide sequence.
#' @param downstream.seq Downstream nucleotide sequence.
#' @param mutated.seq Indel nucleotide sequence
#'
#' @return A data.frame with the following columns:
#' \item{Microhomology_Seq}{Microhomology nucleotide sequence.}
#' \item{Microhomology_Size}{Size (number of bases) of microhomology.}
#' \item{Microhomology_Direction}{Direction of microhomology.}
#'
#' @export
GetDeletionMicrohomologySize <- function(upstream.seq, downstream.seq, mutated.seq) {
  # Prepare some data
  upstream.seq.bases <- strsplit(x = upstream.seq, "")[[1]]
  downstream.seq.bases <- strsplit(x = downstream.seq, "")[[1]]
  mutated.seq.bases <- strsplit(x = mutated.seq, "")[[1]]

  # Get microhomology size and sequence
  microhomology.seqs <- c() # both upstream and downstream sequences might be the same
  microhomology.sizes <- c()
  microhomology.directions <- c()
  for (curr.window.size in 1:(length(mutated.seq.bases) - 1)) {
    # Forward (mutated.seq --> downstream.seq)
    curr.mutated.seq <- paste0(mutated.seq.bases[1:curr.window.size], collapse = "")
    curr.downstream.seq <- paste0(downstream.seq.bases[1:curr.window.size], collapse = "")
    if (identical(curr.mutated.seq, curr.downstream.seq)) {
      microhomology.seqs <- c(microhomology.seqs, curr.downstream.seq)
      microhomology.sizes <- c(microhomology.sizes, curr.window.size)
      microhomology.directions <- c(microhomology.directions, "downstream")
    }

    # Reverse (upstream.seq <-- mutated.seq)
    curr.mutated.seq <- paste0(mutated.seq.bases[(length(mutated.seq.bases) - curr.window.size + 1):length(mutated.seq.bases)], collapse = "")
    curr.upstream.seq <- paste0(upstream.seq.bases[(length(upstream.seq.bases) - curr.window.size + 1):length(upstream.seq.bases)], collapse = "")
    if (identical(curr.mutated.seq, curr.upstream.seq)) {
      microhomology.seqs <- c(microhomology.seqs, curr.upstream.seq)
      microhomology.sizes <- c(microhomology.sizes, curr.window.size)
      microhomology.directions <- c(microhomology.directions, "upstream")
    }
  }

  # Get the microhomology with the largest size
  df <- data.frame(Microhomology_Seq = microhomology.seqs,
                   Microhomology_Size = microhomology.sizes,
                   Microhomology_Direction = microhomology.directions,
                   stringsAsFactors = FALSE)
  if (nrow(df) > 0) {
    df <- df[df$Microhomology_Size == max(df$Microhomology_Size),]
  } else {
    df <- data.frame(Microhomology_Seq = c(),
                     Microhomology_Size = c(),
                     Microhomology_Direction = c(),
                     stringsAsFactors = FALSE)
  }
  return(df)
}

#' @title Returns PCAWG classification of an indel
#'
#' @description This function classifies an indel according to the PCAWG classification.
#'
#' @param type Indel type. Either 'INS' or 'DEL'.
#' @param mutated.seq Indel nucleotide sequence.
#' @param upstream.seq Upstream nucleotide sequence.
#' @param downstream.seq Downstream nucleotide sequence.
#'
#' @return A data.frame with the following columns:
#' \item{PCAWG_Indel_Class}{PCAWG indel mutation type (e.g. 'DEL_C_1_0').}
#' \item{Microhomology_Class}{PCAWG microhomology mutation type (e.g. 'DEL_MH_2_1').}
#' \item{Microhomology_Size}{Microhomology size.}
#' \item{Microhomology_Seq}{Microhomology sequence.}
#' \item{Microhomology_Direction}{Microhomology direction.}
#'
#' @export
#' @import bedr
#' @import dplyr
GetIndelPCAWGClassification <- function(type, mutated.seq, upstream.seq, downstream.seq) {
  # Get repeat size
  repeat.size <- GetIndelRepeatSize(downstream.seq = downstream.seq,
                                    mutated.seq = mutated.seq,
                                    type = type)

  # Reword repeat size
  if (repeat.size >= 5 && type == "DEL") {
    repeat.size <- "5+"
  }
  if (repeat.size >= 5 && type == "INS") {
    repeat.size <- "5+"
  }

  # Mutation type
  if (nchar(mutated.seq) == 1) { # 1bp DEL or INS
    # Get mutation type base
    if (mutated.seq %in% c("T", "A")) {
      base.pair <- "T"
    }
    if (mutated.seq %in% c("C", "G")) {
      base.pair <- "C"
    }
    mutation.type <- paste0(type, "_", base.pair, "_1_", repeat.size)
  } else { # 2+bp DEL or INS
    if (nchar(mutated.seq) >= 5) {
      mutation.type <- paste0(type, "_repeats_5+_", repeat.size)
    } else {
      mutation.type <- paste0(type, "_repeats_", nchar(mutated.seq), "_", repeat.size)
    }
  }

  # Microhomology
  microhomology.group <- ""
  microhomology.seq <- ""
  microhomology.direction <- ""
  microhomology.size <- 0
  if ((type == "DEL") && (nchar(mutated.seq) > 1)) {
    df.microhomology <- GetDeletionMicrohomologySize(upstream.seq = upstream.seq,
                                                     downstream.seq = downstream.seq,
                                                     mutated.seq = mutated.seq)
    if (nrow(df.microhomology) > 0) {
      microhomology.seq <- paste0(df.microhomology$Microhomology_Seq, collapse = ",")
      microhomology.direction <- paste0(df.microhomology$Microhomology_Direction, collapse = ",")
      microhomology.size <- unique(df.microhomology$Microhomology_Size)

      # Construct label
      if (microhomology.size >= 5) {
        microhomology.size.label <- "5+"
      } else {
        microhomology.size.label <- microhomology.size
      }
      if (nchar(mutated.seq) >= 5) {
        mutated.seq.label <- "5+"
      } else {
        mutated.seq.label <- nchar(mutated.seq)
      }
      microhomology.group <- paste0("DEL_MH_", mutated.seq.label, "_", microhomology.size.label)
    }
  }

  df.temp <- data.frame(PCAWG_Indel_Class = mutation.type,
                        Microhomology_Class = microhomology.group,
                        Microhomology_Size = microhomology.size,
                        Microhomology_Seq = microhomology.seq,
                        Microhomology_Direction = microhomology.direction,
                        stringsAsFactors = FALSE)

  return(df.temp)
}

#' @title Assigns PCAWG classification of an indel
#'
#' @description This function returns PCAWG classification data of an indel.
#'
#' @param x One row of data.frame the function PreprocessIndel.
#'
#' @return A data.frame with the following columns:
#' \item{PCAWG_Indel_Class}{PCAWG indel mutation type (e.g. 'DEL_C_1_1').}
#' \item{Microhomology_Class}{PCAWG microhomology mutation type (e.g. 'DEL_MH_2_1').}
#' \item{Microhomology_Size}{Microhomology size.}
#' \item{Microhomology_Seq}{Microhomology sequence.}
#' \item{Microhomology_Direction}{Microhomology direction.}
#'
#' @export
#' @import bedr
#' @import dplyr
AssignIndelPCAWGClassification <- function(x) {
  curr.type <- as.character(x[['Type']])
  curr.mutated.seq <- as.character(x[['Mutated_Seq']])
  curr.upstream.seq <- as.character(x[['Upstream_Seq']])
  curr.downstream.seq <- as.character(x[['Downstream_Seq']])
  df <- GetIndelPCAWGClassification(type = curr.type,
                                    mutated.seq = curr.mutated.seq,
                                    upstream.seq = curr.upstream.seq,
                                    downstream.seq = curr.downstream.seq)
  return(df)
}

#' @title Prepares indel data.frame
#'
#' @description
#' This function prepares an indel data.frame and
#' returns relevant data to run IdentifySignatures function.
#'
#' @param df A data.frame with the following columns: Chr, Pos, Ref, Alt.
#' @param bsg BSgenome object.
#' @param analyze.variants A boolean value that indicates whether variant-level signature analysis should be performed.
#' @param max.len Maximum length of an indel to include. Indels longer than this length will be excluded.
#' @param padding.len Number of bases to pad upstream and downstream of each indel.
#'
#' @return A list with the following items:
#' \item{df.indel}{A data.frame of the indel data.}
#' \item{df.indel.frequencies}{A data.frame of the PCAWG indel class frequencies.}
#'
#' @export
#' @import bedr
#' @import dplyr
PrepareIndelDataFrame <- function(df,
                                  bsg,
                                  analyze.variants,
                                  max.len = 25,
                                  padding.len = 80) {
  # Step 1. Preprocess indels
  df.indels <- PreprocessIndels(df = df,
                                bsg = bsg,
                                padding.len = padding.len,
                                analyze.variants = analyze.variants)

  # Step 2. Filter out indels that have lengths greater than 'max.len'
  df.indels <- df.indels[df.indels$Len <= max.len,]

  if (nrow(df.indels) == 0) {
    df.pcawg.classes <- as.data.frame(matrix(, ncol = 6, nrow = 0))
    colnames(df.pcawg.classes) <- c("PCAWG_Indel_Class",
                                    "Microhomology_Class",
                                    "Microhomology_Size",
                                    "Microhomology_Seq",
                                    "Microhomology_Direction",
                                    "Mutation_Type_Group")
    df.indels <- cbind(df.indels, df.pcawg.classes)

    df.class.frequencies <- data.frame(
      Mutation_Type_Group = INDEL.MUTATION.TYPES.GROUPS,
      Mutation_Type = INDEL.MUTATION.TYPES,
      Frequency = 0,
      stringsAsFactors = FALSE
    )

    return(list(df.indel = df.indels, df.indel.frequencies = df.class.frequencies))
  }

  # Step 3. Classify each indel according to the PCAWG indel classification
  df.mutation.types <- apply(df.indels, MARGIN = 1, FUN = AssignIndelPCAWGClassification)
  df.mutation.types <- bind_rows(df.mutation.types)

  # Step 4. Merge
  df.indels <- cbind(df.indels, df.mutation.types)

  # Step 5. Get mutation type group
  df.id.signatures <- data.frame(Mutation_Type = INDEL.MUTATION.TYPES,
                                 Mutation_Type_Group = INDEL.MUTATION.TYPES.GROUPS,
                                 stringsAsFactors = FALSE)
  GetMutationTypeGroup <- function(x) {
    curr.pcawg.indel.class <- as.character(x[['PCAWG_Indel_Class']])
    curr.mutation.type.group <- as.character(df.id.signatures[df.id.signatures$Mutation_Type == curr.pcawg.indel.class, "Mutation_Type_Group"])
    return(curr.mutation.type.group)
  }
  df.indels$Mutation_Type_Group <- apply(df.indels, MARGIN = 1, FUN = GetMutationTypeGroup)
  df.indels$Mutation_Type <- df.indels$PCAWG_Indel_Class

  # Step 6. Enumerate by PCAWG indel classes
  df.class.frequencies <- data.frame()
  for (i in 1:nrow(df.id.signatures)) {
    curr.mutation.type <- as.character(df.id.signatures[i, "Mutation_Type"])
    curr.mutation.type.group <- as.character(df.id.signatures[i, "Mutation_Type_Group"])
    if (str_detect(string = curr.mutation.type, pattern = "DEL_MH")) {
      df.temp <- data.frame(Mutation_Type_Group = curr.mutation.type.group,
                            Mutation_Type = curr.mutation.type,
                            Frequency = nrow(df.indels[df.indels$Microhomology_Class == curr.mutation.type,]),
                            stringsAsFactors = FALSE)
    } else {
      df.temp <- data.frame(Mutation_Type_Group = curr.mutation.type.group,
                            Mutation_Type = curr.mutation.type,
                            Frequency = nrow(df.indels[df.indels$PCAWG_Indel_Class == curr.mutation.type,]),
                            stringsAsFactors = FALSE)
    }
    df.class.frequencies <- rbind(df.class.frequencies, df.temp)

  }
  return(list(df.indel = df.indels, df.indel.frequencies = df.class.frequencies))
}

#' @title Prepares indel VCF file
#'
#' @description
#' This function prepares an indel VCF file and
#' returns relevant data to run IdentifySignatures function.
#'
#' @param vcf.file VCF file including path.
#' @param bsg BSgenome object.
#' @param max.len Maximum length of an indel to include. Indels longer than this length will be excluded.
#' @param padding.len Number of bases to pad upstream and downstream of each indel.
#'
#' @return A list with the following items:
#' \item{df.indel}{A data.frame of the indel data.}
#' \item{df.indel.frequencies}{A data.frame of the PCAWG indel class frequencies.}
#'
#' @export
#' @import bedr
#' @import dplyr
PrepareIndelVcfFile <- function(vcf.file, bsg, max.len = 25, padding.len = 80) {
  df.temp <- PrepareVcfFile(vcf.file = vcf.file)
  prepared.indel.data <- PrepareIndelDataFrame(df = df.temp,
                                               bsg = bsg,
                                               max.len = max.len,
                                               padding.len = padding.len,
                                               analyze.variants = FALSE)
  return(prepared.indel.data)
}
