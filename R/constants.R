# Constants.
#
# Author: Andy Jinseok Lee


#' @export
SBS.PCAWG.KNOWN.ETIOLOGY.MUTATIONAL.SIGNATURES <- c("SBS1","SBS10a","SBS10b","SBS11","SBS13",
                                                    "SBS14","SBS15","SBS16","SBS17b","SBS18",
                                                    "SBS2","SBS20","SBS21","SBS22","SBS24",
                                                    "SBS26","SBS29","SBS3","SBS30","SBS31",
                                                    "SBS32","SBS35","SBS36","SBS38","SBS4",
                                                    "SBS42","SBS44","SBS5","SBS6","SBS7a",
                                                    "SBS7b","SBS7c","SBS7d","SBS84","SBS85","SBS9")

#' @export
SBS.MUTATION.TYPES <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T",
                        "A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T",
                        "A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T",
                        "A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T",
                        "A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T",
                        "A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")

#' @export
SBS.MUTATION.TYPES.GROUPS <- c("C>A","C>A","C>A","C>A","C>A","C>A","C>A","C>A","C>A","C>A","C>A","C>A","C>A","C>A","C>A","C>A",
                               "C>G","C>G","C>G","C>G","C>G","C>G","C>G","C>G","C>G","C>G","C>G","C>G","C>G","C>G","C>G","C>G",
                               "C>T","C>T","C>T","C>T","C>T","C>T","C>T","C>T","C>T","C>T","C>T","C>T","C>T","C>T","C>T","C>T",
                               "T>A","T>A","T>A","T>A","T>A","T>A","T>A","T>A","T>A","T>A","T>A","T>A","T>A","T>A","T>A","T>A",
                               "T>C","T>C","T>C","T>C","T>C","T>C","T>C","T>C","T>C","T>C","T>C","T>C","T>C","T>C","T>C","T>C",
                               "T>G","T>G","T>G","T>G","T>G","T>G","T>G","T>G","T>G","T>G","T>G","T>G","T>G","T>G","T>G","T>G")

#' @export
PLOT.SBS.X.AXIS.LABELS <- c("ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT",
                            "ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT",
                            "ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT",
                            "ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT",
                            "ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT",
                            "ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT")

#' @export
PLOT.SBS.STRIP.LABELS <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

#' @export
PLOT.SBS.STRIP.COLORS <- c("#03BCEE", "#000000", "#E32926", "#999999", "#A1CE63", "#EBC6C4")

#' @export
PLOT.SBS.STRIP.TEXT.COLORS <- c("#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#000000", "#000000")

#' @export
PLOT.SBS.LEGEND.LABELS <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")


#' @export
DBS.MUTATION.TYPES <- c("AC>CA","AC>CG","AC>CT","AC>GA","AC>GG","AC>GT","AC>TA","AC>TG","AC>TT",
                        "AT>CA","AT>CC","AT>CG","AT>GA","AT>GC","AT>TA",
                        "CC>AA","CC>AG","CC>AT","CC>GA","CC>GG","CC>GT","CC>TA","CC>TG","CC>TT",
                        "CG>AT","CG>GC","CG>GT","CG>TA","CG>TC","CG>TT",
                        "CT>AA","CT>AC","CT>AG","CT>GA","CT>GC","CT>GG","CT>TA","CT>TC","CT>TG",
                        "GC>AA","GC>AG","GC>AT","GC>CA","GC>CG","GC>TA",
                        "TA>AT","TA>CG","TA>CT","TA>GC","TA>GG","TA>GT",
                        "TC>AA","TC>AG","TC>AT","TC>CA","TC>CG","TC>CT","TC>GA","TC>GG","TC>GT",
                        "TG>AA","TG>AC","TG>AT","TG>CA","TG>CC","TG>CT","TG>GA","TG>GC","TG>GT",
                        "TT>AA","TT>AC","TT>AG","TT>CA","TT>CC","TT>CG","TT>GA","TT>GC","TT>GG")

#' @export
DBS.MUTATION.TYPES.REVERSE.COMPLEMENTED <- c("GT>TG","GT>CG","GT>AG","GT>TC","GT>CC","GT>AC","GT>TA","GT>CA","GT>AA",
                                             "AT>TG","AT>GG","AT>CG","AT>TC","AT>GC","AT>TA",
                                             "GG>TT","GG>CT","GG>AT","GG>TC","GG>CC","GG>AC","GG>TA","GG>CA","GG>AA",
                                             "CG>AT","CG>GC","CG>AC","CG>TA","CG>GA","CG>AA",
                                             "AG>TT","AG>GT","AG>CT","AG>TC","AG>GC","AG>CC","AG>TA","AG>GA","AG>CA",
                                             "GC>TT","GC>CT","GC>AT","GC>TG","GC>CG","GC>TA",
                                             "TA>AT","TA>CG","TA>AG","TA>GC","TA>CC","TA>AC",
                                             "GA>TT","GA>CT","GA>AT","GA>TG","GA>CG","GA>AG","GA>TC","GA>CC","GA>AC",
                                             "CA>TT","CA>GT","CA>AT","CA>TG","CA>GG","CA>AG","CA>TC","CA>GC","CA>AC",
                                             "AA>TT","AA>GT","AA>CT","AA>TG","AA>GG","AA>CG","AA>TC","AA>GC","AA>CC")

#' @export
DBS.MUTATION.TYPES.GROUPS <- c("AC>NN","AC>NN","AC>NN","AC>NN","AC>NN","AC>NN","AC>NN","AC>NN","AC>NN",
                               "AT>NN","AT>NN","AT>NN","AT>NN","AT>NN","AT>NN",
                               "CC>NN","CC>NN","CC>NN","CC>NN","CC>NN","CC>NN","CC>NN","CC>NN","CC>NN",
                               "CG>NN","CG>NN","CG>NN","CG>NN","CG>NN","CG>NN",
                               "CT>NN","CT>NN","CT>NN","CT>NN","CT>NN","CT>NN","CT>NN","CT>NN","CT>NN",
                               "GC>NN","GC>NN","GC>NN","GC>NN","GC>NN","GC>NN",
                               "TA>NN","TA>NN","TA>NN","TA>NN","TA>NN","TA>NN",
                               "TC>NN","TC>NN","TC>NN","TC>NN","TC>NN","TC>NN","TC>NN","TC>NN","TC>NN",
                               "TG>NN","TG>NN","TG>NN","TG>NN","TG>NN","TG>NN","TG>NN","TG>NN","TG>NN",
                               "TT>NN","TT>NN","TT>NN","TT>NN","TT>NN","TT>NN","TT>NN","TT>NN","TT>NN")

#' @export
PLOT.DBS.X.AXIS.LABELS <- c("CA","CG","CT","GA","GG","GT","TA","TG","TT","CA","CC","CG","GA","GC","TA",
                            "AA","AG","AT","GA","GG","GT","TA","TG","TT","AT","GC","GT","TA","TC","TT",
                            "AA","AC","AG","GA","GC","GG","TA","TC","TG","AA","AG","AT","CA","CG","TA",
                            "AT","CG","CT","GC","GG","GT","AA","AG","AT","CA","CG","CT","GA","GG","GT",
                            "AA","AC","AT","CA","CC","CT","GA","GC","GT","AA","AC","AG","CA","CC","CG",
                            "GA","GC","GG")

#' @export
PLOT.DBS.STRIP.LABELS <- c("AC>NN", "AT>NN", "CC>NN", "CG>NN", "CT>NN",
                           "GC>NN", "TA>NN", "TC>NN", "TG>NN", "TT>NN")

#' @export
PLOT.DBS.STRIP.COLORS <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                           "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")

#' @export
PLOT.DBS.STRIP.TEXT.COLORS <- c("#000000", "#FFFFFF", "#000000", "#FFFFFF", "#FFFFFF",
                                "#FFFFFF", "#000000", "#FFFFFF", "#000000", "#FFFFFF")

#' @export
PLOT.DBS.LEGEND.LABELS <- c("AC>NN", "AT>NN", "CC>NN", "CG>NN", "CT>NN",
                            "GC>NN", "TA>NN", "TC>NN", "TG>NN", "TT>NN")


#' @export
INDEL.MUTATION.TYPES <- c("DEL_C_1_0","DEL_C_1_1","DEL_C_1_2","DEL_C_1_3","DEL_C_1_4","DEL_C_1_5+","DEL_T_1_0","DEL_T_1_1","DEL_T_1_2","DEL_T_1_3","DEL_T_1_4","DEL_T_1_5+",
                          "INS_C_1_0","INS_C_1_1","INS_C_1_2","INS_C_1_3","INS_C_1_4","INS_C_1_5+","INS_T_1_0","INS_T_1_1","INS_T_1_2","INS_T_1_3","INS_T_1_4","INS_T_1_5+",
                          "DEL_repeats_2_0","DEL_repeats_2_1","DEL_repeats_2_2","DEL_repeats_2_3","DEL_repeats_2_4","DEL_repeats_2_5+",
                          "DEL_repeats_3_0","DEL_repeats_3_1","DEL_repeats_3_2","DEL_repeats_3_3","DEL_repeats_3_4","DEL_repeats_3_5+",
                          "DEL_repeats_4_0","DEL_repeats_4_1","DEL_repeats_4_2","DEL_repeats_4_3","DEL_repeats_4_4","DEL_repeats_4_5+",
                          "DEL_repeats_5+_0","DEL_repeats_5+_1","DEL_repeats_5+_2","DEL_repeats_5+_3","DEL_repeats_5+_4","DEL_repeats_5+_5+",
                          "INS_repeats_2_0","INS_repeats_2_1","INS_repeats_2_2","INS_repeats_2_3","INS_repeats_2_4","INS_repeats_2_5+",
                          "INS_repeats_3_0","INS_repeats_3_1","INS_repeats_3_2","INS_repeats_3_3","INS_repeats_3_4","INS_repeats_3_5+",
                          "INS_repeats_4_0","INS_repeats_4_1","INS_repeats_4_2","INS_repeats_4_3","INS_repeats_4_4","INS_repeats_4_5+",
                          "INS_repeats_5+_0","INS_repeats_5+_1","INS_repeats_5+_2","INS_repeats_5+_3","INS_repeats_5+_4","INS_repeats_5+_5+",
                          "DEL_MH_2_1","DEL_MH_3_1","DEL_MH_3_2","DEL_MH_4_1","DEL_MH_4_2","DEL_MH_4_3","DEL_MH_5+_1","DEL_MH_5+_2","DEL_MH_5+_3","DEL_MH_5+_4","DEL_MH_5+_5+")

#' @export
INDEL.MUTATION.TYPES.GROUPS <- c("1_BP_C_DELETION","1_BP_C_DELETION","1_BP_C_DELETION","1_BP_C_DELETION","1_BP_C_DELETION","1_BP_C_DELETION",
                                 "1_BP_T_DELETION","1_BP_T_DELETION","1_BP_T_DELETION","1_BP_T_DELETION","1_BP_T_DELETION","1_BP_T_DELETION",
                                 "1_BP_C_INSERTION","1_BP_C_INSERTION","1_BP_C_INSERTION","1_BP_C_INSERTION","1_BP_C_INSERTION","1_BP_C_INSERTION",
                                 "1_BP_T_INSERTION","1_BP_T_INSERTION","1_BP_T_INSERTION","1_BP_T_INSERTION","1_BP_T_INSERTION","1_BP_T_INSERTION",
                                 "2_BP_DELETION","2_BP_DELETION","2_BP_DELETION","2_BP_DELETION","2_BP_DELETION","2_BP_DELETION",
                                 "3_BP_DELETION","3_BP_DELETION","3_BP_DELETION","3_BP_DELETION","3_BP_DELETION","3_BP_DELETION",
                                 "4_BP_DELETION","4_BP_DELETION","4_BP_DELETION","4_BP_DELETION","4_BP_DELETION","4_BP_DELETION",
                                 "5_PLUS_BP_DELETION","5_PLUS_BP_DELETION","5_PLUS_BP_DELETION","5_PLUS_BP_DELETION","5_PLUS_BP_DELETION","5_PLUS_BP_DELETION",
                                 "2_BP_INSERTION","2_BP_INSERTION","2_BP_INSERTION","2_BP_INSERTION","2_BP_INSERTION","2_BP_INSERTION",
                                 "3_BP_INSERTION","3_BP_INSERTION","3_BP_INSERTION","3_BP_INSERTION","3_BP_INSERTION","3_BP_INSERTION",
                                 "4_BP_INSERTION","4_BP_INSERTION","4_BP_INSERTION","4_BP_INSERTION","4_BP_INSERTION","4_BP_INSERTION",
                                 "5_PLUS_BP_INSERTION","5_PLUS_BP_INSERTION","5_PLUS_BP_INSERTION","5_PLUS_BP_INSERTION","5_PLUS_BP_INSERTION","5_PLUS_BP_INSERTION",
                                 "2_BP_DELETION_MICROHOMOLOGY",
                                 "3_BP_DELETION_MICROHOMOLOGY","3_BP_DELETION_MICROHOMOLOGY",
                                 "4_BP_DELETION_MICROHOMOLOGY","4_BP_DELETION_MICROHOMOLOGY","4_BP_DELETION_MICROHOMOLOGY",
                                 "5_PLUS_BP_DELETION_MICROHOMOLOGY","5_PLUS_BP_DELETION_MICROHOMOLOGY","5_PLUS_BP_DELETION_MICROHOMOLOGY","5_PLUS_BP_DELETION_MICROHOMOLOGY","5_PLUS_BP_DELETION_MICROHOMOLOGY")

#' @export
PLOT.INDEL.X.AXIS.LABELS <- c("1","2","3","4","5","6+",
                              "1","2","3","4","5","6+",
                              "0","1","2","3","4","5+",
                              "0","1","2","3","4","5+",
                              "1","2","3","4","5","6+",
                              "1","2","3","4","5","6+",
                              "1","2","3","4","5","6+",
                              "1","2","3","4","5","6+",
                              "0","1","2","3","4","5",
                              "0","1","2","3","4","5",
                              "0","1","2","3","4","5",
                              "0","1","2","3","4","5",
                              "1","1","2","1","2","3",
                              "1","2","3","4","5+")

#' @export
PLOT.INDEL.STRIP.LABELS <- c("C", "T", "C", "T",
                             "2", "3", "4", "5+",
                             "2", "3", "4", "5+",
                             "2", "3", "4", "5+")

#' @export
PLOT.INDEL.STRIP.COLORS <- c("#FDBE6F", "#FF8001", "#B0DD8B", "#36A12E",
                             "#FDCAB5", "#FC8A6A", "#F14432", "#BC141A",
                             "#D0E1F2", "#94C4DF", "#4A98C9", "#1764AB",
                             "#E2E2EF", "#B6B6D8", "#8683BD", "#61409B")

#' @export
PLOT.INDEL.STRIP.TEXT.COLORS <- c("#000000", "#FFFFFF", "#000000", "#FFFFFF",
                                  "#000000", "#000000", "#000000", "#FFFFFF",
                                  "#000000", "#000000", "#000000", "#FFFFFF",
                                  "#000000", "#000000", "#000000", "#FFFFFF")

#' @export
PLOT.INDEL.LEGEND.LABELS <- c("Deletion of C (length)",
                              "Deletion of T (length)",
                              "Insertion of C (length)",
                              "Insertion of T (length)",
                              "2bp deletion at repeats (count)",
                              "3bp deletion at repeats (count)",
                              "4bp deletion at repeats (count)",
                              "5+bp deletion at repeats (count)",
                              "2bp insertion at repeats (count)",
                              "3bp insertion at repeats (count)",
                              "4bp insertion at repeats (count)",
                              "5+bp insertion at repeats (count)",
                              "2bp deletion with microhomology (length)",
                              "3bp deletion with microhomology (length)",
                              "4bp deletion with microhomology (length)",
                              "5+bp deletion with microhomology (length)")
