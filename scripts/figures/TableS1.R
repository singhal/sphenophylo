rm(list = ls())

library(ape)
library(tidyverse)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_Phylogeny2/")
t = read.tree("../Sphenomorphine_speciation/phylogeny/dating/synthetic_tree.all.rooted.dated_smooth0.1.18April2023.tre")

d = read.csv("../Sphenomorphine_speciation/metadata/sphenomorphine_all_v13.csv")
outs = d[which(d$INGROUP == FALSE), "SAMPLE_ID"]

d1 = d[!d$SAMPLE_ID %in% outs, ]
d1 = d1[d1$SAMPLE_ID %in% t$tip.label, ]

# nominal taxon
d1$NOMINAL_TAXON = d1$MORPHOLOGICAL_TAXON

d1[which(d1$OPERATIONAL_TAXON == "Ctenotus_euclae"), "NOMINAL_TAXON"] = "Ctenotus_euclae"
d1[which(d1$OPERATIONAL_TAXON == "Ctenotus_kutjupa"), "NOMINAL_TAXON"] = "Ctenotus_kutjupa"
d1[which(d1$OPERATIONAL_TAXON == "Ctenotus_maryani"), "NOMINAL_TAXON"] = "Ctenotus_maryani"
d1[which(d1$OPERATIONAL_TAXON == "Ctenotus_olympicus"), "NOMINAL_TAXON"] = "Ctenotus_olympicus"
d1[which(d1$OPERATIONAL_TAXON == "Ctenotus_orientalis"), "NOMINAL_TAXON"] = "Ctenotus_orientalis"
d1[which(d1$OPERATIONAL_TAXON == "Ctenotus_pallasotus"), "NOMINAL_TAXON"] = "Ctenotus_pallasotus"
d1[which(d1$OPERATIONAL_TAXON == "Ctenotus_rhabdotus"), "NOMINAL_TAXON"] = "Ctenotus_rhabdotus"
d1[which(d1$OPERATIONAL_TAXON %in% c("Ctenotus superciliaris_W", "Ctenotus_superciliaris_E")), "NOMINAL_TAXON"] = "Ctenotus_superciliaris"
d1[which(d1$OPERATIONAL_TAXON == "Eremiascincus_musivus"), "NOMINAL_TAXON"] = "Eremiascincus_musivus"
d1[which(d1$OPERATIONAL_TAXON == "Lerista_emmotti"), "NOMINAL_TAXON"] = "Lerista_emmotti"
d1[which(d1$OPERATIONAL_TAXON == "Lerista_micra"), "NOMINAL_TAXON"] = "Lerista_micra"

# threshold taxon
x = read.csv("data/threshold_taxonomy/threshold_taxa.csv")
d1$THRESHOLD_TAXON = x[match(d1$SAMPLE_ID, x$ind), "SPECIES"]
  
# change ddrad
d1[d1$ddRAD == TRUE & d1$RECENT_RUN == TRUE, "ddRAD"] = FALSE
  
# select columns
d2 = d1 %>% select(SAMPLE_ID, NOMINAL_TAXON, LAT, LON,
             MORPHOLOGICAL_TAXON, OPERATIONAL_TAXON, 
           INCIPIENT_TAXON, THRESHOLD_TAXON, CYTB, SQCL, ddRAD)

# ddrad sra
s1 = do.call("rbind", lapply(list.files("Table_S1/", pattern = "ddRAD", full.names = T), read.csv))
d2$ddRAD_SRA = NA
s1$number = parse_number(s1$sample)
for (i in 1:nrow(d2)) {
  if (d2[i, "ddRAD"]) {
    name = d2[i, "SAMPLE_ID"]
    if (name %in% s1$sample) {
      d2[i, "ddRAD_SRA"] = s1[s1$sample == name, "sra"]
    } else {
      s2 = s1[s1$number == parse_number(name), ]
      if (nrow(s2) == 1) {
        d2[i, "ddRAD_SRA"] = s2$sra
        # message(name, " ", s2$sample[1])
      }
    }
  }
}

d2[d2$SAMPLE_ID == "NA_DLR0504_Le_bipe", "ddRAD_SRA"] = "SRS3431581"
d2[d2$SAMPLE_ID == "NA_DLR0506_Le_bipe", "ddRAD_SRA"] = "SRS3431580"
d2[d2$SAMPLE_ID == "NA_CCM1582_Le_bore", "ddRAD_SRA"] = "SRS3431692"
d2[d2$SAMPLE_ID == "CUMV_14639_Le_dese", "ddRAD_SRA"] = "SRS3431764"
d2[d2$SAMPLE_ID == "CUMV_14666_Le_dese", "ddRAD_SRA"] = "SRS3431756"
d2[d2$SAMPLE_ID == "NA_CCM0766_Le_kalu", "ddRAD_SRA"] = "SRS3431722"
d2[d2$SAMPLE_ID == "NA_CCM0334_Le_sp", "ddRAD_SRA"] = "SRS3431623"
d2[d2$SAMPLE_ID == "NA_CCM1203_Le_sp", "ddRAD_SRA"] = "SRS3431819"

# sqcl sra
s1 = do.call("rbind", lapply(list.files("Table_S1/", pattern = "SqCL", full.names = T), read.csv))
d2$sqcl_SRA = NA
s1$number = parse_number(s1$sample)
for (i in 1:nrow(d2)) {
  if (d2[i, "SQCL"]) {
    name = d2[i, "SAMPLE_ID"]
    if (name %in% s1$sample) {
      d2[i, "sqcl_SRA"] = s1[s1$sample == name, "sra"]
    } else {
      s2 = s1[which(s1$number == parse_number(name)), ]
      if (nrow(s2) == 1) {
        d2[i, "sqcl_SRA"] = s2$sra
        # message(name, " ", s2$sample[1])
      } 
    }
  }
}

# undo wrong ones
d2[d2$SAMPLE_ID == "NA_CCM2913_Ct_essi", "sqcl_SRA"] = NA
d2[d2$SAMPLE_ID == "NTMR_22387_Ct_mili", "sqcl_SRA"] = NA
d2[d2$SAMPLE_ID == "NA_K041_Ct_stor", "sqcl_SRA"] = NA

# genbank acc
d2$CYTB_genbank = NA
x = readxl::read_xlsx("Table_S1/Table_S2_samples_cytb.xlsx")
d2$CYTB_genbank = pull(x[match(d2$SAMPLE_ID, x$`Sample ID`), "GenBank Accession"])

x = list.files("Table_S1/", pattern = "genbank_seqids", full.names = T)
x = lapply(x, function(y) read.csv(y, header = F, sep = "\t"))
x = do.call("rbind", x)

d2$CYTB_genbank = coalesce(d2$CYTB_genbank, x[match(d2$SAMPLE_ID, x$V1), "V2"])

s1 = do.call("rbind", lapply(list.files("Table_S1/", pattern = "mtDNA", full.names = T), read.csv))
s1$number = parse_number(s1$sample)
for (i in 1:nrow(d2)) {
  if (d2[i, "CYTB"]) {
    if (!complete.cases(d2[i, "CYTB_genbank"])) {
      name = d2[i, "SAMPLE_ID"]
      if (name %in% s1$sample) {
        d2[i, "CYTB_genbank"] = pull(s1[s1$sample == name, "name"])
      } else {
        s2 = s1[s1$number == parse_number(name), ]
        if (nrow(s2) == 1) {
          d2[i, "CYTB_genbank"] = s2$name[1]
          # message(name, " ", s2$sample[1], " ", s2$name[1])
        }
      }
    }
  }
}

d2[d2$SAMPLE_ID == "WAMR_158107_Ct_duri", "CYTB_genbank"] = "KJ505341.1"
d2[d2$SAMPLE_ID == "WAMR_163003_Ct_duri", "CYTB_genbank"] = "KJ505356.1"
d2[d2$SAMPLE_ID == "WAMR_170880_Ct_duri", "CYTB_genbank"] = "KJ505368.1"
d2[d2$SAMPLE_ID == "WAMR_159964_Ct_duri", "CYTB_genbank"] = "KJ505343.1"
d2[d2$SAMPLE_ID == "WAMR_160163_Ct_duri", "CYTB_genbank"] = "KJ505346.1"
d2[d2$SAMPLE_ID == "WAMR_161378_Ct_duri", "CYTB_genbank"] = "KJ505350.1"
d2[d2$SAMPLE_ID == "WAMR_161625_Ct_duri", "CYTB_genbank"] = "KJ505353.1"
d2[d2$SAMPLE_ID == "WAMR_170219_Ct_duri", "CYTB_genbank"] = "KJ505361.1"
d2[d2$SAMPLE_ID == "WAMR_170489_Ct_duri", "CYTB_genbank"] = "KJ505364.1"
d2[d2$SAMPLE_ID == "WAMR_170638_Ct_duri", "CYTB_genbank"] = "KJ505365.1"
d2[d2$SAMPLE_ID == "WAMR_127728_Ct_scho", "CYTB_genbank"] = "OM966748.1"
d2[d2$SAMPLE_ID == "WAMR_152013_Ct_tant", "CYTB_genbank"] = "KJ505819.1"

# to check
# NA_MILISP092_Ct_mili SP092 KJ505514.1 - OKAY
# NA_PANTSP044_Ct_pant CCM:0044 OQ091801.1 - NOT OKAY
d2[d2$SAMPLE_ID == "NA_PANTSP044_Ct_pant", "CYTB_genbank"] = NA
# NA_K045_Ct_essi EBUA45 KJ505047.1 - NOT OKAY
d2[d2$SAMPLE_ID == "NA_K045_Ct_essi", "CYTB_genbank"] = NA
# NA_ANGUW103931_Ct_angu W103931 KJ505252.1 - OKAY

# also need to look at cases where there are multiple matches

write.csv(d2, "Table_S1/Table_S1.2October24.csv", row.names = F)

nrow(d2[d2$ddRAD == TRUE & !complete.cases(d2$ddRAD_SRA), ])
nrow(d2[d2$SQCL == TRUE & !complete.cases(d2$sqcl_SRA), ])
nrow(d2[d2$CYTB == TRUE & !complete.cases(d2$CYTB_genbank), ])
