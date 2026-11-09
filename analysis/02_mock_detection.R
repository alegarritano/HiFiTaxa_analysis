#!/usr/bin/env Rscript
.SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd())
})
source(file.path(.SCRIPT_DIR, "00_common.R"))
score_detect <- function(pred, exp_set, rank) {
  pred <- pred %>% mutate(taxon=if(rank=="Species") norm_species(taxon) else norm_genus(taxon)) %>%
          filter(!is.na(taxon)) %>% group_by(taxon) %>% summarise(relab=sum(relab),.groups="drop")
  detected  <- pred$taxon[pred$relab>0]
  recovered <- intersect(exp_set, detected); missed <- setdiff(exp_set, detected)
  foreign <- pred %>% filter(!(taxon %in% exp_set), relab>FP_THRESH)
  if (rank=="Species") foreign <- foreign %>% filter(!is_placeholder_sp(taxon)) else
                       foreign <- foreign %>% filter(!str_detect(taxon,"(ceae$| sp$|^Fungi|Incertae|_fam_|_gen_)"))
  missed_str <- paste(sort(missed), collapse="; ")
  tibble(rank=rank, n_expected=length(exp_set), recovered=length(recovered), missed=length(missed),
         false_positive=nrow(foreign), missed_taxa=missed_str,
         fp_relab=paste(sprintf("%s(%.1f%%)",foreign$taxon,100*foreign$relab),collapse="; "))
}

run_marker <- function(marker) {
  if (marker=="16s") {
    dir<-file.path(BR,"results_atcc_16s"); cnt<-read_feature_table(dir)
    bl<-read_blca(dir)$tab; nb<-read_nb(dir)
    res<-bind_rows(
      score_detect(relab_from_counts(bl$id,norm_genus(bl$Genus),cnt),   atcc_genus,  "Genus")  %>% mutate(classifier="BLCA"),
      score_detect(relab_from_counts(bl$id,norm_species(bl$Species),cnt),atcc_species,"Species")%>% mutate(classifier="BLCA"),
      score_detect(relab_from_counts(nb$id,norm_genus(nb$Genus),cnt),   atcc_genus,  "Genus")  %>% mutate(classifier="NB"),
      score_detect(relab_from_counts(nb$id,norm_species(nb$Species),cnt),atcc_species,"Species")%>% mutate(classifier="NB"),
      score_detect(read_emu(dir,"genus"),                              atcc_genus,  "Genus")  %>% mutate(classifier="Emu"),
      score_detect(read_emu(dir,"species"),                            atcc_species,"Species") %>% mutate(classifier="Emu"))
    res$marker<-"16S (ATCC MSA-1003)"
  } else {
    dir<-file.path(BR,"results_fungi_its"); cnt<-read_feature_table(dir)
    bl<-read_blca(dir)$tab; nb<-read_nb(dir)
    res<-bind_rows(
      score_detect(relab_from_counts(bl$id,norm_genus(bl$Genus),cnt),   fung_genus,  "Genus")  %>% mutate(classifier="BLCA"),
      score_detect(relab_from_counts(bl$id,norm_species(bl$Species),cnt),fung_species,"Species")%>% mutate(classifier="BLCA"),
      score_detect(relab_from_counts(nb$id,norm_genus(nb$Genus),cnt),   fung_genus,  "Genus")  %>% mutate(classifier="NB"),
      score_detect(relab_from_counts(nb$id,norm_species(nb$Species),cnt),fung_species,"Species")%>% mutate(classifier="NB"),
      score_detect(read_emits(dir,"genus"),                            fung_genus,  "Genus")  %>% mutate(classifier="EMITS"),
      score_detect(read_emits(dir,"species"),                          fung_species,"Species") %>% mutate(classifier="EMITS"))
    res$marker<-"Fungal ITS (103-sp mock)"
  }
  res
}

det <- bind_rows(run_marker("16s"), run_marker("its"))
write_csv(det%>%select(marker,classifier,rank,n_expected,recovered,missed,false_positive,missed_taxa,fp_relab),
          file.path(OUT,"tables","mock_detection.csv"))

cat("Wrote tables/mock_detection.csv\n")
