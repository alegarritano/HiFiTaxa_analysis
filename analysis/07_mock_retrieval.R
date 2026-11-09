#!/usr/bin/env Rscript
.SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd())
})
source(file.path(.SCRIPT_DIR, "00_common.R"))
detect <- function(relab, exp_set, rank) {
  d <- relab %>% mutate(taxon=if(rank=="Species") norm_species(taxon) else norm_genus(taxon)) %>%
       filter(!is.na(taxon)) %>% group_by(taxon) %>% summarise(relab=sum(relab),.groups="drop")
  det <- length(intersect(exp_set, d$taxon[d$relab>0]))
  fo <- d %>% filter(!(taxon %in% exp_set), relab>FP_THRESH)
  fo <- if (rank=="Species") fo %>% filter(!is_placeholder_sp(taxon)) else
        fo %>% filter(!str_detect(taxon,"(ceae$| sp$|^Fungi|Incertae|_fam_|_gen_)"))
  tibble(rank=rank, Detected=det, Missed=length(exp_set)-det, `False pos.`=nrow(fo),
         nonred=det+nrow(fo), n_expected=length(exp_set))
}

marker_detect <- function(marker) {
  if (marker=="16s") { dir<-file.path(BR,"results_atcc_16s"); eg<-atcc_genus; es<-atcc_species; rl<-"Emu" }
  else                { dir<-file.path(BR,"results_fungi_its"); eg<-fung_genus; es<-fung_species; rl<-"EMITS" }
  cnt<-read_feature_table(dir); bl<-read_blca(dir)$tab; b8<-read_blca_at80(dir); nb<-read_nb(dir)
  rd <- function(tab,col) relab_from_counts(tab$id, tab[[col]], cnt)
  em <- function(rk) if(marker=="16s") read_emu(dir,tolower(rk)) else read_emits(dir,tolower(rk))
  rows <- bind_rows(
    detect(rd(bl,"Genus"),    eg,"Genus")   %>% mutate(classifier="BLCA"),
    detect(rd(bl,"Species"),  es,"Species") %>% mutate(classifier="BLCA"),
    detect(rd(b8,"Genus"),    eg,"Genus")   %>% mutate(classifier="BLCA@80"),
    detect(rd(b8,"Species"),  es,"Species") %>% mutate(classifier="BLCA@80"),
    detect(em("Genus"),       eg,"Genus")   %>% mutate(classifier=rl),
    detect(em("Species"),     es,"Species") %>% mutate(classifier=rl),
    detect(rd(nb,"Genus"),    eg,"Genus")   %>% mutate(classifier="NB"),
    detect(rd(nb,"Species"),  es,"Species") %>% mutate(classifier="NB"))
  rows$marker <- marker; rows
}

d16 <- marker_detect("16s"); dIT <- marker_detect("its")
THR_TAG <- if (FP_THRESH == 0.001) "" else if (FP_THRESH == 0) "_threshold_any_read" else sprintf("_threshold_%gpct", 100 * FP_THRESH)
write_csv(bind_rows(d16,dIT) %>% select(marker,classifier,rank,Detected,Missed,`False pos.`,nonred,n_expected),
          file.path(OUT,"tables",paste0("mock_retrieval",THR_TAG,".csv")))

cat(sprintf("Wrote tables/mock_retrieval%s.csv (FP_THRESH = %g)\n", THR_TAG, FP_THRESH))
