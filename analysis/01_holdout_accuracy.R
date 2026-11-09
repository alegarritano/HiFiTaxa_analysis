#!/usr/bin/env Rscript
.SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd())
})
source(file.path(.SCRIPT_DIR, "00_common.R"))
train_rank_sets <- function(ref90_path) {
  ref <- read_lineage_tsv(ref90_path)
  setNames(lapply(RANKS, function(r) unique(na.omit(ref[[r]]))), RANKS)
}

score_marker <- function(marker) {
  hd   <- file.path(BR, sprintf("benchmark_holdout_%s_seed42", marker))
  rdir <- file.path(hd, "results")
  truth <- read_lineage_tsv(file.path(hd,"test_10.taxonomy"))
  train <- train_rank_sets(file.path(hd,"reference_90.taxonomy"))
  bl <- read_blca(rdir); nb <- read_nb(rdir)

  ids <- intersect(truth$id, bl$tab$id)
  T  <- truth   %>% filter(id %in% ids) %>% arrange(id)
  B  <- bl$tab  %>% filter(id %in% ids) %>% arrange(id)
  Bc <- bl$conf %>% filter(id %in% ids) %>% arrange(id)
  N  <- nb      %>% filter(id %in% ids) %>% arrange(id)
  stopifnot(all(T$id==B$id), all(T$id==Bc$id), all(T$id==N$id))

  excl <- setNames(lapply(RANKS, function(r) rep(FALSE, nrow(T))), RANKS)
  if (marker=="its") {
    excl$Species <- !is.na(T$Species) & is_placeholder_sp(T$Species)
    excl$Genus   <- !is.na(T$Genus)   & is_placeholder_genus(T$Genus)
  }

  Bmask <- B
  for (r in RANKS) { cc<-suppressWarnings(as.numeric(Bc[[paste0(r,"_conf")]])); cc[is.na(cc)]<-0
                     Bmask[[r]][cc<80] <- NA }

  acc <- bind_rows(lapply(RANKS, function(r){
    tv <- T[[r]]
    keep <- !is.na(tv) & (tv %in% train[[r]]) & !excl[[r]]
    denom <- sum(keep)
    bind_rows(lapply(c("BLCA","BLCA@80","NB"), function(cl){
      pv <- switch(cl,"BLCA"=B[[r]],"BLCA@80"=Bmask[[r]],"NB"=N[[r]])
      cc <- sum(keep & !is.na(pv) & pv==tv); ic <- sum(keep & !is.na(pv) & pv!=tv); un <- sum(keep & is.na(pv))
      tibble(marker=marker, Rank=r, classifier=cl, n=denom, correct=cc,
             accuracy=ifelse(denom>0,100*cc/denom,NA_real_),
             pct_incorrect=ifelse(denom>0,100*ic/denom,NA_real_),
             pct_unassigned=ifelse(denom>0,100*un/denom,NA_real_))
    }))
  }))

  sp_orphan <- !is.na(T$Species) & !(T$Species %in% train$Species) & !excl$Species
  g_present <- !is.na(T$Genus)   &  (T$Genus   %in% train$Genus)   & !excl$Genus
  g_orphan  <- !is.na(T$Genus)   & !(T$Genus   %in% train$Genus)   & !excl$Genus
  f_present <- !is.na(T$Family)  &  (T$Family  %in% train$Family)
  backoff <- bind_rows(
    bind_rows(lapply(c("BLCA","BLCA@80","NB"), function(cl){
      idx <- sp_orphan & g_present; pv<-switch(cl,"BLCA"=B$Genus,"BLCA@80"=Bmask$Genus,"NB"=N$Genus)[idx]; tv<-T$Genus[idx]; n<-length(tv)
      tibble(marker=marker,classifier=cl,level="Species absent, scored at genus",n=n,
             Correct=100*sum(!is.na(pv)&pv==tv)/n,Incorrect=100*sum(!is.na(pv)&pv!=tv)/n,Unassigned=100*sum(is.na(pv))/n)})),
    bind_rows(lapply(c("BLCA","BLCA@80","NB"), function(cl){
      idx <- g_orphan & f_present; pv<-switch(cl,"BLCA"=B$Family,"BLCA@80"=Bmask$Family,"NB"=N$Family)[idx]; tv<-T$Family[idx]; n<-length(tv)
      tibble(marker=marker,classifier=cl,level="Genus absent, scored at family",n=n,
             Correct=100*sum(!is.na(pv)&pv==tv)/n,Incorrect=100*sum(!is.na(pv)&pv!=tv)/n,Unassigned=100*sum(is.na(pv))/n)})))
  sp_present <- !is.na(T$Species) & (T$Species %in% train$Species) & !excl$Species
  denom <- tibble(marker=marker,
    quantity=c("Test sequences scored",
               "Test sequences with a placeholder species label",
               "Test sequences with an Incertae sedis genus label",
               "Species accuracy denominator (species present in training)",
               "Genus accuracy denominator (genus present in training)",
               "Species absent, genus present (scored at genus)",
               "Genus absent, family present (scored at family)"),
    n=c(length(ids), sum(excl$Species), sum(excl$Genus), sum(sp_present), sum(g_present),
        sum(sp_orphan & g_present), sum(g_orphan & f_present)))
  list(acc=acc, backoff=backoff, denom=denom)
}

s16 <- score_marker("16s"); sits <- score_marker("its")
acc<-bind_rows(s16$acc,sits$acc); backoff<-bind_rows(s16$backoff,sits$backoff)
mlab<-c("16s"="16S (GTDB r232 holdout)","its"="Fungal ITS (UNITE holdout)")
acc$marker_lab<-mlab[acc$marker]; backoff$marker_lab<-mlab[backoff$marker]

write_csv(acc,     file.path(OUT,"tables","holdout_accuracy.csv"))
write_csv(backoff, file.path(OUT,"tables","holdout_orphan_backoff.csv"))
write_csv(sits$denom, file.path(OUT,"tables","holdout_its_placeholder_denominators.csv"))
cat("Wrote tables/holdout_accuracy.csv, holdout_orphan_backoff.csv and holdout_its_placeholder_denominators.csv\n")
