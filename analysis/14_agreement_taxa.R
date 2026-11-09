#!/usr/bin/env Rscript
.SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd())
})
source(file.path(.SCRIPT_DIR, "00_common.R"))
CATS <- c("Both correct","BLCA correct only","NB correct only",
          "Agree, both wrong","Disagree, both wrong","Unassigned")
CAT_PAL <- c("Both correct"="#06d6a0","BLCA correct only"="#118ab2",
             "NB correct only"="#ef476f","Agree, both wrong"="#ffd166",
             "Disagree, both wrong"="#073b4c","Unassigned"="#9aa0a6")
MINN <- 5
TOPN <- 100
LBL  <- 8

hd <- file.path(BR,"benchmark_holdout_16s_seed42"); rd <- file.path(hd,"results")
truth <- read_lineage_tsv(file.path(hd,"test_10.taxonomy"))
train <- read_lineage_tsv(file.path(hd,"reference_90.taxonomy"))
bl <- read_blca(rd)$tab; nb <- read_nb(rd)
ids <- Reduce(intersect, list(truth$id, bl$id, nb$id))
o <- function(x) x[match(ids,x$id),]; T<-o(truth); B<-o(bl); N<-o(nb)

tr_sp <- unique(na.omit(train$Species)); tr_ge <- unique(na.omit(train$Genus))
sp_in_ref <- !is.na(T$Species) & T$Species %in% tr_sp
ge_in_ref <- !is.na(T$Genus)   & T$Genus   %in% tr_ge

tv<-T$Genus; bv<-B$Genus; nv<-N$Genus
bc<-!is.na(bv)&bv==tv; nc<-!is.na(nv)&nv==tv
same<-!is.na(bv)&!is.na(nv)&bv==nv
out<-rep("Unassigned",length(tv))
out[!bc&!nc & !is.na(bv)&!is.na(nv)&bv!=nv]<-"Disagree, both wrong"
out[same & !bc & !nc]<-"Agree, both wrong"
out[bc&!nc]<-"BLCA correct only"; out[nc&!bc]<-"NB correct only"; out[bc&nc]<-"Both correct"

d <- tibble(taxon=tv, outcome=out,
            population = ifelse(sp_in_ref, "in_reference", "orphan"),
            genus_in_ref = ge_in_ref) %>% filter(!is.na(taxon))

genus_stats <- function(dd) {
  dd %>% group_by(taxon) %>%
    summarise(n = n(),
              genus_in_reference = first(genus_in_ref),
              both_correct       = sum(outcome=="Both correct"),
              blca_correct_only  = sum(outcome=="BLCA correct only"),
              nb_correct_only    = sum(outcome=="NB correct only"),
              agree_both_wrong   = sum(outcome=="Agree, both wrong"),
              disagree_both_wrong= sum(outcome=="Disagree, both wrong"),
              unassigned         = sum(outcome=="Unassigned"),
              .groups="drop") %>%
    mutate(agreement_pct    = round(100*both_correct/n, 2),
           disagreement_pct = round(100 - 100*both_correct/n, 2)) %>%
    arrange(desc(disagreement_pct), desc(n))
}

build_fig <- function(dd, pop) {
  per <- genus_stats(dd) %>%
    mutate(agree_lbl = ifelse(both_correct==n, 100L,
                              pmin(99L, as.integer(round(agreement_pct)))))
  scored <- per %>% filter(n>=MINN)
  sel   <- scored %>% filter(disagreement_pct > 0) %>%
           arrange(desc(disagreement_pct), desc(n)) %>% slice_head(n=TOPN)
  ny  <- nrow(sel)

  plotdat <- dd %>% filter(taxon %in% sel$taxon) %>%
    count(taxon, outcome, name="k") %>%
    group_by(taxon) %>% mutate(pct=100*k/sum(k)) %>% ungroup() %>%
    tidyr::complete(taxon, outcome=CATS, fill=list(k=0L, pct=0))
  plotdat$outcome <- factor(plotdat$outcome, levels=rev(CATS))
  ord <- sel %>% arrange(disagreement_pct) %>% pull(taxon)
  plotdat$taxon <- factor(plotdat$taxon, levels=ord)
  nlab <- sel %>% mutate(taxon=factor(taxon, levels=ord))

  PAGE_W <- 8.27
  CHROME <- 1.0
  nlbl   <- LBL * 0.82
  row_pt <- LBL / 0.95
  HEAD   <- 1.5
  PAGE_H <- CHROME + (ny + HEAD - 0.6) * row_pt / 72
  ytop   <- ny + 0.5

  p <- ggplot(plotdat, aes(pct, taxon, fill=outcome)) +
    geom_col(width=0.9, colour="white", linewidth=0.1) +
    geom_text(data=nlab, aes(x=103, y=taxon, label=format(n, big.mark=",")),
              inherit.aes=FALSE, hjust=0, size=nlbl/.pt, colour="grey35") +
    geom_text(data=nlab, aes(x=118, y=taxon, label=paste0(agree_lbl, "%")),
              inherit.aes=FALSE, hjust=0, size=nlbl/.pt, colour="grey35") +
    annotate("text", x=103, y=ytop, label="n", hjust=0, vjust=0, size=nlbl/.pt,
             fontface="bold", colour="grey30") +
    annotate("text", x=118, y=ytop, label="agree", hjust=0, vjust=0, size=nlbl/.pt,
             fontface="bold", colour="grey30") +
    coord_cartesian(clip="off") +
    scale_fill_manual(values=CAT_PAL, breaks=CATS, name=NULL, drop=FALSE) +
    scale_x_continuous("% of the genus's test sequences", limits=c(0,122),
                       breaks=seq(0,100,25), expand=expansion(mult=c(0,0))) +
    scale_y_discrete(expand=expansion(add=c(0.6, HEAD))) +
    labs(y=NULL) +
    guides(fill=guide_legend(nrow=1)) +
    theme_ms(12) +
    theme(legend.position="bottom", legend.text=element_text(size=8),
          legend.margin=margin(1,1,1,1), legend.box.spacing=unit(2,"pt"),
          legend.key.size=unit(9,"pt"),
          panel.grid=element_blank(), axis.ticks.y=element_blank(),
          plot.margin=margin(4, 8, 4, 4),
          axis.text.y=element_text(face="italic", size=LBL, margin=margin(r=1)))

  f <- file.path(OUT,"figures",sprintf("fig_agreement_taxa_16s_%s", pop))
  ggsave(paste0(f,".pdf"), p, width=PAGE_W, height=PAGE_H, bg="white")
}

d_ref <- d %>% filter(population=="in_reference")
write_csv(genus_stats(d_ref), file.path(OUT,"tables","agreement_genus_16s_in_reference.csv"))
build_fig(d_ref, "in_reference")

d_orp     <- d %>% filter(population=="orphan")
d_orp_fig <- d_orp %>% filter(genus_in_ref)
write_csv(genus_stats(d_orp), file.path(OUT,"tables","agreement_genus_16s_orphan.csv"))
build_fig(d_orp_fig, "orphan")

cat("Wrote tables/agreement_genus_16s_{in_reference,orphan}.csv and figures/fig_agreement_taxa_16s_{in_reference,orphan}.pdf\n")
