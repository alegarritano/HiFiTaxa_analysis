#!/usr/bin/env Rscript
.SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd())
})
source(file.path(.SCRIPT_DIR, "00_common.R"))
acc <- read_csv(file.path(OUT,"tables","holdout_accuracy.csv"),show_col_types=FALSE,progress=FALSE)

panel <- function(mk, title, leg=TRUE) {
  d <- acc %>% filter(marker==mk) %>%
       transmute(Rank, classifier,
                 Correct=accuracy, Incorrect=pct_incorrect, Unassigned=pct_unassigned) %>%
       pivot_longer(c(Correct,Incorrect,Unassigned), names_to="cat", values_to="pct")
  d$Rank <- factor(d$Rank, levels=RANKS)
  if (mk=="its") levels(d$Rank)[levels(d$Rank)=="Domain"] <- "Kingdom"
  d$cat <- factor(d$cat, levels=c("Correct","Incorrect","Unassigned"))
  d$classifier <- factor(d$classifier, levels=c("BLCA","BLCA@80","NB"))

  ggplot(d, aes(Rank, pct, fill=cat)) +
    geom_col(position=position_dodge(0.8), width=0.75) +
    geom_text(aes(label=ifelse(pct >= 1, round(pct), "")),
              position=position_dodge(0.8), angle=90, hjust=-0.2, vjust=0.5,
              size=7/.pt, colour="grey20") +
    facet_wrap(~classifier) +
    scale_fill_manual(values=CAT_COL, name=NULL) +
    scale_y_continuous(limits=c(0,116), breaks=seq(0,100,25),
                       expand=expansion(mult=c(0,0.02))) +
    labs(x=NULL, y="% of test sequences", title=title) +
    theme_ms(12) +
    theme(legend.position=if(leg) "bottom" else "none", legend.text=element_text(size=10),
          panel.grid=element_blank(),
          axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
}

fig <- panel("16s","a", FALSE) /
       panel("its","b", TRUE)
ggsave(file.path(OUT,"figures","fig_holdout_decomposition.pdf"), fig, width=10, height=9)

cat("Wrote figures/fig_holdout_decomposition.pdf\n")
