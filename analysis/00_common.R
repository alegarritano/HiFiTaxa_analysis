suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(stringr); library(readr)
  library(patchwork); library(scales)
})

BR   <- Sys.getenv("HIFITAXA_DATA", unset = "")
OUT  <- Sys.getenv("HIFITAXA_OUT",  unset = file.path(.SCRIPT_DIR, "regenerated"))
REPO <- Sys.getenv("HIFITAXA_REPO", unset = "")
if (REPO == "") for (.cand in c(normalizePath(file.path(.SCRIPT_DIR, ".."), mustWork = FALSE),
                                normalizePath(file.path(.SCRIPT_DIR, "..", ".."), mustWork = FALSE))) {
  if (file.exists(file.path(.cand, "data", "fungi_103_species.txt"))) { REPO <- .cand; break }
}
if (BR == "" || !dir.exists(BR)) stop(sprintf(
  "Benchmark data not found (HIFITAXA_DATA = '%s').
Set HIFITAXA_DATA to the directory holding the run directories listed in analysis/README.md:
  HIFITAXA_DATA=/path/to/run_dirs Rscript <script>.R", BR))
if (REPO == "") stop(
  "HIFITAXA_REPO is not set and this script is not running from analysis/ of a HiFiTaxa_analysis checkout.
Set HIFITAXA_REPO to a clone of https://github.com/alegarritano/HiFiTaxa_analysis; it provides the
fungal mock truth list data/fungi_103_species.txt.")

dir.create(file.path(OUT, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "tables"),  recursive = TRUE, showWarnings = FALSE)
RANKS <- c("Domain","Phylum","Class","Order","Family","Genus","Species")
TRUTH_KEY <- c(Domain="superkingdom",Phylum="phylum",Class="class",Order="order",
               Family="family",Genus="genus",Species="species")
FP_THRESH <- as.numeric(Sys.getenv("FP_THRESH", unset = "0.001"))
if (is.na(FP_THRESH) || FP_THRESH < 0) stop("FP_THRESH must be a non negative number, for example 0.001")

strip_ani <- function(x) str_replace_all(x, "_[A-Z]+(?=(\\s|$))", "")
GENUS_ALIAS <- c(G047199095="Escherichia",ECMA0423="Escherichia",Pauljensenia="Schaalia",
                 Thaxterogaster="Cortinarius",Phlegmacium="Cortinarius",Calonarius="Cortinarius",
                 Rugosomyces="Calocybe")
SPECIES_ALIAS <- c("G047199095 sp047199095"="Escherichia coli",
                   "Pseudomonas aeruginosa"="Pseudomonas paraeruginosa",
                   "Pauljensenia dentiphila"="Schaalia odontolytica",
                   "Schaalia dentiphila"="Schaalia odontolytica",
                   "Bacillus cereus"="Bacillus pacificus",
                   "Cortinarius ochraceobrunneum"="Cortinarius ochraceobrunneus",
                   "Phlegmacium ochraceobrunneum"="Cortinarius ochraceobrunneus",
                   "Rugosomyces obscurissimus"="Calocybe obscurissima")
KEEP_GTDB_NAME_GENUS <- c("ECMA0423","G047199095")
clean <- function(x){x<-str_trim(x);x[x %in% c("","NA","Unclassified","unclassified","s__","g__","NA NA")]<-NA;x}
norm_genus <- function(g){g<-strip_ani(str_trim(g));ifelse(g %in% names(GENUS_ALIAS),GENUS_ALIAS[g],g)}
norm_species <- function(s){
  s<-str_trim(s); s<-ifelse(s %in% names(SPECIES_ALIAS),SPECIES_ALIAS[s],s); s<-strip_ani(s)
  s<-ifelse(s %in% names(SPECIES_ALIAS),SPECIES_ALIAS[s],s)
  t<-str_split_fixed(s," ",2); g<-ifelse(t[,1] %in% names(GENUS_ALIAS) & !(t[,1] %in% KEEP_GTDB_NAME_GENUS),GENUS_ALIAS[t[,1]],t[,1])
  o<-ifelse(t[,2]=="",g,paste(g,t[,2])); ifelse(o %in% names(SPECIES_ALIAS),SPECIES_ALIAS[o],o)}
is_placeholder_sp <- function(s){ep<-str_split_fixed(s," ",2)[,2]
  ep=="" | str_detect(ep,"^(sp|spp|unidentified|uncultured|indet|incertae)\\.?($|\\s)") |
    str_detect(s,regex("(\\ssp\\.?$|incertae[ _]sedis)",ignore_case=TRUE))}
is_placeholder_genus <- function(g) str_detect(g,regex("incertae[ _]sedis",ignore_case=TRUE))

read_lineage_tsv <- function(path) {
  raw <- read_tsv(path, col_names=c("id","lin"), show_col_types=FALSE, progress=FALSE)
  out <- tibble(id=raw$id)
  for (r in RANKS) out[[r]] <- clean(str_match(raw$lin, paste0("(?:^|;)",TRUTH_KEY[[r]],":([^;]*)"))[,2])
  out
}

read_blca <- function(dir) {
  tab <- read_csv(file.path(dir,"taxonomy_blca","blca_taxonomy_table.csv"), show_col_types=FALSE, progress=FALSE); names(tab)[1]<-"id"
  for (r in RANKS) tab[[r]] <- clean(tab[[r]])
  conf <- read_csv(file.path(dir,"taxonomy_blca","blca_taxonomy_confidence.csv"), show_col_types=FALSE, progress=FALSE); names(conf)[1]<-"id"
  list(tab=tab, conf=conf)
}
read_nb <- function(dir) {
  nb <- read_tsv(file.path(dir,"taxonomy_nb","best_taxonomy.tsv"), show_col_types=FALSE, progress=FALSE); names(nb)[1:2]<-c("id","Taxon")
  prefix <- c(Domain="(?:d|k)__",Phylum="p__",Class="c__",Order="o__",Family="f__",Genus="g__",Species="s__")
  out <- tibble(id=nb$id)
  for (r in RANKS) out[[r]] <- clean(str_match(nb$Taxon, paste0(prefix[[r]],"([^;]*)"))[,2])
  out
}

read_feature_table <- function(dir) {
  f <- file.path(dir,"dada2","feature_table.tsv")
  ft <- read_tsv(f, comment="# Constructed", show_col_types=FALSE, progress=FALSE)
  names(ft)[1] <- "id"
  samp <- setdiff(names(ft), "id")
  tibble(id=ft$id, count=rowSums(sapply(ft[samp], as.numeric), na.rm=TRUE))
}
relab_from_counts <- function(ids, taxa, counts) {
  d <- tibble(id=ids, taxon=taxa) %>% inner_join(counts, by="id") %>%
       filter(!is.na(taxon)) %>% group_by(taxon) %>% summarise(count=sum(count),.groups="drop")
  d$relab <- d$count/sum(d$count); d %>% select(taxon,relab)
}
read_emu <- function(dir, rank) {
  f <- file.path(dir,"taxonomy_emu",sprintf("emu_%s_table.tsv",rank))
  e <- read_tsv(f, show_col_types=FALSE, progress=FALSE)
  samp <- setdiff(names(e), c("species","genus","family","order","class","phylum","superkingdom"))
  tax <- clean(str_replace(e[[rank]], paste0("^",rank,":"), ""))
  em <- as.matrix(e[samp]); em[is.na(em)] <- 0
  m <- rowMeans(em); tibble(taxon=tax, relab=m/sum(m)) %>% filter(!is.na(taxon))
}
read_emits <- function(dir, rank) {
  fs <- list.files(file.path(dir,"taxonomy_emits","per_sample"),
                   pattern=sprintf("_emits_%s\\.tsv$",rank), full.names=TRUE)
  pref <- if (rank=="species") "s__" else "g__"
  rows <- list()
  for (f in fs) {
    d <- suppressWarnings(read_tsv(f, col_names=c("taxstr","relab","reads"),
                                   show_col_types=FALSE, progress=FALSE))
    d <- d[!is.na(d$reads),]
    lin <- str_split_fixed(d$taxstr, "\\|", 5)[,5]
    tax <- str_replace_all(clean(str_match(lin, paste0(pref,"([^;]*)"))[,2]), "_", " ")
    rows[[f]] <- tibble(taxon=tax, reads=as.numeric(d$reads))
  }
  d <- bind_rows(rows) %>% filter(!is.na(taxon)) %>% group_by(taxon) %>%
       summarise(reads=sum(reads),.groups="drop")
  d$relab <- d$reads/sum(d$reads); d %>% select(taxon,relab)
}

atcc_species <- unique(norm_species(c(
  "Acinetobacter baumannii","Bacillus pacificus","Phocaeicola vulgatus","Bifidobacterium adolescentis",
  "Clostridium beijerinckii","Cutibacterium acnes","Deinococcus radiodurans","Enterococcus faecalis",
  "Escherichia coli","Helicobacter pylori","Lactobacillus gasseri","Neisseria meningitidis",
  "Porphyromonas gingivalis","Pseudomonas paraeruginosa","Cereibacter sphaeroides","Schaalia odontolytica",
  "Staphylococcus aureus","Staphylococcus epidermidis","Streptococcus agalactiae","Streptococcus mutans")))
atcc_genus <- unique(norm_genus(str_split_fixed(atcc_species," ",2)[,1]))
ATCC_TRUTH <- c("Acinetobacter baumannii"=0.0018,"Bacillus pacificus"=0.0180,"Phocaeicola vulgatus"=0.0002,
  "Bifidobacterium adolescentis"=0.0002,"Clostridium beijerinckii"=0.0180,"Cutibacterium acnes"=0.0018,
  "Deinococcus radiodurans"=0.0002,"Enterococcus faecalis"=0.0002,"Escherichia coli"=0.1800,
  "Helicobacter pylori"=0.0018,"Lactobacillus gasseri"=0.0018,"Neisseria meningitidis"=0.0018,
  "Porphyromonas gingivalis"=0.1800,"Pseudomonas paraeruginosa"=0.0180,"Cereibacter sphaeroides"=0.1800,
  "Schaalia odontolytica"=0.0002,"Staphylococcus aureus"=0.0180,"Staphylococcus epidermidis"=0.1800,
  "Streptococcus agalactiae"=0.0180,"Streptococcus mutans"=0.1800)
.fung_truth <- file.path(REPO,"data/fungi_103_species.txt")
if (!file.exists(.fung_truth)) stop(sprintf(
  "Mock truth list not found at '%s'.
Set HIFITAXA_REPO to a HiFiTaxa_analysis checkout (git clone https://github.com/alegarritano/HiFiTaxa_analysis).
The list is read whenever 00_common.R is sourced, so it is required even by the holdout scripts.", .fung_truth))
.fr <- readLines(.fung_truth); .fr<-str_trim(.fr)
.fr <- .fr[.fr!="" & !str_starts(.fr,"#")]
fung_species <- unique(norm_species(.fr))
fung_genus   <- unique(norm_genus(str_split_fixed(fung_species," ",2)[,1]))

CLF_COL <- c("Expected"="#888888",
             "BLCA"="#0072B2","BLCA@80"="#56B4E9","Emu"="#D55E00",
             "NB"="#CC79A7","EMITS"="#009E73")
CAT_COL <- c("Correct"="#44AA99","Incorrect"="#CC6677","Unassigned"="#BBBBBB")
TAXA_PAL <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7",
              "#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#882255")
OTHER_COL <- "#cccccc"; UNASSIGNED_COL <- "#888888"
NOT_TAXA <- c("Other","Unassigned")
taxon_labels <- function(x) parse(text=ifelse(x %in% NOT_TAXA, sprintf('plain("%s")', x), sprintf('italic("%s")', x)))
theme_ms <- function(base=12) theme_classic(base_size=base) +
  theme(plot.title=element_text(face="bold",hjust=0),
        strip.background=element_rect(fill="grey95",colour="black"),
        strip.text=element_text(face="bold"))

read_blca_at80 <- function(dir) {
  b <- read_blca(dir); tab <- b$tab; cf <- b$conf
  cf <- cf[match(tab$id, cf$id),]
  for (r in RANKS) { cc<-suppressWarnings(as.numeric(cf[[paste0(r,"_conf")]])); cc[is.na(cc)]<-0
                     tab[[r]][cc<80] <- NA }
  tab
}

comp_reads <- function(ids, taxa, counts) {
  d <- tibble(id=ids, taxon=taxa) %>% inner_join(counts, by="id")
  d$taxon[is.na(d$taxon)] <- "Unassigned"
  d %>% group_by(taxon) %>% summarise(count=sum(count),.groups="drop") %>%
    mutate(relab=count/sum(count)) %>% select(taxon,relab)
}
