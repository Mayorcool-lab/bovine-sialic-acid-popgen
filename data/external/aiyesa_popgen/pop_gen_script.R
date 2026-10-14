library(dplyr)

vcf    <- read.table("vcf_samples.txt", stringsAsFactors = FALSE)
colnames(vcf)[1] <- "vcf_id"
sample <- read.table("sample_info.txt", header = TRUE, stringsAsFactors = FALSE)
sample$vcf_id_orig <- sample$vcf_id
vcf_ids <- vcf$vcf_id

# candidate VCF IDs for each metadata ID (exact match wins, otherwise prefix match)
cands <- lapply(sample$vcf_id_orig, function(id) {
  if (id %in% vcf_ids) return(id)
  vcf_ids[startsWith(vcf_ids, id)]
})
sample$n_cand <- lengths(cands)
table(sample$n_cand)      # 1 = clean, 0 = no match, >1 = ambiguous

# apply only the unambiguous ones
sample$vcf_id <- ifelse(sample$n_cand == 1, unlist(lapply(cands, `[`, 1)), NA)

# look at the ambiguous ones
sample %>% filter(n_cand > 1) %>% select(vcf_id_orig, breed, n_cand)
cands[sample$n_cand > 1]

# checks
anyDuplicated(sample$vcf_id[!is.na(sample$vcf_id)])   # must be 0
sum(!is.na(sample$vcf_id))                            # samples with metadata now

# 1. Confirm what is left
matched   <- sample$vcf_id[!is.na(sample$vcf_id)]
left_vcf  <- setdiff(vcf$vcf_id, matched)
length(left_vcf)                                  # expect 127
table(sub("^([A-Za-z]+).*", "\\1", left_vcf))     # rough prefix summary
left_vcf[grepl("^(NL|PT)", left_vcf)]             # these should be the ambiguous ones (~7)

# 2. Which metadata rows were ambiguous
sample %>% filter(n_cand > 1) %>% select(vcf_id_orig, breed, sex, region)

# 3. Save the matched metadata
meta_ok <- sample %>% filter(!is.na(vcf_id)) %>%
  select(vcf_id, vcf_id_orig, breed, country, sex, region)
write.table(meta_ok, "sample_info_matched.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# 4. Group sizes for the 402
table(meta_ok$region)
sort(table(meta_ok$breed))


# `sample` already has:
#   vcf_id_orig = ID as written in the metadata (truncated)
#   vcf_id      = matched full ID from the VCF (NA if unmatched)

# Optional: recover the duplicate-ID cases, but only when the rows are
# indistinguishable (same breed, sex, region) and #rows == #candidates
for (id in unique(sample$vcf_id_orig[sample$n_cand > 1])) {
  rows <- which(sample$vcf_id_orig == id)
  cand <- cands[[rows[1]]]
  same <- nrow(unique(sample[rows, c("breed", "sex", "region")])) == 1
  if (length(rows) == length(cand) && same) sample$vcf_id[rows] <- cand
}

meta_final <- sample %>%
  transmute(meta_id  = vcf_id_orig,       # original metadata ID
            vcf_id   = vcf_id,            # exact ID in the VCF
            breed, country, sex, region,
            in_vcf   = !is.na(vcf_id))

# Sanity checks
anyDuplicated(meta_final$vcf_id[!is.na(meta_final$vcf_id)])   # must be 0
all(meta_final$vcf_id[!is.na(meta_final$vcf_id)] %in% vcf$vcf_id)  # must be TRUE
table(meta_final$in_vcf)

# Save everything (tab-delimited, no quotes) and the matched subset
write.table(meta_final, "sample_info_with_vcf_id.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(filter(meta_final, in_vcf), "sample_info_matched.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

m <- read.delim("sample_info_with_vcf_id.txt", stringsAsFactors = FALSE)
keep <- m$vcf_id[m$in_vcf %in% c(TRUE, "TRUE") & !is.na(m$vcf_id)]
length(keep); anyDuplicated(keep)
vcf_ids <- readLines("vcf_samples.txt")
setdiff(keep, vcf_ids)                     # must be character(0)

# PLINK2 format: FID IID (VCF import sets FID = 0 unless --double-id is used)
writeLines(paste0("0\t", keep), "keep_matched.txt")

##### MEta data for nu. div ----

library(dplyr); library(purrr)

m <- read.delim("sample_info_with_vcf_id.txt", stringsAsFactors = FALSE) %>%
  filter(in_vcf %in% c(TRUE, "TRUE"), !is.na(vcf_id), vcf_id != "NA") %>%
  distinct(vcf_id, .keep_all = TRUE)

clean <- function(x) gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", x))
dir.create("groups", showWarnings = FALSE)
unlink("groups/*.txt")                       # remove old lists so nothing stale is picked up

make_groups <- function(level) {
  m %>% group_by(label = .data[[level]]) %>%
    summarise(ids = list(vcf_id), n = n(), .groups = "drop") %>%
    mutate(level = level, group = paste0(level, "_", clean(label)),
           file  = paste0("groups/", group, ".txt")) %>%
    { pwalk(list(.$ids, .$file), function(i, f) writeLines(i, f)); . } %>%
    select(level, group, label, n, file)
}

manifest <- bind_rows(map(c("region", "country", "breed"), make_groups))
write.table(manifest, "groups_manifest.txt", sep = "\t", quote = FALSE, row.names = FALSE)

table(manifest$level)
manifest %>% arrange(level, n) %>% print(n = 100)



##### PCA ------


ev <- read.table("pca_matched.eigenvec", header = TRUE, comment.char = "", check.names = FALSE)
colnames(ev)[colnames(ev) == "#IID"] <- "IID"
colnames(ev)[colnames(ev) == "#FID"] <- "FID"
ev$vcf_id <- ev$IID
pca_meta <- merge(ev, m[, c("vcf_id","meta_id","breed","country","sex","region")], by = "vcf_id")
nrow(pca_meta); sum(is.na(pca_meta$breed))

frac <- read.table("pca_matched.eigenval")

library(ggplot2)
library(dplyr)
library(scales)

pve <- NULL

plot_pca <- function(df, colour_by, pcx = "PC1", pcy = "PC2", pve = NULL,
                     point_size = 1.8, alpha = 0.5,
                     jitter_frac = 0.005,   # jitter as a fraction of the axis range
                     legend_cols = 1, seed = 1) {
  
  n <- n_distinct(df[[colour_by]])
  
  pal <- if (n <= 8) {
    RColorBrewer::brewer.pal(max(3, n), "Dark2")[1:n]
  } else {
    hue_pal(l = 55, c = 100)(n)
  }
  
  # axis labels with % variance
  lab <- function(pc) {
    if (is.null(pve)) return(pc)
    sprintf("%s (%.1f%%)", pc, pve[as.integer(sub("PC", "", pc))])
  }
  
  # jitter scaled to the data range, so it is tiny relative to the plot
  jw <- jitter_frac * diff(range(df[[pcx]], na.rm = TRUE))
  jh <- jitter_frac * diff(range(df[[pcy]], na.rm = TRUE))
  
  ggplot(df, aes(x = .data[[pcx]], y = .data[[pcy]], colour = .data[[colour_by]])) +
    geom_point(size = point_size, alpha = alpha,
               position = position_jitter(width = jw, height = jh, seed = seed)) +
    scale_colour_manual(values = pal, name = tools::toTitleCase(colour_by)) +
    guides(colour = guide_legend(ncol = legend_cols,
                                 override.aes = list(size = 3, alpha = 1))) +
    labs(x = lab(pcx), y = lab(pcy)) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = 8))
}


p_region  <- plot_pca(pca_meta, "region",  pve = pve)
p_country <- plot_pca(pca_meta, "country", pve = pve, legend_cols = 2)
p_breed   <- plot_pca(pca_meta, "breed",   pve = pve, legend_cols = 3)

library(data.table)
pi <- fread("pi_res_breed.txt")

library(ggplot2)
library(dplyr)

# Four publication colours (Okabe-Ito, colour-blind safe)
pal4 <- c("#0072B2", "#D55E00", "#009E73", "#E69F00")

# If `pi` has no region column, add it from your metadata first:
# pi <- pi %>% left_join(distinct(pca_meta, country, region), by = "country")

# Order countries by median pi so the plot reads cleanly
pi <- pi %>%
  mutate(country = reorder(country, mean_pi, FUN = median))

p_pi <- ggplot(pi, aes(x = country, y = mean_pi, fill = region)) +
  geom_violin(width = 0.9, alpha = 0.35, colour = NA, trim = TRUE) +
#  geom_boxplot(width = 0.1, outlier.shape = NA,       # points are drawn separately
   #            colour = "grey20", linewidth = 0.4, alpha = 0.9,
    #           position = position_dodge(width = 0.9)) +
  geom_jitter(aes(colour = region), width = 0.06, size = 1.3, alpha = 0.6,
              show.legend = FALSE) +
  scale_fill_manual(values = pal4, name = "Region") +
  scale_colour_manual(values = pal4) +
  labs(x = NULL,
       y = expression("Mean nucleotide diversity ( "*pi*")")) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = c(0.9,0.15))

p_pi


pi[, Country := factor(Country)]
pi[, Group   := factor(Group)]

p <- ggplot(pi, aes(x = Country, y = mean_pi, fill = Group)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.5,
    outlier.shape = NA,          # points are drawn below instead
    alpha = 0.7
  ) +
  geom_point(
    aes(group = Group),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    size = 2, alpha = 0.8, shape = 21, colour = "black"
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = NULL,
    y = "Mean nucleotide diversity (π)",
    fill = "Group"
  ) + lims(y = c(0.12, 0.28)) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(face = "bold")
  )

p


##### Fst Outlier ----- 

library(data.table)
library(ggplot2)

breed <- fread("region_all.weir.fst")
breed_NA <- subset(breed,  WEIR_AND_COCKERHAM_FST == "NaN")
breed_NoNA <- subset(breed,  WEIR_AND_COCKERHAM_FST != "NaN")
setnames(breed_NoNA, "WEIR_AND_COCKERHAM_FST", "FST")

breed_NoNA[FST < 0, FST := 0]    

thr <- breed_NoNA[, .(
  n     = .N,
  mean  = mean(FST),
  median = median(FST),
  q95   = quantile(FST, 0.95),
  q99   = quantile(FST, 0.99),
  q999  = quantile(FST, 0.999)
)]
print(thr)

breed_NoNA[, outlier := FST >= thr$q999]
breed_outliers <- breed_NoNA[outlier == TRUE]

library(data.table)
library(ggplot2)

d <- copy(breed_NoNA)
d[, CHROM := factor(CHROM, levels = as.character(1:29))]

# cumulative x position
chr_len <- d[, .(len = max(POS)), by = CHROM][order(CHROM)]
chr_len[, off := cumsum(as.numeric(shift(len, fill = 0)))]
chr_len[, mid := off + len / 2]
d <- merge(d, chr_len[, .(CHROM, off)], by = "CHROM")
d[, x := POS + off]
d[, chr_col := as.integer(CHROM) %% 2 == 0]

p <- ggplot(d, aes(x, FST)) +
  geom_point(aes(colour = chr_col), size = 0.3, alpha = 0.5, show.legend = FALSE) +
  geom_point(data = d[outlier == TRUE], colour = "red", size = 0.5) +
  geom_hline(yintercept = thr$q999, linetype = "dashed", colour = "blue") +
  scale_colour_manual(values = c("grey30", "grey65")) +
  scale_x_continuous(breaks = chr_len$mid, labels = chr_len$CHROM, expand = c(0.005, 0)) +
  labs(x = "Chromosome", y = expression(F[ST]),
       title = "Breed-level Fst per site") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(size = 8, angle = 60, hjust = 1, vjust = 1))

#ggsave("breed_fst_manhattan.png", p, width = 12, height = 5, dpi = 300)
p

#### conserved SNPs across -----

coord <- read.table("gene_coord.txt", header = T)

library(data.table)
library(ggplot2)

setDT(d)
setDT(coord)

d[, gene := NA_character_]

for (i in 1:nrow(coord)) {
  
  d[
    CHROM == coord$chrom[i] &
      POS >= coord$start[i] &
      POS <= coord$stop[i],
    gene := coord$gene[i]
  ]
}

d[, in_gene := !is.na(gene)]
d[, gene_outlier := in_gene & outlier]

gene_labels[, label_y := max(d$FST, na.rm = TRUE) * 1.05]

gene_labels <- copy(coord)

gene_labels <- merge(
  gene_labels,
  chr_len[, .(
    CHROM_int = as.integer(as.character(CHROM)),
    off
  )],
  by.x = "chrom",
  by.y = "CHROM_int",
  all.x = TRUE
)

gene_labels[, gene_start_x := start + off]
gene_labels[, gene_end_x := stop + off]
gene_labels[, gene_mid_x := (gene_start_x + gene_end_x) / 2]

gene_labels[, label_y := max(d$FST, na.rm = TRUE) * 1.05]

for (i in 1:nrow(coord)) {
  
  d[
    CHROM == coord$chrom[i] &
      POS >= coord$start[i] &
      POS <= coord$stop[i],
    gene := coord$gene[i]
  ]
}

d[!is.na(gene), .(gene, CHROM, POS, FST)]

p <- ggplot(d, aes(x, FST)) +
  
  # All SNPs
  geom_point(
    aes(colour = chr_col),
    size = 0.3,
    alpha = 0.5,
    show.legend = FALSE
  ) +
  
  # SNPs within genes
  geom_point(
    data = d[in_gene == TRUE],
    aes(colour = gene),
    size = 1.2,
    alpha = 0.9
  ) +
  
  # FST outliers
  geom_point(
    data = d[outlier == TRUE],
    colour = "red",
    size = 0.5
  ) +
  
  # Outliers within genes
  geom_point(
    data = d[gene_outlier == TRUE],
    aes(colour = gene),
    size = 2
  ) +
  
  # Gene names
  geom_text(
    data = gene_labels,
    aes(
      x = gene_mid_x,
      y = label_y,
      label = gene,
      colour = gene
    ),
    size = 3,
    fontface = "bold",
    inherit.aes = FALSE,
    angle = 90,
    vjust = 0.5
  ) +
  
  geom_hline(
    yintercept = thr$q999,
    linetype = "dashed",
    colour = "blue"
  ) +
  
  scale_colour_manual(
    values = c(
      "1" = "grey30",
      "2" = "grey65",
      "CMAH" = "#E41A1C",
      "CMAS" = "#377EB8",
      "GNE" = "#4DAF4A",
      "NANP" = "#984EA3",
      "NANS" = "#FF7F00",
      "NEU1" = "#A65628",
      "NEU3" = "#F781BF",
      "SLC35A1" = "#999999",
      "ST3GAL1" = "#66C2A5",
      "ST3GAL4" = "#FC8D62",
      "ST6GAL1" = "#8DA0CB",
      "ST8SIA2" = "#E78AC3",
      "ST8SIA4" = "#A6D854"
    )
  ) +
  
  scale_x_continuous(
    breaks = chr_len$mid,
    labels = chr_len$CHROM,
    expand = c(0.005, 0)
  ) +
  
  labs(
    x = "Chromosome",
    y = expression(F[ST]),
    title = "Region-level FST per site
    "
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(
      size = 8,
      angle = 60,
      hjust = 1,
      vjust = 1
    ),
    legend.position = "none"
  )

p



###### Tajima's D plot -----

library(ggplot2)
library(dplyr)

taj_breed <- read.delim("TajimaD_by_breed.tsv")

taj_breed <- taj_breed %>%
  mutate(
    Region = recode(
      Region,
      "SA" = "South Africa",
      "EU" = "Europe"
    ),
    TajimaD = as.numeric(TajimaD)
  ) %>%
  filter(!is.na(TajimaD))

# Calculate regional means
library(dplyr)
library(ggplot2)

# Calculate mean Tajima's D for each region
region_means <- taj_breed %>%
  group_by(Region) %>%
  summarise(
    mean_TajimaD = mean(TajimaD, na.rm = TRUE)
  )

# Set different y positions for the mean labels
region_means <- region_means %>%
  mutate(
    label_y = ifelse(
      Region == "Europe",
      0.35,
      0.42
    )
  )

ggplot(taj_breed, aes(x = TajimaD, fill = Region)) +
  
  geom_density(alpha = 0.35) +
  
  scale_fill_manual(
    values = c(
      "Europe" = "blue",
      "South Africa" = "red"
    )
  ) +
  
  # Mean lines
  geom_vline(
    data = region_means,
    aes(
      xintercept = mean_TajimaD,
      color = Region
    ),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  # Mean labels at different heights
  geom_text(
    data = region_means,
    aes(
      x = mean_TajimaD,
      y = label_y,
      label = paste0(
        "Mean = ",
        round(mean_TajimaD, 2)
      ),
      color = Region
    ),
    size = 4,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      "Europe" = "blue",
      "South Africa" = "red"
    )
  ) +
  
  labs(
    x = "Tajima's D",
    y = "Density",
    fill = "Region"
  ) +
  
  lims(
    x = c(-2.5, 4),
    y = c(0, 0.45)
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = c(0.88, 0.9)
  )

