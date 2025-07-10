# circular viz

genome_size <- 100000

com_nt <- read.table("outputs/result_blast_competence/blastn_tabular_nt.txt",
                     header = F, sep = "\t") %>% 
  stats::setNames(c("file_name", "qseqid", "temporary_sseqid",
                    "nt_pident", "nt_length", "nt_mismatch",
                    "nt_gapopen", "nt_qstart", "nt_qend",
                    "nt_sstart", "nt_send",
                    "nt_evalue", "nt_bitscore")) %>% 
  dplyr::left_join(
    read.table("inputs/prepare_competence_genes/competence_headers_nt.txt",
               header = F, sep = "\t") %>% 
      dplyr::rename(header = V1) %>% 
      dplyr::mutate(gene    = str_extract(header, "(?<=\\[gene=)[^\\]]+(?=\\])"), # "[gene=" & "]"
                    gene_id = str_extract(header, "^[^ ]+"), # begin & " "
                    prot_id = str_extract(header, "(?<=cds_)[^.]+"), # "cds_" & "_"
                    strain  = str_extract(header, "(?<=\\] \\[)[^\\]]+(?=\\])") # "] [" and "]"
      ) %>% 
      dplyr::mutate(prot_id = paste0(prot_id, ".1")) # correction
    ,
    by = c("temporary_sseqid" = "gene_id")
  ) %>% 
  dplyr::mutate(nt_lengthdiff = abs(abs(nt_sstart-nt_send) - abs(nt_qstart-nt_qend))) %>% 
  dplyr::filter(file_name == "Streptococcus_pneumoniae_RMD131_contigs_from_YM") %>% 
  view() %>% 
  glimpse()

coms <- data.frame(
  start = pmin(com_nt$nt_qstart, com_nt$nt_qend),
  end = pmax(com_nt$nt_qstart, com_nt$nt_qend),
  gene = com_nt$gene
) %>% 
  dplyr::distinct(gene, .keep_all = T) %>% 
  # view() %>% 
  glimpse()



gff <- read.delim("/home/ron/pneumo_2025_serotype1_brief_report/outputs/result_prokka_combined_pathogenWatch_monocle/Streptococcus_pneumoniae_RMD131_contigs_from_YM.gff",
                  comment.char = "#", header = FALSE) %>% 
  stats::setNames(c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")) %>% 
  dplyr::filter(grepl("comA|comB|comC|comD|comE|comX|comEA|comEC|comGC|coiA|=com|gene=com", attributes)) %>% 
  # view() %>% 
  glimpse()




library(circlize)

# Example genome size
genome_size <- 2000000

phages <- readLines("/home/ron/pneumo_2025_serotype1_brief_report/outputs/result_phastest/phastest_contigs/ZZ_9db106ee3a.PHASTEST/phage_regions.fna") %>% 
  grep("^>", ., value = TRUE) %>% 
  glimpse()

phage_df <- do.call(rbind, lapply(phages, function(header) {
  region <- sub("^>\\S+\\s+([^:]+):.*", "\\1", header)
  
  loc_match <- regmatches(header, regexec(":(\\d+)-(\\d+)", header))[[1]]
  start <- as.integer(loc_match[2])
  end <- as.integer(loc_match[3])
  
  data.frame(region = region, start = start, end = end, stringsAsFactors = FALSE)
})) %>% 
  glimpse()

circos.clear()
circos.initialize(factors = "chr", xlim = c(0, genome_size))

# Baseline genome track
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.axis(labels.cex = 0.01)  # Optional: axis labels
}, track.height = 0.05, bg.border = NA)


circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, bg.border = NA,
                       panel.fun = function(x, y) {
                         circos.rect(phage_df$start, 0, phage_df$end, 1, col = "purple", border = "purple")
                       })

circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, bg.border = NA,
                       panel.fun = function(x, y) {
                         circos.rect(coms$start, 0, coms$end, 1, col = "blue", border = "blue")
                       })

# Clear after plot
circos.clear()






regions <- dplyr::bind_rows(
  phage_df %>% 
    dplyr::mutate(category = "Phage")
  ,
  coms %>% 
    dplyr::mutate(category = "Competence")
)

regions <- regions %>%
  mutate(y = as.numeric(factor(category, levels = unique(category))))

# Plot
ggplot(regions) +
  geom_segment(aes(x = start, xend = end, y = y, yend = y, color = category), size = 6) +
  geom_text(aes(x = (start + end)/2, y = y + 0.3, label = gene), size = 3, vjust = 0) +
  scale_y_continuous(breaks = unique(regions$y), labels = unique(regions$category)) +
  scale_x_continuous("Genome Position (bp)", expand = expansion(mult = c(0.01, 0.01))) +
  labs(y = "Feature Type", title = "Linear Genome Annotation") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  )

