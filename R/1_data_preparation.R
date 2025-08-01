library(tidyverse)
# Genome analyses from Prokka
prokka_folder = "outputs/result_prokka_combined_pathogenWatch_monocle/"

# Genome (*.fna)
genome <- Biostrings::readDNAStringSet(paste0(prokka_folder, "Streptococcus_pneumoniae_RMD131.fna"))
cat("Genome width\n")
print(Biostrings::width(genome))

# Predicted genes (*.ffn)
genes <- Biostrings::readDNAStringSet(paste0(prokka_folder, "Streptococcus_pneumoniae_RMD131.ffn"))
gene_lengths <- Biostrings::width(genes)
genes_plot <- ggplot(data = data.frame(gene_lengths), aes(x = gene_lengths)) +
  geom_histogram(binwidth = 50, fill = "steelblue", color = "white") + 
  labs(title = "Gene Length Distribution",
       x = "Gene Length (bp)",
       y = "Count") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()

# Predicted proteins (*.faa)
proteins <- Biostrings::readAAStringSet(paste0(prokka_folder, "Streptococcus_pneumoniae_RMD131.faa"))
protein_lengths <- Biostrings::width(proteins)
prot_plot <- ggplot(data = data.frame(protein_lengths), aes(x = protein_lengths)) +
  geom_histogram(binwidth = 50, fill = "steelblue", color = "white") + 
  labs(title = "Protein Length Distribution",
       x = "Protein Length (aa)",
       y = "Count") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()

# Functional annotation (*.tsv)
annotations <- read.delim(paste0(prokka_folder, "Streptococcus_pneumoniae_RMD131.tsv"), header = TRUE, sep = "\t") %>% 
  dplyr::count(product) %>% 
  dplyr::top_n(10)
# head(annotations)

fun_plot <- ggplot(annotations, aes(x = reorder(product, n), y = n)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = n), hjust = -0.1, size = 1.5) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 6)) +
  scale_y_continuous(limits = c(0, max(annotations$n) + 2)) + 
  labs(title = "Top 10 Gene Functions",
       x = "Function",
       y = "Count",
       )

p_genes <- cowplot::plot_grid(plot(genes_plot),
                              ncol = 1,
                              labels = c("A"))

p_proteins <- cowplot::plot_grid(prot_plot, fun_plot,
                                 ncol = 1,
                                 labels = c("B", "C"))

p_gene_prot <- cowplot::plot_grid(p_genes, p_proteins,
                                  ncol = 2)

png("pictures/genomes.png",
    width = 24, height = 10, unit = "cm", res = 600)
p_gene_prot
dev.off()

# mixed infection, GPSC and MLST profile
kity_folder = "outputs/result_pneumokity/fastq_mix/"
mlst_folder = "outputs/result_mlst/"
poppunk_folder = "outputs/result_poppunk/"

pneumokity_result <- read.csv(paste0(kity_folder, "Collated_result_data.csv"))
cat("Predicted serotype\n")
print(pneumokity_result$predicted_serotype)

mlst_result <- read.csv(paste0(mlst_folder, "mlst_results.csv"), sep = "\t", header = FALSE,
                        col.names = c("name", "species", "mlst",
                                      "aroE", "gdh", "gki", "recP", "spi", "xpt", "ddl")) %>% 
  dplyr::mutate(name = gsub("/home/ron/Strep_Indo/raw_data/serotype_1/prepare_choosen_fasta_pathogenWatch_monocle/", "", name),
                aroE = substr(aroE, 6, nchar(aroE) - 1),
                gdh = substr(gdh, 5, nchar(gdh) - 1),
                gki = substr(gki, 5, nchar(gki) - 1),
                recP = substr(recP, 6, nchar(recP) - 1),
                spi = substr(spi, 5, nchar(spi) - 1),
                xpt = substr(xpt, 5, nchar(xpt) - 1),
                ddl = substr(ddl, 5, nchar(ddl) - 1)) %>% 
  dplyr::select(-species) %>% 
  dplyr::filter(name == "Streptococcus_pneumoniae_RMD131.fasta") %>% 
  glimpse()

poppunk_result <- read.csv(paste0(poppunk_folder, "result_poppunk_clusters.csv"), sep = ",") %>% 
  glimpse()

# antimicrobial profile
cdc_folder = "outputs/result_cdc_amr/"

# Penicillin-related genes mutation type
cdc_result_pen <- read.table(paste0(cdc_folder, "TABLE_Isolate_Typing_results.txt"),
                             header = TRUE) %>% 
  dplyr::select(Sample, Pili, dplyr::contains("PBP"))

cat("Penicillin binding proteins mutation nomenclature\n")
cdc_result_pen

# All AMR result
cdc_result <- read.table(paste0(cdc_folder, "TABLE_Isolate_Typing_results.txt"),
                         header = TRUE) %>% 
  dplyr::select(-WGS_Serotype, -ST, -aroe, -gdh, -gki, -recP, -spi, -xpt, -ddl) %>% 
  dplyr::mutate_all(as.character) %>% 
  tidyr::pivot_longer(cols = everything(),
                      names_to = "Tests", values_to = "Result")

# Generate new df
cdc_names <- cdc_result %>% 
  dplyr::mutate(Antibiotic = case_when(
    grepl("PEN", Tests) ~ "Penicillin (meningitis & non-meningitis)",
    grepl("AMO", Tests) ~ "Amoxicillin",
    grepl("MER", Tests) ~ "Meropenem",
    grepl("TAX", Tests) ~ "Cefotaxime (meningitis & non-meningitis)",
    grepl("CFT", Tests) ~ "Cefotaxime (meningitis & non-meningitis)",
    grepl("CFX", Tests) ~ "Ceftriaxone",
    grepl("AMP", Tests) ~ "Ampicillin",
    grepl("CPT", Tests) ~ "Cefepime",
    grepl("ZOX", Tests) ~ "Ceftizoxime",
    grepl("EC", Tests) ~ "Erythromycin-Clindamycin",
    grepl("ERY", Tests) ~ "Erythromycin",
    grepl("CLI", Tests) ~ "Clindamycin",
    grepl("SYN", Tests) ~ "Synercid",
    grepl("LZO", Tests) ~ "Linezolid",
    grepl("COT", Tests) ~ "Trimethoprim-Sulfamethoxazole",
    grepl("TET", Tests) ~ "Tetracycline",
    grepl("DOX", Tests) ~ "Doxycycline",
    grepl("CIP", Tests) ~ "Ciprofloxacin",
    grepl("LFX", Tests) ~ "Levofloxacin",
    grepl("CHL", Tests) ~ "Chloramphenicol",
    grepl("RIF", Tests) ~ "Rifampin",
    grepl("VAN", Tests) ~ "Vancomycin",
    grepl("DAP", Tests) ~ "Daptomycin",
    TRUE ~ NA_character_
  ),
  Abbreviation = case_when(
    grepl("PEN", Tests) ~ "PEN",
    grepl("AMO", Tests) ~ "AMO",
    grepl("MER", Tests) ~ "MER",
    grepl("TAX", Tests) ~ "TAX",
    grepl("CFT", Tests) ~ "CFT",
    grepl("CFX", Tests) ~ "CFX",
    grepl("AMP", Tests) ~ "AMP",
    grepl("CPT", Tests) ~ "CPT",
    grepl("ZOX", Tests) ~ "ZOX",
    grepl("EC", Tests) ~ "EC",
    grepl("ERY", Tests) ~ "ERY",
    grepl("CLI", Tests) ~ "CLI",
    grepl("SYN", Tests) ~ "SYN",
    grepl("LZO", Tests) ~ "LZO",
    grepl("COT", Tests) ~ "COT",
    grepl("TET", Tests) ~ "TET",
    grepl("DOX", Tests) ~ "DOX",
    grepl("CIP", Tests) ~ "CIP",
    grepl("LFX", Tests) ~ "LFX",
    grepl("CHL", Tests) ~ "CHL",
    grepl("RIF", Tests) ~ "RIF",
    grepl("VAN", Tests) ~ "VAN",
    grepl("DAP", Tests) ~ "DAP",
    TRUE ~ NA_character_
  )
  )

cdc_result_wide <- cdc_names %>% 
  dplyr::mutate(Type = case_when(
    grepl("_SIGN", Tests) ~ "SIGN",
    grepl("_SIR", Tests) ~ "SIR",
    TRUE ~ "MIC"
  )) %>%
  dplyr::select(Antibiotic, Abbreviation, Type, Result) %>%
  tidyr::pivot_wider(names_from = Type, values_from = Result,
                     values_fn = list(Result = ~ paste(unique(.), collapse = ", "))) %>% 
  dplyr::mutate(MIC = paste0(SIGN, " ", MIC),
                MIC = case_when(
                  MIC == "NA NA" ~ NA_character_,
                  TRUE ~ MIC
                ),
                SIR = case_when(
                  SIR == "NA" ~ NA_character_,
                  TRUE ~ SIR
                )) %>% 
  dplyr::filter(!is.na(Antibiotic)) %>% 
  dplyr::select(-SIGN)

cdc_result_wide



# virulence profile
vir_combined_results <- dplyr::left_join(
  read.table("outputs/result_blast_virulences/blastn_tabular_setA_nt_VFDB.txt",
             header = F, sep = "\t") %>% 
    stats::setNames(c("file_name", "qseqid", "temporary_sseqid",
                      "nt_pident", "nt_length", "nt_mismatch",
                      "nt_gapopen", "nt_qstart", "nt_qend",
                      "nt_sstart", "nt_send",
                      "nt_evalue", "nt_bitscore")) %>% 
    dplyr::left_join(
      read.table("inputs/prepare_vfdb_database/VFDB_setA_compiled_headers_nt.txt",
                 header = F, sep = "\t") %>% 
        dplyr::rename(header = V1) %>% 
        dplyr::mutate(gene    = str_extract(header, "(?<= \\()[^\\)]+(?=\\) )"), # " (" & ") "
                      gene_id = str_extract(header, "^[^)]*\\)"), # begin & ")"
                      gene_class = str_extract(header, "(?<= - )[^\\(]+(?= \\()"),  # " - " & " ("
                      species = str_extract(header, "(?<=\\] \\[)[^\\]]+(?=\\])") # "] [" & "]"
        )
      ,
      by = c("temporary_sseqid" = "gene_id")
    ) %>% 
    dplyr::mutate(
      temporary_sseqid = case_when(
        temporary_sseqid == "VFG005653" ~ "VFG005653(gb|WP_142355754.1)", # different ID for nanA from VFDB (gene) and NCBI (amino acids)
        TRUE ~ temporary_sseqid
        ),
      # nt_sdiff = abs(nt_sstart-nt_send),
      # nt_qdiff = abs(nt_qstart-nt_qend),
      nt_lengthdiff = abs(abs(nt_sstart-nt_send) - abs(nt_qstart-nt_qend))
    )%>% 
    dplyr::filter(file_name == "Streptococcus_pneumoniae_RMD131_contigs_from_YM")
  ,
  read.table("outputs/result_blast_virulences/blastx_tabular_setA_aa_VFDB.txt",
             header = F, sep = "\t") %>% 
    stats::setNames(c("file_name", "qseqid", "temporary_sseqid",
                      "aa_pident", "aa_length", "aa_mismatch",
                      "aa_gapopen", "aa_qstart", "aa_qend",
                      "aa_sstart", "aa_send",
                      "aa_evalue", "aa_bitscore")) %>% 
    dplyr::left_join(
      read.table("inputs/prepare_vfdb_database/VFDB_setA_compiled_headers_pro.txt",
                 header = F, sep = "\t") %>% 
        dplyr::rename(header = V1) %>% 
        dplyr::mutate(gene    = str_extract(header, "(?<= \\()[^\\)]+(?=\\) )"), # " (" & ") "
                      gene_id = str_extract(header, "^[^)]*\\)"), # begin & ")"
                      species = str_extract(header, "(?<=\\] \\[)[^\\]]+(?=\\])") # "] [" & "]"
        )
      ,
      by = c("temporary_sseqid" = "gene_id")
    ) %>% 
    dplyr::filter(file_name == "Streptococcus_pneumoniae_RMD131_contigs_from_YM")
  ,
  by = c("file_name", "qseqid", "temporary_sseqid")
) %>% 
  dplyr::arrange(dplyr::desc(aa_pident)) %>% 
  dplyr::distinct(file_name, temporary_sseqid, .keep_all = T) %>% 
  dplyr::select(-contains(".y")) %>% 
  dplyr::mutate(
    gene_present = case_when(
      nt_pident >= 70 ~ "present",
      TRUE ~ "absent"
    ),
    aa_lengthdiff = abs(abs(aa_sstart-aa_send) - abs(aa_qstart-aa_qend)),
    # protein_function = case_when(
    #   aa_pident >= 95 & aa_gapopen == 0 ~ "functional",
    #   aa_pident >= 90 & aa_gapopen <= 5 & aa_mismatch <= 30 ~ "variant",
    #   aa_pident <= 85 & aa_mismatch >= 25 ~ "possibly defective",
    #   TRUE ~ "possibly defective"
    # ),
    bacwgs_check = case_when(
      gene.x %in% c("cbpD", "cps4A", "cps4B", "cps4D", "hysA", "lytA", "lytC", 
                    "nanA", "nanB", "pavA", "pce/cbpE", "pfbA", "ply", "psaA") ~ "detected",
      TRUE ~ "no"
    )
  ) %>% 
  # nanA
  dplyr::mutate(
    gene.x = ifelse(is.na(gene.x), "nanA", gene.x),
    gene_class = ifelse(is.na(gene_class), "Exoenzyme", gene_class),
    species.x = ifelse(is.na(species.x), "Streptococcus pneumoniae NCTC7465 (serotype 1)", species.x)
  ) %>% 
  view() %>%
  glimpse()

vir_blastp <- read.table("outputs/result_blastp_virulences_additional_analysis/Streptococcus_pneumoniae_RMD131_contigs_from_YM_setA_blastp_results.txt",
                         header = F, sep = "\t") %>% 
  stats::setNames(c("qseqid", "temporary_sseqid",
                    "aap_pident", "aap_length",
                    # "aap_mismatch",
                    # "aap_gapopen", 
                    "aap_qstart", "aap_qend",
                    "aap_sstart", "aap_send",
                    "aap_evalue", "aap_bitscore", "qlen")) %>% 
  dplyr::left_join(
    read.table("inputs/prepare_vfdb_database/VFDB_setA_compiled_headers_pro.txt",
               header = F, sep = "\t") %>% 
      dplyr::rename(header = V1) %>% 
      dplyr::mutate(gene    = str_extract(header, "(?<= \\()[^\\)]+(?=\\) )"), # " (" & ") "
                    gene_id = str_extract(header, "^[^)]*\\)"), # begin & ")"
                    gene_class = str_extract(header, "(?<= - )[^\\(]+(?= \\()"),  # " - " & " ("
                    species = str_extract(header, "(?<=\\] \\[)[^\\]]+(?=\\])") # "] [" & "]"
      )
    ,
    by = c("temporary_sseqid" = "gene_id")
  ) %>% 
  dplyr::mutate(
    temporary_sseqid = case_when(
      temporary_sseqid == "VFG005653" ~ "VFG005653(gb|WP_142355754.1)", # different ID for nanA from VFDB (gene) and NCBI (amino acids)
      TRUE ~ temporary_sseqid
    ),
    # aap_sdiff = abs(aap_sstart-aap_send),
    # aap_qdiff = abs(aap_qstart-aap_qend),
    aap_lengthdiff = abs(abs(aap_sstart-aap_send) - abs(aap_qstart-aap_qend)),
    
    alignment_length = abs(aap_qstart-aap_qend) + 1,
    query_coverage = alignment_length / qlen
  ) %>% 
  dplyr::mutate(
    gene = ifelse(is.na(gene), "nanA", gene),
    gene_class = ifelse(is.na(gene_class), "Exoenzyme", gene_class),
    species = ifelse(is.na(species), "Streptococcus pneumoniae NCTC7465 (serotype 1)", species)
  ) %>% 
  # dplyr::filter(
  #   gene == "pspA"
  # ) %>% 
  dplyr::filter(aap_evalue <= 1e-3,
                aap_lengthdiff <= 70,
                aap_bitscore >= 50,
                aap_pident >= 70,
                query_coverage >= 0.7
                ) %>%
  dplyr::group_by(gene, temporary_sseqid, species, query_coverage) %>%
  dplyr::summarise(
    hits = n(),
    mean_coverage = mean(query_coverage)*100,
    min_coverage = min(query_coverage)*100,
    max_coverage = max(query_coverage)*100,
    
    mean_identity = mean(aap_pident),
    min_identity = min(aap_pident),
    max_identity = max(aap_pident),
    
    mean_bitscore = mean(aap_bitscore),
    min_bitscore = min(aap_bitscore),
    max_bitscore = max(aap_bitscore),
    
    mean_evalue = mean(aap_evalue),
    min_evalue = min(aap_evalue),
    max_evalue = max(aap_evalue),
  ) %>%
  dplyr::arrange(desc(hits)) %>%
  dplyr::mutate(
    coverage = ifelse(
      hits > 1, 
      paste0(round(mean_coverage, 3), "% (", min_coverage, "%-", max_coverage, "%)"), 
      paste0(round(mean_coverage, 3), "%")
    ),
    avg_identity = ifelse(
      hits > 1, 
      paste0(round(mean_identity, 3), "% (", min_identity, "%-", max_identity, "%)"), 
      paste0(round(mean_identity, 3), "%")
    ),
    avg_evalue = paste0(round(mean_evalue, 3), " (", min_evalue, "-", max_evalue, ")")
    ,
    bitscore = ifelse(
      hits > 1, 
      paste0(round(mean_bitscore, 3), " (", min_bitscore, "-", max_bitscore, ")"), 
      paste0(round(mean_bitscore, 3), "")
    )
  ) %>% 
  dplyr::select(gene, temporary_sseqid, species, hits,
                coverage,
                avg_identity,
                # avg_evalue,
                bitscore) %>% 
  view() %>%
  glimpse()

# test length
length <- vir_combined_results %>% 
  dplyr::select(contains(c("gene", "start", "end"))) %>% 
  dplyr::transmute(
    q_nt_diff = abs(nt_qstart - nt_qend),
    s_nt_diff = abs(nt_sstart - nt_send),
    q_aa_diff = abs(aa_qstart - aa_qend),
    s_aa_diff = abs(aa_sstart - aa_send),
    
    test_nt_diff = abs(q_nt_diff - s_nt_diff),
    test_aa_diff = abs(q_aa_diff - s_aa_diff)
  ) %>% 
  # view() %>% 
  glimpse()

# nanA is not included in VFDB pro list -_-)
report_blastp <- dplyr::left_join(
  vir_blastp
  ,
  vir_combined_results %>% 
    select(gene_class, gene.x, #species.x,
           # nt_pident, nt_lengthdiff, nt_mismatch, nt_evalue,
           # aa_pident, aa_lengthdiff, aa_mismatch, aa_gapopen, aa_evalue,
           # gene_present, #protein_function, 
           # bacwgs_check
    ) #%>% 
    # dplyr::distinct(gene.x, .keep_all = T)
  ,
  by = c("gene" = "gene.x")
) %>% 
  dplyr::ungroup() %>%
  dplyr::select(gene_class, gene, species,
                hits, coverage, avg_identity, bitscore
                ) %>% 
  dplyr::arrange(gene_class) %>% 
  # view() %>%
  glimpse()

# use blastp result instead
write.csv(report_blastp, "report/report_blastp_virulence.csv", row.names = F)

# crosscheck with panvita result
vir_mtx <- read.csv("outputs/result_panvita/Results_vfdb_25-04-2025_14-48-16/matriz_vfdb.csv", sep = ";") %>% 
  dplyr::select(-X) %>% 
  tidyr::pivot_longer(cols = 2:ncol(.),
                      names_to = "gene",
                      values_to = "aa_percent") %>% 
  dplyr::mutate(gene = gene %>%
                  str_replace("srtC", "srtC-") %>%
                  str_replace("-.", "-") %>%
                  str_replace("\\.", "/"),
                panvita_check = "detected"
  ) %>% 
  dplyr::full_join(
    report_blastp
    ,
    by = "gene"
  ) %>% 
  # view() %>% 
  glimpse()

# competence genes profile
com_combined_results <- dplyr::left_join(
  read.table("outputs/result_blast_competence/blastn_tabular_nt.txt",
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
    dplyr::filter(file_name == "Streptococcus_pneumoniae_RMD131_contigs_from_YM")
  ,
  read.table("outputs/result_blast_competence/blastx_tabular_aa.txt",
             header = F, sep = "\t") %>% 
    stats::setNames(c("file_name", "qseqid", "temporary_sseqid",
                      "aa_pident", "aa_length", "aa_mismatch",
                      "aa_gapopen", "aa_qstart", "aa_qend",
                      "aa_sstart", "aa_send",
                      "aa_evalue", "aa_bitscore")) %>% 
    dplyr::left_join(
      read.table("inputs/prepare_competence_genes/competence_headers_aa.txt",
                 header = F, sep = "\t") %>% 
        dplyr::rename(header = V1) %>% 
        dplyr::mutate(gene    = str_extract(header, "(?<= )[^\"]+?(?= )"), # " " & " "
                      prot_id = str_extract(header, "^[^ ]+"), # begin & " "
                      gene_class = "Competence",
                      strain  = str_extract(header, "(?<=\\[)[^\\]]+(?=\\])") # "[" and "]"
        )
      ,
      by = c("temporary_sseqid" = "prot_id")
    ) %>% 
    dplyr::filter(file_name == "Streptococcus_pneumoniae_RMD131_contigs_from_YM")
  ,
  by = c("file_name", "qseqid", "prot_id" = "temporary_sseqid"),
  relationship = "many-to-many"
) %>% 
  dplyr::arrange(dplyr::desc(aa_pident)) %>% 
  dplyr::distinct(file_name, temporary_sseqid, .keep_all = T) %>% 
  dplyr::select(-contains(".y")) %>% 
  dplyr::mutate(
    gene_present = case_when(
      nt_pident >= 80 ~ "present",
      TRUE ~ "absent"
    ),
    aa_lengthdiff = abs(abs(aa_sstart-aa_send) - abs(aa_qstart-aa_qend))
    # protein_function = case_when(
    #   aa_pident >= 95 & aa_gapopen == 0 ~ "functional",
    #   aa_pident >= 90 & aa_gapopen <= 5 & aa_mismatch <= 30 ~ "variant",
    #   aa_pident <= 85 & aa_mismatch >= 25 ~ "possibly defective",
    #   TRUE ~ "possibly defective"
    # )
  ) %>% 
  # view() %>%
  glimpse()

com_blastp <- read.table("outputs/result_blastp_comptence_additional_analysis/Streptococcus_pneumoniae_RMD131_contigs_from_YM_setA_blastp_results.txt",
                         header = F, sep = "\t") %>% 
  stats::setNames(c("qseqid", "temporary_sseqid",
                    "aap_pident", "aap_length",
                    # "aap_mismatch",
                    # "aap_gapopen", 
                    "aap_qstart", "aap_qend",
                    "aap_sstart", "aap_send",
                    "aap_evalue", "aap_bitscore", "qlen")) %>% 
  dplyr::left_join(
    read.table("inputs/prepare_competence_genes/competence_headers_aa.txt",
               header = F, sep = "\t") %>% 
      dplyr::rename(header = V1) %>% 
      dplyr::mutate(gene    = str_extract(header, "(?<= )[^\"]+?(?= )"), # " " & " "
                    prot_id = str_extract(header, "^[^ ]+"), # begin & " "
                    species  = str_extract(header, "(?<=\\[)[^\\]]+(?=\\])") # "[" and "]"
      )
    ,
    by = c("temporary_sseqid" = "prot_id")
  ) %>%
  dplyr::mutate(
    # aap_sdiff = abs(aap_sstart-aap_send),
    # aap_qdiff = abs(aap_qstart-aap_qend),
    aap_lengthdiff = abs(abs(aap_sstart-aap_send) - abs(aap_qstart-aap_qend)),
    alignment_length = abs(aap_qstart-aap_qend) + 1,
    query_coverage = alignment_length / qlen
  ) %>%  
  dplyr::filter(aap_evalue <= 1e-3,
                aap_lengthdiff <= 70,
                aap_bitscore >= 50,
                aap_pident >= 70,
                query_coverage >= 0.7
  ) %>%
  dplyr::group_by(gene, temporary_sseqid, species, query_coverage) %>%
  dplyr::summarise(
    hits = n(),
    mean_coverage = mean(query_coverage)*100,
    min_coverage = min(query_coverage)*100,
    max_coverage = max(query_coverage)*100,
    
    mean_identity = mean(aap_pident),
    min_identity = min(aap_pident),
    max_identity = max(aap_pident),
    
    mean_bitscore = mean(aap_bitscore),
    min_bitscore = min(aap_bitscore),
    max_bitscore = max(aap_bitscore),
    
    mean_evalue = mean(aap_evalue),
    min_evalue = min(aap_evalue),
    max_evalue = max(aap_evalue),
  ) %>%
  dplyr::arrange(desc(hits)) %>%
  dplyr::mutate(
    coverage = ifelse(
      hits > 1, 
      paste0(round(mean_coverage, 3), "% (", min_coverage, "%-", max_coverage, "%)"), 
      paste0(round(mean_coverage, 3), "%")
    ),
    avg_identity = ifelse(
      hits > 1, 
      paste0(round(mean_identity, 3), "% (", min_identity, "%-", max_identity, "%)"), 
      paste0(round(mean_identity, 3), "%")
    ),
    avg_evalue = paste0(round(mean_evalue, 3), " (", min_evalue, "-", max_evalue, ")")
    ,
    bitscore = ifelse(
      hits > 1, 
      paste0(round(mean_bitscore, 3), " (", min_bitscore, "-", max_bitscore, ")"), 
      paste0(round(mean_bitscore, 3), "")
    )
  ) %>% 
  dplyr::select(gene, temporary_sseqid, species, hits,
                coverage,
                avg_identity,
                # avg_evalue,
                bitscore) %>% 
  view() %>%
  glimpse()

report_blastp_com <- com_blastp %>% 
  dplyr::mutate(
    gene_class = "Competence"
  ) %>% 
  dplyr::select(gene_class, gene, species,
                hits, avg_identity, max_bitscore
  ) %>% 
  # view() %>%
  glimpse()



test <- com_combined_results %>% 
  select(nt_pident, nt_mismatch, nt_lengthdiff,
         aa_pident, aa_mismatch, aa_lengthdiff,
         aa_gapopen, gene.x, strain.x,
         gene_present, #protein_function
         ) %>% 
  view() %>% 
  glimpse()

test_gene_distinction <- com_combined_results %>% 
  dplyr::arrange(dplyr::desc(aa_pident)) %>% 
  dplyr::distinct(gene.x, .keep_all = T) %>% 
  view()
# comCDE are closely related to R6. Interesting.


# re-blast transposon genes
transposon_results <- read.table("outputs/result_blast_transposon/blastn_tabular_nt.txt",
                                   header = F, sep = "\t") %>% 
  stats::setNames(c("file_name", "qseqid", "temporary_sseqid",
                    "nt_pident", "nt_length", "nt_mismatch",
                    "nt_gapopen", "nt_qstart", "nt_qend",
                    "nt_sstart", "nt_send",
                    "nt_evalue", "nt_bitscore")) %>% 
  dplyr::filter(
    file_name == "Streptococcus_pneumoniae_RMD131_contigs_from_YM",
    !str_detect(temporary_sseqid, "Tn5253")
    ) %>% 
  dplyr::mutate(
    nt_lengthdiff = abs(abs(nt_sstart-nt_send) - abs(nt_qstart-nt_qend))
    ) %>% 
  dplyr::select(
    temporary_sseqid,
    nt_qstart,
    nt_qend,
    nt_lengthdiff
  ) %>% 
  view() %>% 
  glimpse()






