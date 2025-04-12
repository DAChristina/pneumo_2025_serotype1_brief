library(tidyverse)
# Genome analyses from Prokka
prokka_folder = "outputs/result_prokka_combined_pathogenWatch_monocle/"

# Genome (*.fna)
genome <- Biostrings::readDNAStringSet(paste0(prokka_folder, "Streptococcus_pneumoniae_RMD131.fna"))
cat("Genome width\n")
print(width(genome))

# Predicted genes (*.ffn)
genes <- Biostrings::readDNAStringSet(paste0(prokka_folder, "Streptococcus_pneumoniae_RMD131.ffn"))
gene_lengths <- width(genes)
genes_plot <- ggplot(data = data.frame(gene_lengths), aes(x = gene_lengths)) +
  geom_histogram(binwidth = 50, fill = "steelblue", color = "black") + 
  labs(title = "Gene Length Distribution",
       x = "Gene Length (bp)",
       y = "Count") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()

# Predicted proteins (*.faa)
proteins <- Biostrings::readAAStringSet(paste0(prokka_folder, "Streptococcus_pneumoniae_RMD131.faa"))
protein_lengths <- width(proteins)
prot_plot <- ggplot(data = data.frame(protein_lengths), aes(x = protein_lengths)) +
  geom_histogram(binwidth = 50, fill = "steelblue", color = "black") + 
  labs(title = "Protein Length Distribution",
       x = "Protein Length (aa)",
       y = "Count") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()

# Functional annotation (*.tsv)
annotations <- read.delim(paste0(prokka_folder, "Streptococcus_pneumoniae_RMD131.tsv"), header = TRUE, sep = "\t")
# head(annotations)
fun_plot <- annotations %>%
  count(product) %>%
  top_n(10) %>%
  ggplot(aes(x = reorder(product, n), y = n)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = n), hjust = -0.1) +  # Add text labels slightly outside the bars
  coord_flip() +
  theme_bw() +
  labs(title = "Top 10\nGene Functions", x = "Function", y = "Count")

p_genes <- cowplot::plot_grid(NULL, plot(genes_plot),
                              ncol = 1,
                              rel_heights = c(0.1, 2),
                              labels = c("", "A"))

p_proteins <- cowplot::plot_grid(prot_plot, fun_plot,
                                 ncol = 1,
                                 labels = c("B", "C"))

p_gene_prot <- cowplot::plot_grid(p_genes, p_proteins,
                                  ncol = 2)

p_gene_prot

# mixed infection, GPSC and MLST profile
kity_folder = "outputs/result_pneumokity/fastq_mix/"
mlst_folder = "outputs/result_mlst/"

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
                      species = str_extract(header, "(?<=\\] \\[)[^\\]]+(?=\\])") # "] [" and "]"
        )
      ,
      by = c("temporary_sseqid" = "gene_id")
    ) %>% 
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
      nt_pident >= 80 ~ "present",
      TRUE ~ "absent"
    ),
    protein_function = case_when(
      aa_pident >= 95 & aa_gapopen == 0 ~ "functional",
      aa_pident >= 90 & aa_gapopen <= 5 & aa_mismatch <= 30 ~ "variant",
      aa_pident <= 85 & aa_mismatch >= 25 ~ "possibly defective",
      TRUE ~ "possibly defective"
    ),
    bacwgs_check = case_when(
      gene.x %in% c("cbpD", "cps4A", "cps4B", "cps4D", "hysA", "lytA", "lytC", 
                    "nanA", "nanB", "pavA", "pce/cbpE", "pfbA", "ply", "psaA") ~ "detected",
      TRUE ~ "no"
    )
  ) %>% 
  # view() %>% 
  glimpse()

test <- vir_combined_results %>% 
  select(nt_pident, aa_pident, aa_mismatch, aa_gapopen, gene.x, species.x,
         gene_present, protein_function, bacwgs_check) %>% 
  view()


# competence genes profile












