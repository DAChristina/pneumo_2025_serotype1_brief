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
