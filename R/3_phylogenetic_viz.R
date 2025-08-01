# generate AllTarget_isolates.list for SKA2 & gubbins
library(tidyverse)
library(ggtree)
library(ggtreeExtra)

metadata <- read.csv("inputs/prepare_gubbins_tree_no_outgroups/metadata_GPSC31_outgroup.csv")

target <- data.frame(ID = metadata$Lane_id) %>% 
  dplyr::mutate(
    across(where(is.character), ~ str_replace_all(., "#", "_")),
    location = paste0("/home/ron/pneumo_2025_serotype1_brief_report/inputs/prepare_gubbins_tree_no_outgroups/", ID, ".fa")
  ) %>% 
  glimpse()

# write.table(target,
#             "inputs/prepare_gubbins_tree_no_outgroups/AllTarget_isolates.list",
#             row.names = F, col.names = F, quote = F,
#             sep = "\t")

################################################################################
# load gubbins tree
tre <- ape::read.tree("outputs/result_tree_no_outgroup/GPSC2.node_labelled.final_tree.tre")
tre$tip.label <- gsub("^Streptococcus_pneumoniae_", "", tre$tip.label)
tre$tip.label <- gsub(".contigs_velvet$", "", tre$tip.label)

ape::write.tree(tre, "outputs/result_tree_no_outgroup/GPSC2.node_labelled.final_tree_shortened_ID.tre")

# test node
ggtree(tre) + 
  geom_tiplab(size = 2) +
  geom_label2(aes(subset=!isTip, label=node), size=2, color="darkred", alpha=0.5)

# analyse weird subtree:
subtre <- ape::extract.clade(tre, node = 135)
ggtree(subtre) + 
  geom_tiplab(size = 2) +
  geom_label2(aes(subset=!isTip, label=node), size=3, color="darkred", alpha=0.5)

# arrange metadata based on tre$tip.label
metadata_final <- dplyr::left_join(
  data.frame(tre$tip.label)
  ,
  read.csv("outputs/result_tree_no_outgroup/metadata_all_GPSC2_monocle.csv") %>% 
    dplyr::mutate(
      across(where(is.character), ~ str_replace_all(., "#", "_")),
      across(where(is.character), ~ str_replace_all(., "Streptococcus_pneumoniae_", ""))
    ) 
  ,
  by = c("tre.tip.label" = "Lane_id")
) %>% 
  dplyr::mutate(
    Continent = case_when(
      Country == "NIGERIA" ~ "AFRICA",
      T ~ Continent
    ),
    Lane_id = tre.tip.label
  ) %>% 
  # factor correction
  dplyr::mutate(
    label_Indo = ifelse(Country == "INDONESIA", "Indonesia", "Others"),
    Continent = factor(Continent,
                       levels = c("AFRICA", "ASIA", "EUROPE")),
    Country = factor(Country,
                     levels = c("THE GAMBIA", "MALAWI", "NIGERIA",
                                "CAMBODIA", "INDIA", "INDONESIA", "NEPAL",
                                "NORWAY")),
    In_silico_ST_regrouped = case_when(
      In_silico_ST == 9529 | In_silico_ST == 14507 | In_silico_ST == 17067 ~ "OTHERS",
      T ~ as.character(In_silico_ST)
    ),
    In_silico_ST_regrouped = factor(In_silico_ST_regrouped,
                          levels <- c("217", "303",
                                      "3081", "5002", "5012", "5672", "8158",
                                      "12197",
                                      "OTHERS")),
    null_data = NA
                          
    
  ) %>% 
  glimpse()
rownames(metadata_final) <- tre$tip.label

write.csv(metadata_final,
          "outputs/result_tree_no_outgroup/metadata_final.csv")

# basic
show_pp <- ggtree(tre,
                  layout = "fan",
                  open.angle=30,
                  size=0.75
                  ) +
  geom_hilight(node=135, fill="pink", alpha=0.5)
show_pp

show_pp2 <- ggtree(tre,
                   layout = "fan",
                   open.angle=30,
                   size=0.75,
                   aes(colour = label_Indo)
) %<+% 
  metadata_final +
  theme(legend.position='none') +
  scale_colour_manual(
    values=c("red","grey40")
  ) + 
  geom_hilight(node=135, fill="pink", alpha=0.5) + 
  ggnewscale::new_scale_colour()
show_pp2

# trial grouping Indo as clade for the clade group
temp_indoLabel <- metadata_final %>% select(c("Lane_id", "label_Indo"))
temp_indoLabel <- aggregate(.~label_Indo, temp_indoLabel, FUN=paste, collapse=",")
label_Indos <- lapply(temp_indoLabel$Lane_id, function(x){unlist(strsplit(x,split=","))})
names(label_Indos) <- temp_indoLabel$label_Indo

tr <- groupOTU(tre, label_Indos, "Clade")
Clade <- NULL
show_pp3 <- ggtree(tr=tr, layout = "fan",
                       open.angle=30,
                       size=0.75,
                       aes(colour = Clade)
                   ) +
  scale_colour_manual(
    values=c("red","grey40"),
    guide = "none"
  ) +
  theme(legend.position = "none") +
  geom_hilight(node=135, fill="pink", alpha=0.5) + 
  ggnewscale::new_scale_colour()
show_pp3

tree_gen_country_tree <- show_pp3 %<+%
  metadata_final +
  # continent
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=metadata_final$Continent),
    width=20,
    offset=0.1
  ) +
  scale_fill_viridis_d(
    name = "Continent",
    option = "C",
    direction = -1,
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3,
                         ncol = 5, order = 1)
  ) +
  theme(
    legend.title=element_text(size=12), 
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm")
  ) +
  # country
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=metadata_final$Country),
    width=20,
    offset=0.1
  ) +
  scale_fill_viridis_d(
    name = "Country",
    option = "C",
    direction = -1,
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3,
                         ncol = 3, order = 2)
  ) +
  theme(
    legend.title=element_text(size=12),
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm")
  ) +
  # sequence type
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=metadata_final$In_silico_ST_regrouped),
    width=20,
    offset=0.1
  ) +
  scale_fill_viridis_d(
    name = "Sequence type",
    option = "C",
    direction = -1,
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3,
                         ncol = 4, order = 3)
  ) +
  theme(
    legend.title=element_text(size=12),
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm")
  ) +
  # null
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=metadata_final$null_data),
    width=20,
    offset=0.1
  ) +
  scale_fill_viridis_d(
    name = "null_data",
    option = "C",
    direction = -1,
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3,
                         ncol = 4, order = 3)
  ) +
  theme(
    legend.title=element_text(size=12),
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm")
  ) 
tree_gen_country_tree

# AMR chosen for transposon-mediated AMR COT, TET, cat, ermB, mefA
filtered_df <- metadata_final %>% 
  dplyr::select(
    contains("autocolour"),
    -PBP1A_2B_2X__autocolour,
    ermB, mefA, cat
  ) %>% 
  dplyr::rename(
    Tetracycline = Tet__autocolour,
    Sulfamethoxazole = folP__autocolour,
  ) %>% 
  dplyr::mutate(
    Tetracycline = case_when(
      Tetracycline == "TETM" ~ "Tet(M)",
      Tetracycline == "NEG" | Tetracycline == "" | is.na(Tetracycline) ~ "Not found",
      TRUE ~ Tetracycline
    ),
    Sulfamethoxazole = case_when(
      Sulfamethoxazole == "FOLP_169_INS" | Sulfamethoxazole == "FOLP_AE007317 INS AT 169" ~ "FolP 169 insertion",
      Sulfamethoxazole == "FOLP_178_INS" ~ "FolP 178 insertion",
    ),
    Erythromycin = case_when(
      ermB == "POS" ~ "ermB",
      mefA == "POS" ~ "mefA",
      ermB == "POS" & mefA == "POS" ~ "ermB & mefA",
      TRUE ~ "Not found"
    ),
    Chloramphenicol = case_when(
      cat == "POS" ~ "cat",
      TRUE ~ "Not found"
    )
  ) %>% 
  dplyr::select(Chloramphenicol,
                Erythromycin,
                Sulfamethoxazole,
                Tetracycline,
                -c(ermB, mefA, cat)) %>% 
  glimpse()

all_labels <- unique(unlist(filtered_df))
manual <- c("Not found" = "palegoldenrod",
            # "NEG" = "palegoldenrod",
            "NA" = "white")
others <- setdiff(all_labels, names(manual))

auto_col <- scales::hue_pal()(length(others))
names(auto_col) <- others

final_col <- c(manual, auto_col)
factor_levels <- c("FolP 169 insertion", 
                   "FolP 178 insertion", 
                   "Tet(M)", 
                   "cat",
                   "mefA",
                   "Not found", 
                   "NA")

filtered_df <- filtered_df %>%
  dplyr::mutate(across(everything(), ~factor(.x, levels = factor_levels)))

# tree_gen_country_tree %<+%
#   filtered_df +
# ggnewscale::new_scale_fill() +
png(file = "pictures/tree.png",
    width = 25, height = 12, unit = "cm", res = 600)
ggtree::gheatmap(tree_gen_country_tree, filtered_df,
                 offset=70, width=0.5, font.size=3, 
                 colnames_angle=-60, hjust=0,
                 legend_title = "Antimicrobial resistance") +
  scale_fill_manual(
    values = final_col,
    na.translate = FALSE,
    name = "Antimicrobial resistance",
    guide = guide_legend(ncol = 2)
  )
dev.off()
