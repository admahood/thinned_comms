# species data

source("R/a_tc_data_prep.R")
usda_list <- read_csv("data/usda/usda_plantlist.csv") %>%
  janitor::clean_names() %>%
  filter(is.na(synonym_symbol)) %>%
  dplyr::select(-synonym_symbol)

invasives <- read_csv("data/usda/invasive.csv", skip =3) %>%
  janitor::clean_names() %>%
  dplyr::select(symbol = accepted_symbol) %>%
  mutate(origin = "introduced")

code_types <- comm_raw %>%
  dplyr::select(CodeType, symbol = SpeciesCode) %>%
  unique()

prevalence <- comm
prevalence[prevalence>0] <- 1
prv <- colSums(prevalence)
prevalence <- tibble(symbol = names(prv), prevalence = prv) %>%
  arrange(prevalence)
prevalence_summary <-
  prevalence %>%
  group_by(prevalence) %>%
  summarise(n = n())

sp_df_raw <- tibble(symbol = colnames(comm)) %>%
  left_join(usda_list) %>%
  tidyr::separate(scientific_name_with_author, 
                  into = c("genus", "species", "author"), remove = F,
                  extra = "merge") %>%
  mutate(scientific_name = str_c(genus, " ", species))
apply(sp_df_raw, 2, function(x)sum(is.na(x)))

sp_df_clean <- sp_df_raw %>%
  dplyr::select(-family, -common_name, -species, -author, -genus) %>%
  left_join(code_types)

print(sp_df_clean %>% arrange(family, genus), n=50)
apply(sp_df_clean, 2, function(x)sum(is.na(x)))





