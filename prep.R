library(tidyverse)
library(here)


all_arms <- read_csv(here("data","all_data.csv"),col_types=cols())

islands <- tribble(
  ~island,~island_name,
  "OAH","Oahu",
  "FFS","French Frigate Shoals",
  "HAW","Hawaii Island",
  "KAU","Kauai",
  "KUR","Kure Atoll",
  "LIS","Lisianski",
  "MAI","Maui",
  "PHR","Pearl and Hermes Reef"
)

metadata <- all_arms %>%
  distinct(unit,.keep_all = T) %>%
  select(unit,region,island,deployment_year,recovery_year,depth,lat,lon) %>%
  inner_join(islands,by="island") %>%
  select(-island) %>%
  rename(island=island_name) %>% 
  select(unit,island,everything()) 

write_csv(metadata,here("data","metadata.csv"))
  
crabs_long <- read_csv(here("data","crabs_long.csv"),col_types = cols())
crabs <- all_arms %>%
  group_by(unit,species,family,genus,taxon_level,region,island,depth,lat,lon) %>%
  summarise(count = sum(count)) %>%
  ungroup() %>%
  select(unit,species,count) %>%
  arrange(unit,species) %>%
  pivot_wider(names_from="species",values_from="count",values_fill = 0)

write_csv(crabs,here("data","crabs.csv"))


taxonomy <- all_arms %>%
  select(species,family,genus,taxon_level) %>%
  distinct(species,.keep_all = T) %>%
  rename(otu = species) %>%
  mutate(species = str_replace_all(otu,"_"," ")) %>%
  select(otu,taxon_level,family,genus,species)
  
write_csv(taxonomy,here("data","taxonomy.csv"))


metadata <- read_csv(here("data","all_data.csv"))
