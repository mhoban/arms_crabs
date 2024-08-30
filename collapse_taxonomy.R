# load taxonomy table and transform it
taxonomy <- read_csv(here("data","taxonomy.csv"))
taxonomy_collapsed <- taxonomy %>%
  mutate(
    species = case_when(
      # fix two misspellings
      species == "Charybdis (Charybdis) hawaiensis" ~ "Charybdis Charybdis hawaiensis",
      species == "Catoptrus inequalis" ~ "Catoptrus inaequalis",
      # make morphospecies into just "genus/family sp."
      str_detect(species,"sp[0-9]+") ~ str_c(coalesce(genus,family)," sp."),
      .default = species
    ),
    # reset known taxon level
    taxon_level = case_when(
      species == str_glue("{genus} sp.") ~ "genus",
      species == str_glue("{family} sp.") ~ "family",
      .default = taxon_level
    )
  ) %>%
  # assign them new otu IDs
  group_by(species) %>%
  mutate( old_otu = otu,  otu = str_replace_all(str_replace(species,fixed(" sp."),"")," ","_")) %>%
  ungroup() %>%
  select(-ends_with("otu"),ends_with("otu")) %>%
  arrange(family,genus,species)

# make a lookup table
tax_lookup <- taxonomy_collapsed %>% select(ends_with("otu"))


# now remap the community matrix
crabs <- read_csv(here("data","crabs.csv"))

crabs_collapsed <- crabs %>%
  # pivot long and get the new otu IDs
  pivot_longer(-unit,names_to = "otu",values_to = "count") %>%
  inner_join(tax_lookup,by=c("otu" = "old_otu")) %>%
  select(-otu,unit,otu=otu.y,count) %>%
  # pivot back to wide and sum values by new otu ID
  pivot_wider(names_from = "otu", values_from = "count", values_fn = sum)

# make sure the taxonomy table is unique
taxonomy_collapsed <- taxonomy_collapsed %>% select(otu,taxon_level:species) %>%
  distinct(otu,taxon_level,family,genus,species)

# save the new collapsed tables
dir_create(here("data","generated"))
write_csv(taxonomy_collapsed,here("data","generated","taxonomy_collapsed.csv"))
write_csv(crabs_collapsed,here("data","generated","crabs_collapsed.csv"))
