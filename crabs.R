# required libraries
# all of these can be directly installed except for phyloseq, which is a bioconductor package
# phyloseq isn't directly necessary, but it's what I used for all my eDNA stuff and it conveniently
# wraps the datasets, so I'm using it here
library(tidyverse)
# to install phyloseq in regular R, do this:

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("phyloseq")
library(phyloseq) 
library(vegan)
library(ggrepel)
library(here)
library(worrms)

# to check model assumptions
library(performance)

# for linear mixed effects models
library(lme4)

# set the random seed
set.seed(31337)

run_things <- FALSE

# functions ---------------------------------------------------------------

# plot an ordination with various options
plot_ord <- function(ps,ord,arrows=TRUE,sig=FALSE,by="term",perm=999,alpha=0.05,scale=1,labelsize=5,pointsize=3,...) {
  
  # map classnames to axis names
  axes <- c(
    "dbrda" = "dbRDA",
    "dbRDA" = "dbRDA",
    "capscale" = "CAP",
    "CAP" = "CAP",
    "rda" = "RDA",
    "RDA" = "RDA" ,
    "cca" = "CCA",
    "CCA" = "CCA"
  )
  # figure out axis names and b0rk out if they don't exist
  axis <- axes[class(ord)[1]] 
  if (is.na(axis) | is.null(axis)) {
    stop("invalid analysis type")
  }
  
  # get the model terms
  terms <- attr(ord$terms,"term.labels")
  if (arrows & sig) {
    # get *significant* model terms, if requested
    terms <- anova(ord,by=by,permutations=perm) %>%
      as_tibble(rownames="term") %>%
      filter(term != "Residual",`Pr(>F)` < alpha) %>%
      pull(term)
  }
  
  # do the initial plot with whatever options were passed in
  plotz <- plot_ordination(ps,ord,...) + 
    geom_point(size=pointsize)
  
  # plot loading arrows, if requested
  if (arrows) {
    # regex to match term names
    pattern <- str_c("^(",str_c(terms,collapse="|"),")")
    # make a tibble of arrow endpoints, optionally scaled
    arrowdf  <- vegan::scores(ord, display = "bp") %>%
      as_tibble(rownames = "labels") %>%
      # this is where we scale the arrow endpoints
      mutate(across(where(is.numeric),~.x*scale))
    if (sig) {
      # if we only want significant terms, this is 
      # where we filter for that
      arrowdf <- arrowdf %>%
        filter(str_detect(labels,pattern))
    }
    
    # Define the arrow aesthetic mappings
    # map for arrowheads
    arrow_map <- aes(xend = .data[[str_c(axis,"1")]], # this is where we use our mapped axis name
                     yend = .data[[str_c(axis,"2")]], # e.g. RDA1, CCA1, etc.
                     x = 0, 
                     y = 0, 
                     shape = NULL, 
                     color = NULL)
    # map for arrow labels
    label_map <- aes(x = 1.3 * .data[[str_c(axis,"1")]], 
                     y = 1.3 * .data[[str_c(axis,"2")]], 
                     shape = NULL, 
                     color = NULL, 
                     label = labels)
    # make a little arrowhead
    arrowhead = arrow(length = unit(0.02, "npc"))
    # add the arrows to the plot
    plotz <- plotz + 
      geom_segment(
        mapping = arrow_map, 
        size = .5, 
        data = arrowdf, 
        color = "black", 
        arrow = arrowhead
      ) + 
      geom_text_repel(
        mapping = label_map, 
        size = labelsize,  
        data = arrowdf, 
        show.legend = FALSE,
        color='black'
      ) 
  }
  return(plotz)
}

# do a pairwise PERMANOVA analysis by some factor
pairwise_adonis <- function(comm, factors, permutations = 1000, correction = "fdr", method = "bray") {
  # get possible pairwise factor combinations
  factor_combos <- combn(unique(factors), 2)
  # map through factor combinations and run model for each pair,
  # dumping the output into a tibble
  model_output <- map_dfr(array_branch(factor_combos,2), ~{
    fact <- factors[factors %in% .x]
    
    # get our factor specific community or distance matrix
    if (inherits(comm,'dist')) {
      dd <- as.matrix(comm)[factors %in% .x, factors %in% .x]
    } else {
      comm <- as(comm,'matrix')
      dd <- vegdist(comm[factors %in% .x,], method = method)
    }
    
    # run the permanova model and return the results in a list
    # to make the data frame rows
    model <- adonis2(dd ~ fact, permutations = permutations) %>%
      as_tibble()
    df <- model$Df[1]
    ss <- model$SumOfSqs[1]
    pseudo_f <- model$F[1]
    p_val <- model$`Pr(>F)`[1]
    
    list(
      "left" = .x[1],
      "right" = .x[2],
      "df" = df,
      "ss" = ss,
      "pseudo_f" = pseudo_f,
      "p_val" = p_val
    )
  })
  
  # return the results with adjusted p-values
  model_output %>%
    mutate(
      p_adj = p.adjust(p_val, method = correction)
    )
}

# a function to standardize a phyloseq object, basically 'decostand' for ps objects
# with a few extras thrown in
ps_standardize <- function(ps, method="total", ...) {
  # supported methods
  vegan_methods <- c("total", "max", "frequency", "normalize", "range", "rank", "standardize", "pa", "chi.square", "hellinger", "log")
  other_methods <- c("wisconsin","sqrt")
  
  # how to interpret the data
  trows <- phyloseq::taxa_are_rows(ps)
  if(trows == TRUE) { 
    marg <- 2 
  } else {
    marg <- 1 
  }
  
  # get the otu table matrix
  otus <- ps %>%
    otu_table() %>%
    as("matrix")
  
  # a special case where we can call "total" "relative" if we want to
  if (method == "relative") {
    method <- "total"
  }
  
  if (method %in% vegan_methods) {
    # if it's supported by decostand, do it
    otus_std <- decostand(otus,method=method,MARGIN=marg,...)
    if (method == "chi.square") {
      # transpose if necessary
      otus_std <- t(otus_std)
    }
  } else if (method %in% other_methods) {
    # otherwise do the specific things
    if (method == "wisconsin") {
      if (trows) {
        otus_std <- wisconsin(t(otus))
        otus_std <- t(otus_std)
      } else {
        otus_std <- wisconsin(otus)
      }
    } else if (method == "sqrt") {
      otus_std <- sqrt(otus)
    }
  } else {
    otus_std <- NULL
  }
  
  if (!is.null(otus_std)) {
    # reassign the transformed otu tables if we succeeded in transforming it
    otu_table(ps) <- otu_table(otus_std,taxa_are_rows = trows)
  }
  return(ps) 
}

# an easy way to load a phyloseq object from otu matrix, sample data, and taxonomy
load_ps <- function(ott, sample_col="sample", tt, otu_col="otu", ssd, taxa_are_rows = FALSE) {
  # read the otu table, convert the appropriate column to rownames, and make it a matrix
  otus <- read_csv(ott,col_types = cols()) %>% 
    column_to_rownames(sample_col) %>%
    as("matrix")
  # read the taxonomy table and do the same thing
  tax <- read_csv(tt,col_types = cols()) %>%
    column_to_rownames(otu_col) %>%
    as("matrix")
  # likewise the sample data
  sd <- read_csv(ssd,col_types = cols()) %>%
    column_to_rownames(sample_col)
  
  return(
    # smash it into a phyloseq object
    phyloseq(
      otu_table(otus,taxa_are_rows = taxa_are_rows),
      tax_table(tax),
      sample_data(sd)
    )
  )
}

# a quick wrapper to get a tibble of a phyloseq's sample data
sample_tibble <- function(ps,sample_col="sample") {
  ps %>%
    sample_data() %>%
    as.data.frame() %>%
    as_tibble(rownames=sample_col)
}

# a wrapper to grab a tibble of the taxonomy
taxa_tibble <- function(ps,otu_col="otu") {
  ps %>%
    tax_table() %>%
    as("matrix") %>%
    as_tibble(rownames=otu_col)
}


# data prep ---------------------------------------------------------------

setup_crabs <- function() {
  cc <- list()
  # load our main phyloseq object and hellinger-transform it
  cc$crabs_untransformed <- load_ps(
    here("data","crabs.csv"),"unit",
    here("data","taxonomy.csv"),"otu",
    here("data","metadata.csv")
  ) 
  
  cc$crabs <- cc$crabs_untransformed %>%
    ps_standardize("hellinger")
    # this is the shallow subset
  cc$crabs_shallow <- cc$crabs %>%
    subset_samples(region != "mce") 
    # we're gonna mess with the sample data a bit
  # first, get rid of rows that have NAs
  cc$crab_data_shallow <- cc$crabs_shallow %>%
    sample_tibble(sample_col = "unit") %>%
    drop_na()
    # make sure we drop samples that had NAs in their sample data
  cc$crabs_shallow <- prune_samples(sample_names(cc$crabs_shallow) %in% cc$crab_data_shallow$unit,cc$crabs_shallow)
    # now we select the columns containing the environmental variables
  # we care about and scale the numeric ones to unit variance
  cc$crab_data_shallow <- cc$crabs_shallow %>%
    sample_tibble(sample_col = "unit") %>%
    select(unit,region,island,depth,chl,sst,slope,coral_cover,closest_island,larval_connectivity,human_impact) %>%
    mutate(across(where(is.numeric),~as.numeric(scale(.x)))) %>%
    column_to_rownames("unit") 
    # reassociate the new sample data
  sample_data(cc$crabs_shallow) <- cc$crab_data_shallow
    # precalculate the bray-curtis distance for all crabs
  cc$crab_dist <- distance(cc$crabs,"bray")
  # pull out the otu table
  cc$crab_otus <- otu_table(cc$crabs)
  # get a nice tibble of the sample data
  cc$crab_data <- sample_tibble(cc$crabs)
    # likewise for the shallow subset
  cc$crab_dist_shallow <- distance(cc$crabs_shallow,"bray")
  cc$crab_otus_shallow <- otu_table(cc$crabs_shallow)
  cc$crab_data_shallow <- sample_tibble(cc$crabs_shallow)
  
  return(cc)
}

sample_summary <- function(ps) {
  ss <- list() 
  
  # basic numbers
  ss$total <- ps %>%
    sample_sums() %>%
    sum()
  
  ss$abundance_range <- ps %>%
    sample_sums() %>%
    range()
  
  ss$abundance_mean <- ps %>%
    sample_sums() %>%
    mean()
  
  ss$abundance_sd <- ps %>%
    sample_sums() %>%
    sd()
  
  # get 5 most common species
  top5 <- ps %>%
    taxa_sums() %>%
    sort() %>%
    rev() %>%
    head(5) 
  
  otus <- top5 %>%
    names() %>%
    str_replace_all("_"," ") 
  
  ss$top5_n <- top5
  
  
  to_keep <- ps %>%
    taxa_tibble() %>%
    filter(species %in% otus) %>%
    pull(otu)
  
  # ss$top5_units <- ps %>% 
  #   subset_taxa(species %in% otus) %>%
  ss$top5_units <- prune_taxa(to_keep,ps) %>%
    ps_standardize("pa") %>%
    taxa_sums() %>% 
    sort() %>%
    rev()
  
  # make it the same order as the top 5 species
  ss$top5_units <- ss$top5_units[match(names(ss$top5_units),names(top5))]
  
  ss$top5 <- wm_records_taxamatch(otus) %>%
    map2(otus,~{
      if (nrow(.x) > 0) {
        .x <- .x %>% 
          filter(rank == "Species") %>%
          select(species = scientificname, authority) %>%
          slice(1)  %>%
          mutate(found = TRUE)
      } 
      if (nrow(.x) == 0) {
        .x <- list(species = .y, authority = "", found = FALSE)
      }
      return(as.list(.x))
    })
  
  
  ts <- ps %>% taxa_sums()
  ss$singleton_count <- sum(ts == 1)
  
  ss$singleton_species <-  ts[ts == 1] %>%
    names() %>%
    sort()
  
  tt <- ps %>%
    taxa_tibble() %>%
    mutate(thing = "everything")
  
  ss$types_breakdown <- tt %>%
    group_by(taxon_level) %>%
    summarise(families = n_distinct(family), genera = n_distinct(genus), species = n_distinct(species)) %>%
    select(taxon_level,species) %>%
    deframe() %>%
    as.list()
  
  ss$taxonomy_breakdown <- tt %>%
    group_by(thing) %>%
    summarise(families = n_distinct(family), genera = n_distinct(genus), species = n_distinct(species)) %>%
    select(-thing) %>%
    as.list()
  
  return(ss)
}

alpha_diversity <- function(ps,measures=c("Observed","Simpson")) {
  
  ad <- list()
  
  # skip fisher because it's broken
  richness <- ps %>%
    estimate_richness(measures = measures) %>%
    as_tibble(rownames="sample") %>%
    mutate(sample=str_replace(sample,"^X",""))
  
  cd <- ps %>%
    sample_tibble()
  
  all_richness <- richness %>%
    inner_join(cd,by="sample") %>%
    mutate(island = fct_reorder(island,-lat))
  
  ad$richness_table <- all_richness
  ad$richness <- mean(all_richness$Observed)
  ad$richness_sd <- sd(all_richness$Observed)
  ad$simpson <- mean(all_richness$Simpson)
  ad$simpson_sd <- sd(all_richness$Simpson)
  
  ri <- all_richness %>%
    group_by(island_group) %>%
    summarise(richness = mean(Observed), richness_sd = sd(Observed), simpson = mean(Simpson), simpson_sd = sd(Simpson))
  ad$richness_main <- ri %>%
    filter(island_group == "main") %>%
    select(-island_group) %>%
    as.list()
  ad$richness_nwhi <- ri %>%
    filter(island_group == "northwest") %>%
    select(-island_group) %>%
    as.list()
  
  try(ad$t_r_island <- t.test(Observed ~ island_group, data=all_richness),silent=TRUE)
  try(ad$t_s_island <- t.test(Simpson ~ island_group, data=all_richness),silent=TRUE)
  
  # m <- lm(Observed ~ island_group, data=all_richness)
  # summary(m)
  # check_model(m)
  # 
  # m <- lm(Simpson ~ island_group, data=all_richness)
  # summary(m)
  # check_model(m)
  # 
  # kruskal.test(Simpson ~ island_group, data = all_richness)
  # wilcox.test(Simpson ~ island_group, data = all_richness)
  # 
  # t.test(Observed ~ shallow_deep, data=all_richness)
  # t.test(Simpson ~ shallow_deep, data=all_richness)
  # f <- lm(Simpson ~ 1 + I(shallow_deep == "deep"),data=all_richness)
  # summary(f)
  # confint(f)
  # 
  # # environmental variables and alpha diversity
  # 
  # scaled_richness <- all_richness %>%
  #   mutate(
  #     across(c(lat,depth,chl,sst,slope,coral_cover,closest_island,larval_connectivity,human_impact),~as.numeric(scale(.x)))
  #   )
  # # lme_richness <- lmer(Observed ~ lat + depth + chl + sst + slope + coral_cover + closest_island + larval_connectivity + human_impact + (1|region/island), data=scaled_richness)
  # lme_richness <- lmer(Observed ~ lat + depth + chl + sst + slope + coral_cover + closest_island + larval_connectivity + human_impact + (1|region) + (1|island), data=scaled_richness)
  return(ad)
}

# start here

if (run_things) {
# alpha diversity / richness ----------------------------------------------
# do the  dbRDA analysis ---------------------------------------------------

# establish our upper and lower bounds for forward model selection
# upper is all the terms, lower is none of the terms
db_upr <- dbrda(crab_otus_shallow ~ depth + chl + sst + slope + coral_cover + closest_island + larval_connectivity + human_impact,data=crab_data_shallow,distance="bray")
db_lwr <- dbrda(crab_otus_shallow  ~ 1,data=crab_data_shallow,distance="bray")

# do our forward selection to find the best model
# set trace to TRUE if you want to see a bunch of output
db_modl <- ordiR2step(db_lwr,db_upr,trace=FALSE)

# see marginal effects of model terms
anova(db_modl,by="margin")


vif.cca(db_upr)
vif.cca(db_modl)

# warning: the following *may* be sketchy:
# dbRDA doesn't automatically give you species scores, so let's assign them
# but the community data needs to be transformed to make sense
# see https://github.com/vegandevs/vegan/issues/254 for information on the transformation below
sppscores(db_modl) <- sqrt(decostand(crab_otus_shallow, "total")/2) 

# this was to grab the 10 most abundant species overall
top_species <- names(head(rev(sort(taxa_sums(crabs_shallow)))))
# but instead, we're just using the same species as the RDA from the initial manuscript
other_species <- c("Carupa_tenuipes","Dynomene_hispida","Percnon_abbreviatum","Xanthias_latifrons","Chlorodiella_laevissima","Garthiella_aberrans","Perinia_tumida")

# extract the species scores for the species we want
species_scores <- scores(db_modl,display="species") %>%
  as_tibble(rownames = "species") %>%
  filter(species %in% other_species)

# plot the ordination without species
dbrda_shallow_env <- plot_ord(crabs_shallow,db_modl,arrows=TRUE,sig=TRUE,scale=3,type="samples",color="region") + 
  stat_ellipse()
dbrda_shallow_env 

# plot the ordination with species
dbrda_shallow_env_spp <- plot_ord(crabs_shallow,db_modl,arrows=TRUE,sig=TRUE,scale=3,type="samples",color="region") + 
  stat_ellipse() +
  geom_point(aes(x=dbRDA1,y=dbRDA2),color="black",data=species_scores) +
  geom_text_repel(aes(x=dbRDA1,y=dbRDA2,label=species),color="black",data=species_scores) 
dbrda_shallow_env_spp

# save the ordinations
ggsave(plot=dbrda_shallow_env_spp,filename=here("..","manuscript","reanalysis figures","dbRDA_shallow_reduced_spp.pdf"),device=cairo_pdf,width=12,height=9,units="in")
ggsave(plot=dbrda_shallow_env,filename=here("..","manuscript","reanalysis figures","dbRDA_shallow_reduced.pdf"),device=cairo_pdf,width=12,height=9,units="in")



# do the CAP analysis for the complete dataset ----------------------------

# this turns out to be roughly the same thing as the dbRDA, but not as good
# the CAP that PRIMER does is different. So let's mostly ignore this section
# we can do this directly through phyloseq

cap_obj <- ordinate(
  crabs,
  method = "CAP",
  distance = "bray",
  binary = FALSE,
  formula = ~region + island_group + sample_group
  # formula = ~region + depth_zone + island_group
  # formula = ~island_group + region + depth_zone
  # formula = ~island_group + depth_zone + region
  # formula = ~depth_zone + island_group + region
  # formula = ~depth_zone + region + depth_zone
)

# see marginal effects of terms
anova(cap_obj,by="margin")

# plot it
cap_region_depth  <- plot_ordination(crabs,cap_obj,type="samples",color="sample_group",shape="depth_zone") +
  geom_point(size=3) + 
  stat_ellipse(aes(group=sample_group))
cap_region_depth 

# and save it
ggsave(plot=cap_region_depth,filename=here("..","manuscript","reanalysis figures","CAP_region_depth.pdf"),device=cairo_pdf,width=12,height=9,units="in")


# a couple little PCoA ordinations ----------------------------------------

# do and plot the pcoa for island region vs depth zone
pcoa <- ordinate(crabs_untransformed,"PCoA","bray")
pcoa_region_depth <- plot_ordination(crabs_untransformed,pcoa,color="sample_group",shape="depth_zone") + 
  geom_point(size=3) + 
  stat_ellipse(aes(group=sample_group),level=0.95) + 
  scale_y_reverse()
pcoa_region_depth 
# save the plot
ggsave(plot=pcoa_region_depth,filename=here("..","manuscript","reanalysis figures","PCoA_region_depth.pdf"),device=cairo_pdf,width=12,height=9,units="in")


# pcoa for just the shallow guys
pcoa <- ordinate(crabs_shallow,"PCoA","bray")
pcoa_shallow <- plot_ordination(crabs_shallow,pcoa,color="region") + 
  stat_ellipse(aes(group=region))
pcoa_shallow 

# do some PERMANOVA analyses ----------------------------------------------

# permanova tests for region and island, shallow crabs
anova_region <- adonis2(crab_dist_shallow ~ region,data=crab_data_shallow,permutations=999)
anova_island <- adonis2(crab_dist_shallow ~ island,data=crab_data_shallow,permutations=999)

# beta dispersion test for region, shallow crabs
pd <- betadisper(crab_dist_shallow,crab_data_shallow$region)
anova(pd,permutations=999)
# plot the confidence interval for the beta dispersion estimate
plot(TukeyHSD(pd))



# here's a pairwise PERMANOVA by region, for some reason
pairwise_adonis(crab_otus,crab_data$region,correction = "fdr")


# combined model by margin of island group (main,northwest) vs sample group (shallow main, deep main, shallow northwest)
combined_model_margin <- adonis2(crab_dist ~ island_group + sample_group, data=crab_data,by="margin")
combined_model_margin 

adonis2(crab_dist ~ sample_group + island_group, data=crab_data)
adonis2(crab_dist ~ island_group * sample_group, data=crab_data)

# the same thing by term, I guess
combined_model_term <- adonis2(crab_dist ~ island_group + sample_group, data=crab_data,by="term")
combined_model_term 

# what if we try this same thing as a dbRDA?
combined_model_db <- dbrda(crab_dist ~ island_group + sample_group, data=crab_data)

# show the significance by margin and term
anova(combined_model_db,by="margin",permutations = 999)
anova(combined_model_db,by="term",permutations = 999)

# plot the dbrda (it won't work becauuse we have the wrong sample data)
plot_ordination(crabs,combined_model_db,color="sample_group",shape="depth_zone") +
  stat_ellipse(aes(group=sample_group))

# summary data for all ARMS -----------------------------------------------
cc <- crabs_untransformed

cat("total individuals:\n")
cc %>%
  sample_sums() %>%
  sum()

cat("abundance range:\n")
cc %>%
  sample_sums() %>%
  range()

cat("abundance mean:\n")
cc %>%
  sample_sums() %>%
  mean()

cat("abundance SD:\n")
cc %>%
  sample_sums() %>%
  sd()

# get 5 most common species


cat("most abundant 5 species:\n")
top_5 <- cc %>%
  taxa_sums() %>%
  sort() %>%
  rev() %>%
  head(5) 
top_5

otus <- top_5 %>%
  names() %>%
  str_replace_all("_"," ") -> otus

cat("number of units where top 5 species occur:\n")
cc %>% 
  subset_taxa(species %in% otus) %>%
  ps_standardize("pa") %>%
  taxa_sums() 

cat("number of singletons:\n")
ts <- cc %>% taxa_sums()
sum(ts == 1)

shallow_singletons <- enframe(ts[ts == 1]) %>%
  arrange(name)


tt <- cc %>%
  taxa_tibble() %>%
  mutate(thing = "everything")

cat("taxonomy breakdown by id type:\n")
tt %>%
  group_by(taxon_level) %>%
  summarise(families = n_distinct(family), genera = n_distinct(genus), species = n_distinct(species))

cat("general taxonomy breakdown:\n")
tt %>%
  group_by(thing) %>%
  summarise(families = n_distinct(family), genera = n_distinct(genus), species = n_distinct(species))

# summary data for deep ARMS ----------------------------------------------

cc <- crabs_untransformed %>%
  subset_samples(depth_zone != "shallow")

cc <- prune_taxa(taxa_sums(cc) > 0,cc)

cat("total individuals:\n")
cc %>%
  sample_sums() %>%
  sum()

cat("abundance range:\n")
cc %>%
  sample_sums() %>%
  range()

cat("abundance mean:\n")
cc %>%
  sample_sums() %>%
  mean()

cat("abundance SD:\n")
cc %>%
  sample_sums() %>%
  sd()

# get 5 most common species


cat("most abundant 5 species:\n")
top_5 <- cc %>%
  taxa_sums() %>%
  sort() %>%
  rev() %>%
  head(5) 
top_5

otus <- top_5 %>%
  names() %>%
  str_replace_all("_"," ") -> otus

cat("number of units where top 5 species occur:\n")
cc %>% 
  subset_taxa(species %in% otus) %>%
  ps_standardize("pa") %>%
  taxa_sums() 

cat("number of singletons:\n")
ts <- cc %>% taxa_sums()
sum(ts == 1)
sum(ts == 1) / sum(ts > 0)

deep_singletons <- enframe(ts[ts == 1]) %>%
  arrange(name)


tt <- cc %>%
  taxa_tibble() %>%
  mutate(thing = "everything")

cat("taxonomy breakdown by id type:\n")
tt %>%
  group_by(taxon_level) %>%
  summarise(families = n_distinct(family), genera = n_distinct(genus), species = n_distinct(species))

cat("general taxonomy breakdown:\n")
tt %>%
  group_by(thing) %>%
  summarise(families = n_distinct(family), genera = n_distinct(genus), species = n_distinct(species))

}