# set the random seed
set.seed(31337)

# collapse taxonomy and crabs tables by unknowns
source("collapse_taxonomy.R")
# functions ---------------------------------------------------------------


# plot an ordination with various options
plot_ord <- function(ps,ord,arrows=TRUE,sig=FALSE,by="term",perm=999,alpha=0.05,scale=1,labelsize=5,pointsize=3,term_map=NULL,...) {

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

    arrowdf <- arrowdf %>%
      mutate(
        hjust = if_else(.data[[str_c(axis,"1")]] < 0,1,0),
        xnudge = if_else(.data[[str_c(axis,"1")]] < 0,-0.02,0.02)
      )

    if (!is.null(term_map)) {
      arrowdf <- arrowdf %>%
        inner_join(term_map,by=c("labels" = "variable")) %>%
        select(-labels) %>%
        rename(labels=display)
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
    label_map <- aes(x = .data[[str_c(axis,"1")]]+xnudge,
                     y = .data[[str_c(axis,"2")]],
                     shape = NULL,
                     color = NULL,
                     label = labels,
                     hjust = hjust)
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
      # geom_text_repel(
      # geom_text(
      # geom_label(
      geom_richtext(
        mapping = label_map,
        size = labelsize,
        data = arrowdf,
        show.legend = FALSE,
        color='black',
        fill="grey96"
      )
  }
  return(plotz)
}

# do a pairwise PERMANOVA analysis by some factor
pairwise_adonis <- function(comm, factors, permutations = 1000, correction = "fdr", method = "bray") {
  # get possible pairwise factor combinations
  factor_combos <- combn(unique(as.character(factors)), 2)
  # map through factor combinations and run model for each pair,
  # dumping the output into a tibble
  model_output <- map_dfr(array_branch(factor_combos,2), ~{
    fact <- factors[factors %in% .x]

    # get our factor specific community or distance matrix
    if (inherits(comm,'dist')) {
      dd <- as.dist(as.matrix(comm)[factors %in% .x, factors %in% .x])
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
    here("data","generated","crabs_collapsed.csv"),"unit",
    here("data","generated","taxonomy_collapsed.csv"),"otu",
    # here("data","crabs.csv"),"unit",
    # here("data","taxonomy.csv"),"otu",
    here("data","metadata.csv")
  )

  island_transform <- tribble(
    ~island,~hawaiian_name,
    "Kure Atoll","Hōlanikū",
    "Pearl and Hermes Reef","Manawai",
    "Lisianski","Kapou",
    "French Frigate Shoals","Lalo",
    "Kauai","Kaua‘i",
    "Oahu","O‘ahu",
    "Maui","Maui",
    "Hawaii Island","Hawai‘i Island"
  )

  cd <- sample_tibble(cc$crabs_untransformed) %>%
    left_join(island_transform,by="island") %>%
    mutate(
      island = fct_reorder(island,-lat),
      hawaiian_name = fct_reorder(hawaiian_name,-lat)
    ) %>%
    column_to_rownames("sample")

  sample_data(cc$crabs_untransformed) <- cd

  cc$crabs <- cc$crabs_untransformed %>%
    ps_standardize("hellinger")

  # shallow subset
  cc$crabs_shallow <- prune_samples(
    sample_tibble(cc$crabs)$shallow_deep == "shallow",
    cc$crabs
  )
  # prune taxa
  cc$crabs_shallow <- prune_taxa(taxa_sums(cc$crabs_shallow) > 0,cc$crabs_shallow)

  # deep subset
  cc$crabs_deep <- prune_samples(
    sample_tibble(cc$crabs)$shallow_deep == "deep",
    cc$crabs
  )
  # prune taxa
  cc$crabs_deep <- prune_taxa(taxa_sums(cc$crabs_deep) > 0,cc$crabs_deep)

  # untransformed deep subset
  cc$crabs_deep_untransformed <- prune_samples(
    sample_tibble(cc$crabs_untransformed)$shallow_deep == "deep",
    cc$crabs_untransformed
  )
  # prune deep crabs taxa
  cc$crabs_deep_untransformed <- prune_taxa(taxa_sums(cc$crabs_deep_untransformed) > 0,cc$crabs_deep_untransformed)

  # untransformed shallow subset
  cc$crabs_shallow_untransformed <- prune_samples(
    sample_tibble(cc$crabs_untransformed)$shallow_deep == "shallow",
    cc$crabs_untransformed
  )
  # prune shallow crabs taxa
  cc$crabs_shallow_untransformed <- prune_taxa(taxa_sums(cc$crabs_shallow_untransformed) > 0,cc$crabs_shallow_untransformed)


  # save an unscaled version of the shallow crab metadata
  cc$crab_data_shallow_unscaled <- cc$crabs_shallow %>%
    sample_tibble(sample_col = "sample")

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
    select(unit,region,island_group,deployment_year,recovery_year,soak_time,island,hawaiian_name,lat,lon,depth,contains("chl"),contains("sst"),slope,coral_cover,closest_island,larval_immigration,contains('impact')) %>%
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


  cc$crabs_oahu <- prune_samples(
    sample_tibble(cc$crabs)$island == "Oahu",
    cc$crabs
  )
  # prune taxa
  cc$crabs_oahu <- prune_taxa(taxa_sums(cc$crabs_oahu) > 0,cc$crabs_oahu)
  cc$crab_dist_oahu <- distance(cc$crabs_oahu,"bray")

  return(cc)
}

# make a species be cf.
cf_species <- function(spp,deitalicize=TRUE) {
  if (deitalicize) {
    str_replace(spp," ","* cf. *")
  } else {
    str_replace(spp," "," cf. ")
  }
}

sample_summary <- function(ps,top_n=5,taxonomy=NULL) {
  ss <- list()

  # make otu lookup table
  otu_lookup <- ps %>%
    taxa_tibble(otu_col = "otu") %>%
    mutate(species = as.character(case_when(
      as.logical(complex) ~ cf_species(species),
      .default = species
    ))) %>%
    select(otu,species) %>%
    deframe()

  # basic numbers

  # total number of individuals
  ss$total <- ps %>%
    sample_sums() %>%
    sum()

  # overall range of abundance
  ss$abundance_range <- ps %>%
    sample_sums() %>%
    range()

  # overall abundance mean
  ss$abundance_mean <- ps %>%
    sample_sums() %>%
    mean()

  # overall abundance sd
  ss$abundance_sd <- ps %>%
    sample_sums() %>%
    sd()

  # get 5 most common species
  top <- ps %>%
    taxa_sums() %>%
    sort() %>%
    rev() %>%
    head(top_n)

  ss$top_otus <- names(top)

  top <- top %>%
    set_names(otu_lookup[names(.)])

  ss$top_taxa <- names(top)
  ss$top_n <- top

  ss$top_units <- prune_taxa(ss$top_otus,ps) %>%
    ps_standardize("pa") %>%
    taxa_sums() %>%
    sort() %>%
    rev() %>%
    set_names(otu_lookup[names(.)])

  # make it the same order as the top 5 species
  ss$top_units <- ss$top_units[match(names(ss$top_units),names(top))]

  ss$top <- taxonomy %>%
    filter(otu %in% ss$top_otus) %>%
    slice(match(ss$top_otus,otu)) %>%
    mutate(
      disp = case_when(
        as.logical(complex) & taxon_level == "species" ~ str_glue("{cf_species(display_species)} {authority}"),
        taxon_level == "species" ~ str_glue("{display_species} {authority}"),
        .default = display_species
      )
    ) %>%
    pull(disp)

  ts <- ps %>% taxa_sums()

  ss$singleton_count <- sum(ts == 1)

  ss$singleton_species <-  ts[ts == 1] %>%
    names() %>%
    sort() %>%
    map_chr(~otu_lookup[.x])

  tt <- ps %>%
    taxa_tibble()

  ss$types_breakdown <- tt %>%
    count(taxon_level) %>%
    deframe() %>%
    as.list()

  ss$taxonomy_breakdown <- list(
    families = n_distinct(tt$family),
    genera = n_distinct(tt$genus),
    species = n_distinct(grep(" sp\\.$",tt$species,invert = TRUE,value=TRUE)),
    all_species = n_distinct(tt$species)
  )

  return(ss)
}

alpha_diversity <- function(ps,measures=c("Observed","InvSimpson")) {

  ad <- list()

  # skip fisher because it's broken
  richness <- ps %>%
    estimate_richness(measures = measures) %>%
    as_tibble(rownames="sample") %>%
    rename(Simpson = InvSimpson) %>%
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
  try(ad$t_r_depth <- t.test(Observed ~ shallow_deep, data=all_richness),silent=TRUE)
  try(ad$t_s_depth <- t.test(Simpson ~ shallow_deep, data=all_richness),silent=TRUE)

  try(ad$k_r_island <- kruskal.test(Observed ~ island_group, data=all_richness),silent=TRUE)
  try(ad$k_s_island <- kruskal.test(Simpson ~ island_group, data=all_richness),silent=TRUE)
  try(ad$k_r_depth <- kruskal.test(Observed ~ shallow_deep, data=all_richness),silent=TRUE)
  try(ad$k_s_depth <- kruskal.test(Simpson ~ shallow_deep, data=all_richness),silent=TRUE)
  try(ad$k_r_depth_main <- kruskal.test(Observed ~ shallow_deep, data=all_richness %>% filter(island_group == "main")),silent=TRUE)
  try(ad$k_s_depth_main <- kruskal.test(Simpson ~ shallow_deep, data=all_richness  %>% filter(island_group == "main")),silent=TRUE)
  try(ad$k_r_depth_oahu <- kruskal.test(Observed ~ shallow_deep, data=all_richness %>% filter(island == "Oahu")),silent=TRUE)
  try(ad$k_s_depth_oahu <- kruskal.test(Simpson ~ shallow_deep, data=all_richness  %>% filter(island == "Oahu")),silent=TRUE)

  try(ad$w_r_island <- wilcox.test(Observed ~ island_group, data=all_richness),silent=TRUE)
  try(ad$w_s_island <- wilcox.test(Simpson ~ island_group, data=all_richness),silent=TRUE)
  try(ad$w_r_depth <- wilcox.test(Observed ~ shallow_deep, data=all_richness),silent=TRUE)
  try(ad$w_s_depth <- wilcox.test(Simpson ~ shallow_deep, data=all_richness),silent=TRUE)

  return(ad)
}
