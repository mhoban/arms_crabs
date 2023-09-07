library(BiodiversityR)

otus <- otu_table(cc$crabs) %>%
  as.data.frame()
cd <- cc$crab_data %>%
  mutate(sample_group = factor(sample_group,ordered=FALSE)) %>%
  as.data.frame()



# cap <- CAPdiscrim(otus ~ sample_group, data=cd, dist="bray", axes=2, m=0, add=TRUE,permutations = 10)
cap <- CAPdiscrim(otus ~ island_group + depth_zone, data=cd, dist="bray", axes=2, m=0, add=TRUE)
cap <- CAPdiscrim(otus ~ sample_group, data=cd, dist="bray", axes=2, m=0, add=TRUE)

var <- cap$manova$Eigenvalues / sum(cap$manova$Eigenvalues)
var <- list(x = var[1], y = var[2])

nice <- function(s) {
  
}

pcoa <- cap$x %>%
  as_tibble(rownames="sample") %>%
  left_join(cd,by="sample") %>%
  mutate(
    sample_group = str_replace(sample_group,"shallow","Shallow"),
    sample_group = str_replace(sample_group,"main","(MHI)"),
    sample_group = str_replace(sample_group,"nwhi","(NWHI)"),
    sample_group = str_split(sample_group,"((?<=[0-9])(?=[a-zA-Z]))|(_)") %>%
      map_chr(~str_c(.x,collapse=" ")),
    sample_group = fct_reorder(sample_group,depth)
  )

pal <- c("Shallow (MHI)" = "#ffb400", "Shallow (NWHI)" = "#9080ff","30 m" = "#d87939", "60 m" = "#04c3c8", "90 m" = "#04738d")

# ggplot(pcoa,aes(x=LD1,y=LD2,color=sample_group)) + 
#   stat_ellipse(geom="polygon",aes(fill=sample_group,color=NULL),level=0.95,alpha=0.1) +
#   geom_point()  +
#   scale_color_manual(values=pal) + 
#   scale_fill_manual(values=pal)

ggplot(pcoa,aes(x=LD1,y=LD2,fill=sample_group)) + 
  stat_ellipse(geom="polygon",level=0.95,alpha=0.1,show.legend = FALSE) +
  geom_point(color="black",shape=21,size=5)  +
  # scale_color_manual(values=pal) + 
  scale_fill_manual(values=pal,name="Sample group") +
  xlab(str_glue("Linear discriminant 1 [{scales::percent(var$x,0.1)}]")) +
  ylab(str_glue("Linear discriminant 2 [{scales::percent(var$y,0.1)}]")) +
  theme_bw() +
  theme(
    legend.key = element_blank(),
    panel.grid = element_blank()
  )



