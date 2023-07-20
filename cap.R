library(BiodiversityR)


crabdata <- sample_data(crabs) %>%
  as("data.frame")
crabmatrix <- otu_table(crabs) %>%
  as("matrix") %>%
  as.data.frame()

crabdata <- crabdata %>%
  dplyr::select(dplyr::where(~!any(is.na(.x)))) %>%
  mutate(across(where(is.character),as.factor))


cap <- CAPdiscrim(
  crabmatrix ~ island_group + sample_group,
  data=crabdata,
  dist="bray",
  axes=2,
  m=10,
  mmax=10,
  add=FALSE#,
  # permutations = 99
)

cd <- crabdata %>%
  as_tibble(rownames="sample")

capscores <- scores(cap) %>%
  as_tibble(rownames="sample")

capdata <- cap$PCoA %>%
  as_tibble(rownames="sample",.name_repair = ~c("pcoa1","pcoa2")) %>%
  left_join(cd,by="sample") %>%
  left_join(capscores,by="sample")

# quartz()
ggplot(capdata,aes(x=pcoa1,y=pcoa2)) + 
  geom_point(aes(shape=island_group,color=sample_group),size=3) +
  stat_ellipse(aes(group=sample_group,color=sample_group),level=0.95)

ggplot(capdata,aes(x=LD1,y=LD2)) + 
  geom_point(aes(shape=island_group,color=sample_group)) +
  stat_ellipse(aes(group=sample_group,color=sample_group),level=0.95) +
  scale_y_reverse()


plot1 <- ordiplot(cap, type="none")
ordisymbol(plot1, crabdata, "sample_group", legend=TRUE)

# plot change in classification success against m
plot(seq(1:14), rep(-1000, 14), xlim=c(1, 14), ylim=c(0, 100), xlab="m", 
     ylab="classification success (percent)", type="n")
for (mseq in 1:14) {
  CAPdiscrim.result <- CAPdiscrim(dune~Management, data=dune.env, 
                                  dist="bray", axes=2, m=mseq)
  points(mseq, CAPdiscrim.result$percent)
}


crab_data <- crab_data %>%
  mutate(sample_group = case_when(
    startsWith(sample_group,"deep") ~ depth_zone,
    TRUE ~ sample_group
  ))
