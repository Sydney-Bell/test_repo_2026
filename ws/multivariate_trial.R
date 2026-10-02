## ---- libraries
library(tidyverse)
library(vegan)
library(GGally)
library(corrplot)
library(car)
library(ggvegan)
library(ggrepel)
library(scales)
## ----end


## ---- read dune
dune <- read_csv(file = "./data/dune.csv", trim_ws = TRUE)
glimpse(dune)
head(dune)
#do farming practices impact plant community?

#MANAGEMENT = cat
# plant species = count data

#go with a PCoA
## ----end

## ---- PCoA
## principal coordinates analysis or metric multidimensional scaling

#need to first remove categorical variable from dune data
man <- dune |> 
  dplyr::select(where(is.character) | where(is.factor)) |> 
  mutate(MANAGEMENT = factor(MANAGEMENT, levels = c("NM", "BF", "HF", "SF")))

plants <- dune |> 
  dplyr::select(where(is.numeric))


dune.dist <- vegdist(plants, method='bray')
dune.capscale <- capscale(dune.dist~1, data=man) #runs analysis
summary(dune.capscale, display=NULL) #go with first 4 axes
plot(dune.capscale)
autoplot(dune.capscale, geom='text')

# Distance based redundancy analysis
dune.capscale2 <- capscale(dune.dist ~
                              MANAGEMENT,
                            data = man)
summary(dune.capscale2, display=NULL) #not explaining well, needs 6 axes

## ---- Try MDS
dune <- dune |> 
  mutate(MANAGEMENT = factor(MANAGEMENT, levels = c("NM", "BF", "HF", "SF")))

dune.mds <- metaMDS(dune[,-1], k=2,  plot=TRUE) #k = number of axes, -1 removes the first column because its not a bird
dune.mds
#good stress level ~0.11

dune.mds$stress #want stress <0.1
stressplot(dune.mds)

plot(dune.mds)

dune.mds.scores <- dune.mds |>
  fortify() |>
  full_join(dune |>
              rownames_to_column(var='label'),
            by =  'label') #turns row names in to columns

g <-
  ggplot(data = NULL, aes(y=nmds2, x=nmds1)) +
  geom_hline(yintercept=0, linetype='dotted') +
  geom_vline(xintercept=0, linetype='dotted') +
  geom_point(data=dune.mds.scores %>% filter(score=='sites'),
             aes(color=MANAGEMENT)) +
  geom_text(data=dune.mds.scores %>% filter(score=='sites'),
            aes(label=label, color=MANAGEMENT), hjust=-0.2, show.legend = FALSE) +
  geom_segment(data=dune.mds.scores %>% filter(score=='species'),
               aes(y=0, x=0, yend=nmds2, xend=nmds1),
               arrow=arrow(length=unit(0.3,'lines')), color='red',
               alpha =  0.2) +
  geom_text(data=dune.mds.scores %>% filter(score=='species'),
            aes(y=nmds2*1.1, x=nmds1*1.1, label=label), color='red',
            alpha =  0.2)
g

g2 <- g +
  stat_density_2d(
    data = dune.mds.scores %>% filter(score == "sites"),
    geom = "polygon", #this makes kernel density
    aes(y = nmds2, x = nmds1, fill = MANAGEMENT),
    contour_var = "ndensity",
    breaks = c(0.05, 0.1),
    alpha = 0.3,
    position = "identity",
    show.legend = FALSE
  )
g2

centroids <- dune.mds.scores |>
  filter(score == "sites") |>
  group_by(MANAGEMENT) |>
  summarise(across(c(nmds1, nmds2), list(c = mean)))

dune.mds.scores <- dune.mds.scores |>
  full_join(centroids)

dune.mds.scores.centroids <- dune.mds.scores |>
  filter(score == "sites") |>
  group_by(MANAGEMENT) |>
  summarise(across(c(nmds1, nmds2), list(c = mean)))
dune.mds.scores <- dune.mds.scores |>
  full_join(dune.mds.scores.centroids)

#spider plot

g1 <-
     ggplot(data = NULL, aes(y=nmds2, x=nmds1)) +
     geom_hline(yintercept=0, linetype='dotted') +
     geom_vline(xintercept=0, linetype='dotted') +
     geom_point(data=dune.mds.scores %>% filter(score=='sites'),
                               aes(color=MANAGEMENT))
g1

g.plot <- g1 +
  stat_density_2d(
    data = dune.mds.scores %>% filter(score == "sites"),
    geom = "polygon",
    aes(y = nmds2, x = nmds1, fill = MANAGEMENT),
    contour_var = "ndensity",
    breaks = c(0.2),
    alpha = 0.1,
    position = "identity",
    show.legend = FALSE) +
  geom_segment(data = dune.mds.scores,
               aes(x = nmds1_c, xend = nmds1, y = nmds2_c, yend = nmds2, colour = MANAGEMENT)) +
  scale_x_continuous(limits = c(-2, 2)) +
  scale_y_continuous(limits = c(-2, 2)) +
  theme_classic()
g.plot

Xmat <- model.matrix(~1+MANAGEMENT, data = dune)
colnames(Xmat) <-gsub("MANAGEMENT","",colnames(Xmat))
envfit <- envfit(dune.mds, env=Xmat)
envfit #how different is centroid from NM? Only SF (standard farming) is different from NM

#pretty plot
dune.env.scores <- envfit |> fortify()
g3 <- g.plot +
  geom_segment(data=dune.env.scores,
               aes(y=0, x=0, yend=NMDS2, xend=NMDS1),
               arrow=arrow(length=unit(0.3,'lines')), color='blue') +
  geom_text(data=dune.env.scores,
            aes(y=NMDS2*1.1, x=NMDS1*1.1, label=label), color='blue')
g3

#beta dispersion to check if dispersions are equal
dune.disp <- betadisper(dune.dist, dune$MANAGEMENT)
boxplot(dune.disp)
plot(dune.disp)
anova(dune.disp) #evidence of equal dispersion
permutest(dune.disp, pairwise = TRUE)
TukeyHSD(dune.disp)
#Dispersion is equal across all groups

#Figure out what management groups are different from each other
mm <-  model.matrix(~-1 + MANAGEMENT, data=dune)
head(mm)
colnames(mm) <-gsub("MANAGEMENT","",colnames(mm))
mm <- data.frame(mm)

dune.adonis<-adonis2(dune.dist ~
                           NM + BF + HF + SF,
                         data=mm,
                         by = "terms",
                         perm=9999)
print(dune.adonis)
#NM has a significantly different community than the other management types

#pairwise tests
permutest(dune.disp, pairwise = TRUE)

#SIMPER test: similarity test
#for the different management comparisons, which species are the most responsive?
dune.std <- wisconsin(dune[,c(-1)])
simper(dune.std, dune$MANAGEMENT) |> summary()
#Lolper/Antodo was the most responsive in all of the NM comparisons

## ----end





