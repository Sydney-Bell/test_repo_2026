############################################################
##             RDA (Correlations, Euclidean distance)
##            /  \
##Unconstrained   Constrained ---> ordination
##   (PCA)            (RDA)   ---> anova
##
##              CA (Chisq distance)
##             /  \
##Unconstrained   Constrained ---> ordination
## (CA)             (CCA)     ---> anova
##
##             PCoA (any distance)
##             /  \
##Unconstrained   Constrained ---> ordination
##                            ---> anova
##
##             dbRDA (any distance)
##             /  \
##Unconstrained   Constrained ---> ordination
##                            ---> anova
##
##Unconstrained  ---> ordination
##               ---> envfit (overlay enviromental data) (permutation test) important that dispersion is equal
##               ---> lm/glm etc (response or predictor)
#############################################################
##     Dissimilarity
##            --> MDS      ---> ordination
##            --> bioenv   ---> variable importance       (perm test) helps select important variables, but not the best method
##            --> adonis2* ---> anova                     (perm test) takes dist matrix and tries to look for patterns in that, important that dispersion is equal
##            --> simper   ---> similarity percentages
##            --- betadisp ---> homogeneity of dispersion (perm test)
#############################################################
##     Model based ordination
##            ---> glmmTMB (via reduced rank / latent variable)
##            ---> gllvm (generalised latent variable models)
#############################################################
##     Model based
##            ---> manyglm ---> anova
##            ---> gllvm (generalized latent variable models)
#############################################################

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


## ---- read Spider
spider.abund <- read_csv(file = "./data/spider.abund.csv", trim_ws = TRUE)
spider.env <- read_csv(file = "./data/spider.env.csv", trim_ws = TRUE)
glimpse(spider.abund)
glimpse(spider.env)
## ----end

## ---- EDA spider
#want to see how correlated species are, if none of them are...then don't do multivariate stats
spider.abund |>
  cor() |>
  corrplot(type = 'upper',
    diag = FALSE)
## And now with axes arrange according to first principal component axis
spider.abund |>
  cor() |>
  corrplot(type = 'upper',
    order = 'FPC',
    diag = FALSE)
## ----end

## ---- EDA 2
#assumptions for PCA are normal, linearity, homogeneity of variance
spider.abund |>
  ggpairs(lower = list(continuous = "smooth"),
    diag = list(continuous = "density"),
    axisLabels = "show")
#not normal because the distribution graphs in the middle aren't a bell curve
#this causes problems for our other assumptions too

#can use square root or fourth root transformations here because we have 0s and dont need to back transform
spider.abund^0.25 |>
  ggpairs(lower = list(continuous = "smooth"),
          diag = list(continuous = "density"),
          axisLabels = "show")
#not perfect but looks better
#still not meeting linearity assumption
#this data wouldn't be good for PCA but just using it as an example/test
## ----end

## PCA ------------------------------------------------------
#can use standardization to even out the contribution of species and sites
#a typical way to standardize data is to take a column and divide it by its maximum
#then divide rows by sum
#this processes is so common its called the wisconsin double standardization
spider.std <- spider.abund |>
  mutate(across(everything(), ~.x^0.25)) |>
  wisconsin()
#OR
spider.std <- (spider.abund^0.25) |>
  wisconsin()
spider.std

spider.rda <- rda(spider.abund, scale=TRUE) #this is the pca, scale = TRUE is saying yes to correlation
#correlation is -1 to 1, it is covariance that has been standardized to a set range
#covariance is infinite 
#but murray recommends doing your own standardization and then saying scale=FALSE
summary(spider.rda, display=NULL)
#inertia = units of variation
#first axis is explaining 5 units of variation across all 12 species and 43% of variation
#need to decide how many axes to keep --> anything >1 unit of variance (more than 1 eigenvalue rule) IF scale = TRUE
#if scale=FALSE, then need to add up Eigenvalues and divide by # of axes
#OR just keep cumulative proportion of 80%

screeplot(spider.rda)
abline(a=1,b=0)

#now look at standardized data with scale=FALSE
spider.rda <- rda(spider.std, scale=FALSE)
summary(spider.rda, display=NULL)

## skip
scores(spider.rda, choices=1:3, display='sites')
scores(spider.rda, choices=1:3, display='species')

## Quick and nasty ordination plots
biplot(spider.rda, scaling='species')
biplot(spider.rda, scaling='sites')
#Alopacce has high correlation with PC1, but very low correlation with PC2

## Quick and nasty ordination plots
pl <- vegan::ordiplot(spider.rda)
points(pl, "sites", pch=21, col="red", bg="yellow")
text(pl, "sites", col="red", cex=0.9)
text(pl, "species", col="blue", cex=0.9)

autoplot(spider.rda)
autoplot(spider.rda) + theme_bw()
autoplot(spider.rda,geom='text') + theme_bw()


spider.rda.scores <- spider.rda |>
  fortify() #fortify makes the data able to be put into ggplot, extracts scores and component loadings
spider.rda.scores #values are the arrowheads and corresponds to the biplot above

g <-
  ggplot(data = NULL, aes(y=pc2, x=pc1)) +
  geom_hline(yintercept=0, linetype='dotted') +
  geom_vline(xintercept=0, linetype='dotted') +
  geom_point(data=spider.rda.scores |> filter(score=='sites')) +
  geom_text(data=spider.rda.scores |> filter(score=='sites'),
    aes(label=label), hjust=-0.2) +
  geom_segment(data=spider.rda.scores |> filter(score=='species'),
    aes(y=0, x=0, yend=pc2, xend=pc1),
    arrow=arrow(length=unit(0.3,'lines')), color='red') +
  geom_text_repel(data=spider.rda.scores |> filter(score=='species'),
    aes(y=pc2*1.1, x=pc1*1.1, label=label), color='red') +
  theme_bw()
g

## Nice axes titles
eig <- eigenvals(spider.rda)

g <- g +
  scale_y_continuous(paste(names(eig[2]),
    sprintf('(%0.1f%% explained var.)', #sprintf is used to combine numbers and text
    100 * eig[2]/sum(eig))))+
  scale_x_continuous(paste(names(eig[1]),
    sprintf('(%0.1f%% explained var.)',
    100 * eig[1]/sum(eig)))) #doing this to convert eigs to %
g

#put a circle
circle.prob <- 0.95 #would normally be 0.68 or 0.89 (represent stdev)
r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(spider.rda$CA$u[,1:2]^2))^(1/4)
theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
circle <- data.frame(PC1 = r * cos(theta), PC2 = r * sin(theta))
g <- g +
  geom_path(
    data = circle,
    aes(y = PC2, x = PC1),
    color = muted("white"), size = 1 / 2, alpha = 1 / 3
  )
g

## ----end


## ---- Envfit
spider.env |>
  cor() |>
  corrplot(type = 'upper',
    order = 'FPC',
    diag = FALSE)

spider.env |>
  ggpairs(lower = list(continuous = "smooth"),
    diag = list(continuous = "density"),
    axislabels = "show")
spider.envfit <- envfit(spider.rda, env = spider.env)
spider.envfit #community composition does change with habitat
#permutation test = shuffles predictor or response to get a bunch of values, makes a distribution, and you count how many values are greater than the original value
#does not assume variance, independence, normality


spider.env.scores <- spider.envfit |>
  fortify() |>
  mutate(Flag = factor(ifelse(sqrt(PC1^2 + PC2^2) > r, 1, 0)))
g <- g +
  geom_segment(data=spider.env.scores,
    aes(y=0, x=0, yend=PC2, xend=PC1, alpha = Flag, show.legend = FALSE),
    arrow=arrow(length=unit(0.3,'lines')), color='blue') +
  geom_text(data=spider.env.scores,
    aes(y=PC2*1.1, x=PC1*1.1, label=label, alpa = Flag),
    color='blue', show.legend = FALSE)
g

## ----end

## ---- lm
pc1 <- spider.rda.scores |> filter(score=='sites') |> pull(pc1)
pc2 <- spider.rda.scores |> filter(score=='sites') |> pull(pc2)

#could do something like this: glm(pc1 ~soil.dry, data = spider.env, family = "gaussian")
#how does community composition (pc1) relate to soil.dry

lm(1:nrow(spider.env) ~ soil.dry + bare.sand + fallen.leaves +
     moss + herb.layer + reflection, data =  spider.env) |>
  vif()
lm(1:nrow(spider.env) ~ herb.layer + fallen.leaves + bare.sand + moss, data=spider.env) |>
  vif()
lm(pc1 ~ herb.layer + fallen.leaves + bare.sand + moss, data=spider.env) |>
  summary()
lm(pc2 ~ herb.layer + fallen.leaves + bare.sand + moss, data=spider.env) |>
  summary()

## ----end

## ---- RDA spiders

## ---- RDA
spider.rda <- rda(
  spider.std ~
    scale(herb.layer) +
    scale(fallen.leaves) +
    scale(bare.sand) +
    scale(moss),
  data = spider.env,
  scale = FALSE
)
vif.cca(spider.rda)
#constrained = how much do our proposed predictors explain? 69%
summary(spider.rda, display=NULL)
# 4 constrained axes (for the 4 predictors we added in), the rest are unconstrained

#biplot when you have predictor variables included is looking at how things relate in relation to predictor variables
#RDA is basically optimized for certain predictors you suspect or know are important

## ----end

## ---- goodness of fit
#if one species is ALWAYS present, you could remove it to make sure its not diluting effects
#if any were lower than 0.1 across all axes, you would remove it, maybe unresponsive to all predictors
goodness(spider.rda)
goodness(spider.rda, display = "sites")
inertcomp(spider.rda)
inertcomp(spider.rda, proportional = TRUE)
## ----end

## ---- Anova
anova(spider.rda)
anova(spider.rda, by='axis') #which axes relate to our predictors?
anova(spider.rda, by='margin') #which env variables have impacts on communities?
## ----end

## ---- Variance inflation factors
#how correlated this predictor is to all of the others
#values >5 are considered bad --> strongly correlated
vif.cca(spider.rda)
## ----end

## ---- other parameters
coef(spider.rda)
RsquareAdj(spider.rda)
## ----end

## ---- ordination plot
screeplot(spider.rda)
autoplot(spider.rda, geom='text')
## ----end

## ---- CA
#standardize columns for relative abundance
spider.std <- (spider.abund^0.25) |>
  wisconsin()
spider.std
spider.ca <- cca(spider.std, scale=FALSE) #unconstrained because no predictors
## ----end

## ---- CA summary
summary(spider.ca, display=NULL)
#keep CA1 and CA2
## ----end

## ---- CA ordination plot
screeplot(spider.ca)
sum(eigenvals(spider.ca))/length(eigenvals(spider.ca))
eigenvals(spider.ca)/sum(eigenvals(spider.ca))
plot(spider.ca, scaling='species')

autoplot(spider.ca)
autoplot(spider.ca) + theme_bw()
autoplot(spider.ca, geom='text') + theme_bw()
## ----end

## ---- CA ordination plot pretty
spider.ca.scores <- spider.ca |>
  fortify()
spider.ca.scores |> head()

g <-
  ggplot(data = NULL, aes(y=ca2, x=ca1)) +
  geom_hline(yintercept=0, linetype='dotted') +
  geom_vline(xintercept=0, linetype='dotted') +
  geom_point(data=spider.ca.scores %>% filter(score=='sites')) +
  geom_text(data=spider.ca.scores %>% filter(score=='sites'),
    aes(label=label), hjust=-0.2) +
  geom_segment(data=spider.ca.scores %>% filter(score=='species'),
    aes(y=0, x=0, yend=ca2, xend=ca1),
    arrow=arrow(length=unit(0.3,'lines')), color='red') +
  ## geom_text(data=spider.rda.scores %>% filter(score=='species'),
  ##           aes(y=PC2*1.1, x=PC1*1.1, label=label), color='red') +
  geom_text_repel(data=spider.ca.scores %>% filter(score=='species'),
    aes(y=ca2*1.1, x=ca1*1.1, label=label), color='red') +
  theme_bw()
g

## ----end

## ---- CA envfit
spider.envfit <- envfit(spider.ca, env=spider.env)
spider.envfit
autoplot(spider.envfit)

spider.env.scores <- spider.envfit |> fortify()
g <- g +
  geom_segment(data=spider.env.scores,
    aes(y=0, x=0, yend=CA2, xend=CA1),
    arrow=arrow(length=unit(0.3,'lines')), color='blue') +
  geom_text(data=spider.env.scores,
    aes(y=CA2*1.1, x=CA1*1.1, label=label), color='blue')
g

## ---- PCoA
## principal coordinates analysis or metric multidimensional scaling
# CANNOT use these scores in other analyses, they are not independent
#data reduction can only be done for PCA and CA
spider.dist <- vegdist(spider.std, method='bray')
spider.capscale <- capscale(spider.dist~1, data=spider.env) #runs analysis
summary(spider.capscale, display=NULL)
plot(spider.capscale)
autoplot(spider.capscale, geom='text')

# Distance based redundancy analysis
spider.capscale <- capscale(spider.dist ~
    scale(herb.layer) +
    scale(fallen.leaves) +
    scale(bare.sand) +
    scale(moss),
  data = spider.env)
summary(spider.capscale, display=NULL)
plot(spider.capscale)

summary(spider.capscale, display=NULL)
anova(spider.capscale)

anova(spider.capscale, by='margin')
screeplot(spider.capscale)
sum(eigenvals(spider.capscale))/length(eigenvals(spider.capscale))
eigenvals(spider.capscale)/sum(eigenvals(spider.capscale))


## ---- MDS macnally
macnally <- read.csv('./data/macnally_full.csv',strip.white=TRUE)
head(macnally)
macnally <- macnally |>
  mutate(HABITAT = factor(HABITAT, levels = c(
    "Mixed", "Gipps.Manna",
    "Montane Forest", "Foothills Woodland", "Box-Ironbark", "River Red Gum"
  )))

macnally.mds <- metaMDS(macnally[,-1], k=2,  plot=TRUE) #k = number of axes, -1 removes the first column because its not a bird
macnally.mds
#purely about make an ordination that is optimized for the number of dimension (ideally 2)

macnally.mds$stress #want stress <0.1
stressplot(macnally.mds)

plot(macnally.mds)

macnally.mds.scores <- macnally.mds |>
  fortify() |>
  full_join(macnally |>
             rownames_to_column(var='label'),
    by =  'label') #turns row names in to columns

g <-
    ggplot(data = NULL, aes(y=nmds2, x=nmds1)) +
    geom_hline(yintercept=0, linetype='dotted') +
    geom_vline(xintercept=0, linetype='dotted') +
    geom_point(data=macnally.mds.scores %>% filter(score=='sites'),
               aes(color=HABITAT)) +
    geom_text(data=macnally.mds.scores %>% filter(score=='sites'),
              aes(label=label, color=HABITAT), hjust=-0.2, show.legend = FALSE) +
    geom_segment(data=macnally.mds.scores %>% filter(score=='species'),
                 aes(y=0, x=0, yend=nmds2, xend=nmds1),
                 arrow=arrow(length=unit(0.3,'lines')), color='red',
      alpha =  0.2) +
    geom_text(data=macnally.mds.scores %>% filter(score=='species'),
      aes(y=nmds2*1.1, x=nmds1*1.1, label=label), color='red',
      alpha =  0.2)
g


g1 <-
    ggplot(data = NULL, aes(y=nmds2, x=nmds1)) +
    geom_hline(yintercept=0, linetype='dotted') +
    geom_vline(xintercept=0, linetype='dotted') +
    geom_point(data=macnally.mds.scores %>% filter(score=='sites'),
               aes(color=HABITAT))
g1

g2 <- g +
  stat_density_2d(
    data = macnally.mds.scores %>% filter(score == "sites"),
    geom = "polygon", #this makes kernel density
    aes(y = nmds2, x = nmds1, fill = HABITAT),
    contour_var = "ndensity",
    breaks = c(0.05, 0.1),
    alpha = 0.3,
    position = "identity",
    show.legend = FALSE
  )
g2

centroids <- macnally.mds.scores |>
  filter(score == "sites") |>
  group_by(HABITAT) |>
  summarise(across(c(nmds1, nmds2), list(c = mean)))

macnally.mds.scores <- macnally.mds.scores |>
  full_join(centroids)

macnally.mds.scores.centroids <- macnally.mds.scores |>
  filter(score == "sites") |>
  group_by(HABITAT) |>
  summarise(across(c(nmds1, nmds2), list(c = mean)))
macnally.mds.scores <- macnally.mds.scores |>
  full_join(macnally.mds.scores.centroids)

#spider plot
g1 <- g1 +
  stat_density_2d(
    data = macnally.mds.scores %>% filter(score == "sites"),
    geom = "polygon",
    aes(y = nmds2, x = nmds1, fill = HABITAT),
    contour_var = "ndensity",
    breaks = c(0.05, 0.1),
    alpha = 0.1,
    position = "identity",
    show.legend = FALSE) +
  geom_segment(data = macnally.mds.scores,
  aes(x = nmds1_c, xend = nmds1, y = nmds2_c, yend = nmds2, colour = HABITAT)) +
  theme_classic()
g1

#pretty plot
macnally.env.scores <- envfit |> fortify()
g3 <- g1 +
  geom_segment(data=macnally.env.scores,
               aes(y=0, x=0, yend=NMDS2, xend=NMDS1),
               arrow=arrow(length=unit(0.3,'lines')), color='blue') +
  geom_text(data=macnally.env.scores,
            aes(y=NMDS2*1.1, x=NMDS1*1.1, label=label), color='blue')
g3

Xmat <- model.matrix(~-1+HABITAT, data = macnally) #-1 = no intercept, columns are just the means of each categorical level
colnames(Xmat) <-gsub("HABITAT","",colnames(Xmat))
envfit <- envfit(macnally.mds, env=Xmat)
envfit #how different is centroid from the general community?
#this analysis assumes within a habitat, points are equally dispersed
#BUT it cant tell whether habitats are different because of true differences or because of disperison
#envfit overlays env data onto your data (response and predictors)

Xmat2 <- model.matrix(~1+HABITAT, data = macnally)
colnames(Xmat2) <-gsub("HABITAT","",colnames(Xmat2))
envfit2 <- envfit(macnally.mds, env=Xmat2)
envfit2 #how different is centroid from Mixed community?

#beta dispersion to check if dispersions are equal
macnally.disp <- betadisper(macnally.dist, macnally$HABITAT)
boxplot(macnally.disp)
plot(macnally.disp)
anova(macnally.disp) #evidence of unequal dispersion
permutest(macnally.disp, pairwise = TRUE)
TukeyHSD(macnally.disp)
#you can say that mixed community is potentially more diverse

#PERMANOVA
#make dist matrix
macnally.dist <- vegdist(macnally[,-1], 'bray')

adonis2(macnally.dist ~ HABITAT, data=macnally)
#community composition does differ between habitats BUT could be from dispersion

#Figure out what habitats are different from each other
mm <-  model.matrix(~-1 + HABITAT, data=macnally)
head(mm)
colnames(mm) <-gsub("HABITAT","",colnames(mm))
mm <- data.frame(mm)

macnally.adonis<-adonis2(macnally.dist ~
                           Mixed + Box.Ironbark + Foothills.Woodland + Gipps.Manna +
                           Montane.Forest + River.Red.Gum,
  data=mm,
  by = "terms",
  perm=9999)
print(macnally.adonis)

#now try to remove Mixed
macnally.adonis<-adonis2(macnally.dist ~ Box.Ironbark + Foothills.Woodland + Gipps.Manna +
                           Montane.Forest + River.Red.Gum,
  data=mm,
  by = "margin",
  perm=9999)
print(macnally.adonis)
#doesn't have Mixed anymore. We likely have equal dispersion and its a fairer comparison

#pairwise tests
permutest(macnally.disp, pairwise = TRUE)

#SIMPER test: similarity test
#for the different habitat comparisons, which species are the most responsive?
macnally.std <- wisconsin(macnally[,c(-1)]^0.25)
simper(macnally.std, macnally$HABITAT) |> summary()
#when comparing Mixed and Gipps.Manna, the White throated honeyeater (WPHE) was the biggest difference
#only use this test to compare habitats that truly ARE DIFFERENT

## ----end
