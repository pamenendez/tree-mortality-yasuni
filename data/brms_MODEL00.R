# Remove prior commands in R
rm(list=ls()) 

require(tidyr)
require(dplyr)
require(ggplot2)
require(bayesplot)
require(ape)
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
require(brms)
require(ggmcmc)
require(knitr)
library(gdata)
library(shiny)#PUEDE SER ESTO POR LO QUE NO SALE COMO CON JONATAHN?
library(shinystan)#PUEDE SER ESTO POR LO QUE NO SALE COMO CON JONATAHN?
library(parallel)

#DETECT NUMBER OD CORE
parallel::detectCores()


load("yasuni.census1.Rdata")
load("yasuni.census2A.Rdata")
yasuni.census2A00 <- merge(yasuni.census2A, yasuni.census1[, c("tag", "sp", "dbh")], by = c("tag", "sp"), all.x = TRUE)
# Completa los valores faltantes de dbh en yasuni.census2A con los valores de dbh de yasuni.census1
yasuni.census2A00$dbh.x[is.na(yasuni.census2A00$dbh.x)] <- yasuni.census2A00$dbh.y[is.na(yasuni.census2A00$dbh.x)]
# Elimina la columna dbh.y y renombra la columna dbh.x si es necesario
yasuni.census2A00 <- yasuni.census2A00[, !(names(yasuni.census2A00) %in% c("dbh.y"))]
colnames(yasuni.census2A00)[colnames(yasuni.census2A00) == "dbh.x"] <- "dbh"


load("yasuni.census2B.Rdata")
load("yasuni.census3A.Rdata")
yasuni.census3A00 <- merge(yasuni.census3A, yasuni.census2B[, c("tag", "sp", "dbh")], by = c("tag", "sp"), all.x = TRUE)
# Completa los valores faltantes de dbh en yasuni.census2A con los valores de dbh de yasuni.census1
yasuni.census3A00$dbh.x[is.na(yasuni.census3A00$dbh.x)] <- yasuni.census3A00$dbh.y[is.na(yasuni.census3A00$dbh.x)]
# Elimina la columna dbh.y y renombra la columna dbh.x si es necesario
yasuni.census3A00 <- yasuni.census3A00[, !(names(yasuni.census3A00) %in% c("dbh.y"))]
colnames(yasuni.census3A00)[colnames(yasuni.census3A00) == "dbh.x"] <- "dbh"

#############
load("yasuni.census3B.Rdata")
load("yasuni.census4A.Rdata")   
yasuni.census4A00 <- merge(yasuni.census4A, yasuni.census3B[, c("tag", "sp", "dbh")], by = c("tag", "sp"), all.x = TRUE)
# Completa los valores faltantes de dbh en yasuni.census2A con los valores de dbh de yasuni.census1
yasuni.census4A00$dbh.x[is.na(yasuni.census4A00$dbh.x)] <- yasuni.census4A00$dbh.y[is.na(yasuni.census4A00$dbh.x)]
# Elimina la columna dbh.y y renombra la columna dbh.x si es necesario
yasuni.census4A00 <- yasuni.census4A00[, !(names(yasuni.census4A00) %in% c("dbh.y"))]
colnames(yasuni.census4A00)[colnames(yasuni.census4A00) == "dbh.x"] <- "dbh"



yasuni.census2A00$period<-"first"
yasuni.census3A00$period<-"second"
yasuni.census4A00$period<-"third"

yasuni.census2A00$dry<-1
yasuni.census3A00$dry<-0
yasuni.census4A00$dry<-1

data<-rbind(yasuni.census2A00, yasuni.census3A00,yasuni.census4A00)

data <- data %>%
  mutate(status = ifelse(status %in% c("A", "AB", "ABD", "AD"), 0, 1))

sppnames<-read.csv("yas_spp_traits2.csv")
# Merge para agregar la columna 'especie' a df_grande
data <- merge(data, sppnames[, c("sp", "species00", "wsg")], by = "sp", all.x = TRUE)

#data <- read.csv("ANALYSES_TOTAL_VANCOUVER_2019_2_T1-T3.csv")
#data<-data[c(4,5,7,12:16, 18,20,21,22,23)]
#data<-data[c(4,5,6,7,13,14,22,24,26,27,28)]
#data<-data[c(4,5,6,7,13,14,22,23,26,27,28,29,30,31,32,33,34,35,36,37,38,39)]
#data<-data[!(data$SPECIES=="Pristimantis_w-nigrum"),]#no w-nigrum
#data<-data[!(data$SPECIES=="Rhinella_margaritifer"),]#no r-margaritifer, when I
#exclude R_margaritifer of analyses, i can inclde locality as random factor, but in that case
#I have p_change non singnificant.

#names(data) <- tolower(names(data))

#Fill all "blanks cells" to NA's
#data[data==""] <- NA
#Eliminate NA's
#data<-data[complete.cases(data),]

#data<-data[!(data$locality=="CASHCA"),]

###############
###########################
###########################
tree <- read.tree("myspecies997_in_Qiang.tre")
tree <- drop.tip(tree, setdiff(tree$tip.label,as.character(data$species00)))

# some species missing...
setdiff(data$species, tree$tip.label)#no species missing???

# remove these rows in data
data <- data[data$species00%in%tree$tip.label,]

# one row has NA in response - remove
#sum(is.na(data$ind_mite_load))

#data <- data[!is.na(data$class_response),]

# phylogenetic correlation structure
phy_cov <- vcv(tree, corr=TRUE)

ggplot(gather(data%>% purrr::keep(is.numeric)), aes(value)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~key, scales = 'free')#shouls I transform volumes00????

# looks that in terms of predictors, lake_sa could be log transformed
# perhaps sqrt transform would make more sense, but I will log for now
#data$lake_sa <- log(data$lake_sa)


##########################################################
#THIS IS THE FINAL MODEL FOR CHAPTER THREE AFTER MANY
#ANALYSES, INCLUDING THOSE WITH J. DAVIES IN VANCOUVER AND MCGILL
##INCLUDES NO W-NIGRUM!
#THIS IS THE SAME as brms MODEL
#I HAVE TO EXPLAIN WHY I DID NOT INCLUDE LOCALITY AS RANDOM
#FACTOR, I SHOULD SAY THAT I ALREADY INCLUDE ALTITUDE AND 
#P_CHANGE AND T INTENSITY
#RESULTS ARE SIMILAR TO GLMER, EXCEPT FOR INTENSITY
#MANY OF THE RUNNING MODELS ARE IN THE FILES:
#"NO INTENSITY MODELS.docx" AND "baysian models.docx"
#TALK ABOUT WWHY INTENSITY IS DEGATIVE...THAT MAYBE 
#P_CHANGE SET OF BLA BLA

# Establecer priors básicos
#priors <- c(
 # prior(normal(0, 10), class = Intercept),  # Prior para el intercepto
  #prior(student_t(3, 0, 2.5), class = sd),  # Prior para la desviación estándar
  #prior(student_t(3, 0, 2.5), class = sd, coef = "Intercept"),  # Prior para la desviación estándar del intercepto
  #prior(student_t(3, 0, 2.5), class = b),  # Priors para los coeficientes b
  #prior(student_t(3, 0, 2.5), class = b, coef = "habsecondary"),  # Prior para el coeficiente de habsecondary
  #prior(student_t(3, 0, 2.5), class = b, coef = "habvalley"),  # Prior para el coeficiente de habvalley
  #prior(student_t(3, 0, 2.5), class = b, coef = "period"),  # Prior para el coeficiente de period
  #prior(student_t(3, 0, 2.5), class = b, coef = "scaledbh")  # Prior para el coeficiente de scaledbh
#)


# Define priors
  priors <- c(
    set_prior("normal(0, 1)", class = "b"),  # Normal prior for all regression coefficients
    set_prior("student_t(3, 0, 2)", class = "Intercept"),  # Specific prior for the Intercept
    # Remove the global sd prior and specify only where needed
    set_prior("student_t(3, 0, 2)", class = "sd", coef = "Intercept", group = "q20"),
    set_prior("student_t(3, 0, 2)", class = "sd", coef = "Intercept", group = "species00"),
    set_prior("student_t(3, 0, 2)", class = "sd", coef = "Intercept", group = "tag"))
  
  

# Muestrear aleatoriamente el 1% de los datos
#data <- data %>% sample_frac(size = 0.01)

#BEST
fit_mortality <- brm(status ~ scale(dbh) + hab  + scale(wsg) +
                     (1|tag) + (1|q20) + dry +
                     (1|gr(species00, cov = phy_cov)), #ver si incluyo esto!
                     data=data, family=bernoulli(), #bernoulli
                     data2 = list(phy_cov=phy_cov),
                     prior=priors,
                     iter=2000, 
                     warmup = 1000,cores = 4,
                     control=list(adapt_delta=0.95, max_treedepth = 10))



summary(fit_mortality)
bayes_R2(fit_mortality)
library(sjstats)
tidy_stan(fit_mortality)#calculate highest posterior density interval
#tidy_stan(fit_extinction, prob = 0.95)#other option with other probability
r2(fit_mortality)
saveRDS(fit_extinction, "extinction.rds")###see what it does.

plot(marginal_effects(fit_mortality, mean = TRUE , spaghetti = TRUE, nsamples = 100))#, ask = TRUE)

marginal_effects(fit_mortality, "niche.unfit", nsamples=1000, mean = TRUE,spaghetti = TRUE)
plot(me, plot = FALSE)[[1]] + 
  scale_color_grey() +
  scale_fill_grey()
launch_shinystan(fit_mortality)
######################################################
#########################################################
#################################################
####################################################################################
#BEST
#fit_extinction <- brm(response ~ scale(volumes00) + scale(distance_volume00) + 
# scale(p_change) +  scale(ed_equal) + scale(change_tmax_before_mean) 
#  + scale(altitude) +
#   (1|species), 
#  dat=data, family=bernoulli(), #bernoulli
# iter=2000, cov_ranef = list(species = phy_cov),
# control=list(adapt_delta=0.98))#good, except for the sign "T intensity"
#esto sin w-nigrum...
#summary(fit_extinction)
#bayes_R2(fit_extinction)
#tidy_stan(fit_extinction)#calculate highest posterior density interval
#saveRDS(fit_extinction, "extinction.rds")###see what it does.

#plot(marginal_effects(fit_extinction), ask = TRUE)
########################################################
###BEST EVER!!!!!!BEST EVER#####BEST####THE BEST##############WITH JONATHAN

#fit_extinction<- brm (response ~ scale(volumes00)  +  scale(distance_volume00)
#+ (1|locality) + (1|species) + scale(p_change)  
#+ scale(change_tmax_before_mean) + scale(altitude) + scale(ed_equal), 
#   dat=data, family=bernoulli(), #bernoulli
#  iter=2000, cov_ranef = list(species = phy_cov),
# control=list(adapt_delta=0.98, max_treedepth = 15))#THESE PARATMETERS CHANGED

# summary(fit_extinction)

####TRY THIS ONE WITH W-NIGRUM
#with altitude?
#with ed?.....NO ESTA SALIENDO COMO CUANDO CORRI CON JONATHAN

#################################################################################3

#by using the cov_ranef argument, we make sure that species are correlated as specified by the covariance matrix A
#saveRDS(fit_extinction, "extinction_.rds")###see what it does.

#summary(fit_extinction)
#library(sjstats)
#tidy_stan(fit_extinction)#calculate highest posterior density interval
#r2(fit_extinction)#computes either the Bayes 
#r-squared value, or - if loo = TRUE - a LOO-adjusted
#r-squared value (which comes conceptionally closer 
#to an adjusted r-squared measure).


#Bayesian version of R2
#bayes_R2(fit_extinction)
#plot(fit_extinction)

############################Posterior predictive checks############################
#Posterior predictive checks
# first make posterior predictions
pp_extinct <- brms::posterior_predict(fit_extinction)
#Density overlay
# define log+1 transformation for visualizations
log1 <- scale_x_continuous(trans="log1p")

ppc_dens_overlay(y=data$response, pp_extinct[1:50, ]) + log1 + 
  coord_cartesian(xlim = c(0, 100))#change this value?

#SOME LINES OF MAX HTML FILE ARE MISSINIG: see the html file


#Posterior Estimate of Phylogenetic Signal
source("icc.R")
#fit_extint_icc <- icc(fit_extinction, is_negbin=TRUE)
fit_extinct_icc00 <- icc(fit_extinction, is.binary=TRUE)
plot(density((fit_extinct_icc00$icc_species)), main="")#DOES NOT WORK
abline(v=median(fit_nb_icc$icc_genus_species), col=4, lty=3)#DOES NOT WORK
abline(v=quantile(fit_nb_icc$icc_genus_species, c(0.025,0.975)),col=2, lty=3)#DOES NOT WORK
#In short, the ICC can be interpreted as "the proportion of the variance explained by the grouping
#structure in the population" (Hox 2002: 15).



#Plot of species effects on the phylogeny
# Extract species level effects
ranef_spec <- ranef(fit_extinction)$species
ranef_spec <- as.data.frame(unlist(ranef_spec))

# re-order efects to match tree
ranef_spec <- ranef_spec[tree$tip.label,]

# Plot the beside the phylogeny# DOES NOT WORK THE GRAPH
cols <- c("#7570b3","#d95f02")
plot(tree,type="p",TRUE, cex=1,no.margin=FALSE, show.tip.label = TRUE, label.offset=10)
tiplabels(pch=22,bg=cols[ifelse(ranef_spec$Estimate.Intercept>0,1,2)],col="black",cex=abs(ranef_spec$Estimate.Intercept),adj=1.5)
legend("topright",legend=c("2","1","0","-1","-2")),pch=21,
pt.bg=cols[c(1,1,1,2,2)],bty="n",
text.col="gray32",cex=1, pt.cex=c(2,1,0.1,1,2))

########################################################################
#other ways to plot tarits on a phylogenetic tree#

library(phytools)
x<-ranef_spec[1]

#The simplest class of visualization method simply plots the data at the tips of the tree. For instance:
dotTree(tree,x,length=10,ftype="i")

plotTree.barplot(tree,x)

x1<-as.matrix(x)
plotTree.barplot(tree,x1,args.barplot=list(col=sapply(x1,
                                                      function(x1) if(x1>0) "blue" else "red"),
                                           xlim=c(-20,15)))








