##########################
# POP77054 Final Project #
##########################

# Preparing the environment
rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Set options
options(scipen = 100)

# Loading required packages
library(tidyverse)
library(psych)
library(GDAtools)
library(dplyr)

# Loading the data to the environment
data <- read_csv("ISSP_Religion.csv")

# Subsetting the data to keep only the specified variables
data <- data[, c("v1", "v2", "v3", "v4", "v5", "v6", "v17", "v20", "v21", "v22", "v23", "v24", "v25", "v26", 
                 "AGE", "DEGREE", "MARITAL", "SEX", "URBRURAL", "country"), drop = FALSE]

# Selecting only EU countries
data_eu <- data[!data$country %in% c("152. CL-Chile", "158. TW-Taiwan", "352. IS-Iceland", "376. IL-Israel", "392. JP-Japan", "410. KR-Korea (South)", "554. NZ-New Zealand", "578. NO-Norway", "608. PH-Philippines", "643. RU-Russia", "710. ZA-South Africa", "740. SR-Suriname", "756. CH-Switzerland", "764. TH-Thailand", "792. TR-Turkey", "826. GB-Great Britain and/or United Kingdom", "840. US-United States"), ]
data_eu <- na.omit(data_eu)

# Getting rid of numerics, punctuations and apostophe from value labels
data_eu <- as.data.frame(lapply(data_eu, function(x) gsub("^\\d+\\.\\s*", "", x)))
data_eu <- as.data.frame(lapply(data_eu, function(x) gsub("'", "", x)))

# Converting age variable to numeric and then into a categorical variable
data_eu$AGE <- as.numeric(as.character(data_eu$AGE))
na_index <- which(is.na(data_eu$AGE))
data_eu <- data_eu[-na_index,]
data_eu$AGE <- cut(data_eu$AGE, 
                breaks = c(15, 35, 50, 65, Inf),
                labels = c("15-34", "35-49", "50-64", "65+"),
                include.lowest = TRUE)

# Renaming variables for easier reference
names(data_eu) <- c(
  "happy",            # v1, Respondent is happy
  "satisfied",        # v2, Satisfied with family life
  "extra_m",          # v3, Extramarital relations
  "same_sex",         # v4, Same-sex relation
  "abort",            # v5, Abortion views
  "gender_roles",     # v6, Gender roles
  "diff_rel",         # v17, Accept someone from other religion
  "god_belief",       # v20, Extent of beliefs about God
  "god_belief_pastnow",# v21, Beliefs about God in past and now
  "belief_afterlife", # v22, Belief in afterlife
  "heaven",           # v23, Belief in heaven
  "hell",             # v24, Belief in hell
  "miracles",         # v25, Belief in miracles
  "supernatural_powers", # v26, Belief in supernatural powers
  "age",              # AGE, Age of the respondent
  "degree",           # DEGREE, Degree of the respondent
  "marital",          # MARITAL, Marital status of the respondent
  "sex",              # SEX, Sex of the respondent
  "urbrural",         # URBRURAL, Urban-rural location of the respondent
  "country"           # country, Country of the respondent
)

# Exploring the reduced dataset
describe(data_eu)

## MCA
# define active variables
active.variables <- c(
  "happy", # v1, Respondent is happy
  "same_sex", # v4, Same-sex relation
  "abort", # v5, Abortion views
  "gender_roles", # v6, Gender roles
  "god_belief", # v20, Extent of beliefs about God
  "belief_afterlife", # v22, Belief in afterlife
  "heaven", # v23, Belief in heaven
  "hell", # v24, Belief in hell
  "miracles", # v25, Belief in miracles
  "supernatural_powers" # v26, Belief in supernatural powers
)

# examine frequency of categories in active variables
lapply(data_eu[,active.variables], ftable) %>% # get a flat table for each var
  lapply(., prop.table) %>% # calculate proportions
  lapply(., round, 2) # round the proportions to 2 decimal places

## passive categories (specific MCA)
# create vector of categories to include as passive (specific MCA)
passive.categories <- c("happy.Cant choose", "happy.No answer; GE: Refused",
                        "same_sex.Cant choose", "same_sex.No answer; GE: Refused",
                        "abort.Cant choose", "abort.No answer; GE: Refused",
                        "gender_roles.Cant choose", "gender_roles.No answer; GE: Refused",
                        "god_belief.Dont know", "god_belief.No answer",
                        "belief_afterlife.Cant choose", "belief_afterlife.No answer; GE: Refused",
                        "heaven.Cant choose", "heaven.No answer; GE: Refused",
                        "hell.Cant choose", "hell.No answer; GE: Refused",
                        "miracles.Cant choose", "miracles.No answer; GE: Refused",
                        "supernatural_powers.Cant choose", "supernatural_powers.No answer; GE: Refused")

## perform MCAspe
mca.spe <- speMCA(data_eu[,active.variables], excl = passive.categories)

## summarise results
# extract eigenvalues
e.vals <- lapply(mca.spe$eig, round, 2) # we round the eigenvalue output to 2 d.p.
data.frame(eigen = e.vals$eigen[1:10], # we extract the first 10 eigenvalues
           mod.rate = e.vals$mrate,
           cumul.rate = e.vals$cum.mrate)

# extract variable statistics
dimdescr(mca.spe, 
         dim = c(1,2,3)
         ) # this provides us with summary statistics for all vars and cats

# find contributions (and critical thresholds)
dimcontrib(mca.spe, c(1:3))$var # note, this output provides negative contributions

# to calculate absolute values of contribution first, extract the data for variables
dim.contrib <- dimcontrib(mca.spe, c(1:3), best = FALSE)$var
# second, convert to absolute values (remove negative signs)
dim.contrib.abs <- lapply(dim.contrib, abs)
# third (and finally) arrange the variables in descending order by dimension.
# note: this requires us to use a lambda (anonymous) function inside lapply()
dim.contrib.abs <- lapply(dim.contrib.abs, function(x) arrange(x, desc(ctr)))
dim.contrib.abs

# calculate cutoff
K <- length(mca.spe$call$marge.col) # number of active categories
c.val <- 1/K * 100

# apply cutoff to contributions object
lapply(dim.contrib.abs, function(x) x[x$ctr >= c.val,])

# we can also calculate contributions by plane (in case of Guttman effect, etc.)
# what cut-off would you apply? 2*K?
plane1.contrib <- planecontrib(mca.spe, axes = c(1,2))$var 
sort(plane1.contrib$ctr, decreasing = TRUE)

plane2.contrib <- planecontrib(mca.spe, axes = c(1,3))$var
sort(plane2.contrib$ctr, decreasing = TRUE)

#####################################
#### GDA: plotting with GDAtools ####
#####################################

# GDAtools has markedly better integration with ggplot2, as well as improved
# functions for plotting in general.

## Plotting basic clouds using GDAtools and ggplot2
# cloud of variables
ggcloud_variables(mca.spe, 
                  axes = c(1,2), 
                  points = "best") +
  ggtitle("First factorial plane, top contributions")

ggcloud_variables(mca.spe, 
                  axes = c(1,3), 
                  points = "best") +
  ggtitle("Second factorial plane, top contributions")

# cloud of individuals
ggcloud_indiv(mca.spe,
              axes = c(1,2)) +
  ggtitle("Cloud of individuals, first factorial plane")

ggcloud_indiv(mca.spe,
              axes = c(1,3)) +
  ggtitle("Cloud of individuals, second factorial plane")

## Adding supplementary variables
# Because we are plotting with ggplot2, we can build up our plots one layer at a time:

# Preparing supplementary variables
data_eu$sex <- factor(data_eu$sex)
# Remove "no answer" category
data_eu$sex <- droplevels(data_eu$sex, exclude = "No answer")
levels(data_eu$sex)
data_eu$age <- factor(data_eu$age)

data_eu$degree <- factor(data_eu$degree)
levels(data_eu$degree)
data_eu$degree <- droplevels(data_eu$degree, exclude = c("GE: Dont know", "No answer"))
# Recode the levels of the degree variable
data_eu$degree <- recode_factor(data_eu$degree,
                                "No formal education" = "No formal education",
                                "Lower secondary (secondary completed that does not allow entry to university: end of obligatory school)" = "Lower education",
                                "Primary school (elementary education)" = "Lower education",
                                "Upper secondary (programs that allow entry to university)" = "Secondary education",
                                "Lower level tertiary, first stage (also technical schools at a tertiary level)" = "Higher education",
                                "Post secondary, non-tertiary (other upper secondary programs toward the labour market or technical formation)" = "Higher education",
                                "Upper level tertiary (Master, Doctor)" = "Higher education")

data_eu$marital <- factor(data_eu$marital)
levels(data_eu$marital)
data_eu$marital <- droplevels(data_eu$marital, exclude = c("No answer", "Refused"))
# Recode the levels of the marital status variable
data_eu$marital <- recode_factor(data_eu$marital,
                                 "Civil partnership; IS: In registered cohabition" = "Married or in civil partnership",
                                 "Divorced from spouse/ legally separated from civil partner; AT: Div./Sep." = "Divorced or separated",
                                 "Married" = "Married or in civil partnership",
                                 "Never married/ never in a civil partnership, single" = "Never married or single",
                                 "Separated from spouse/ civil partner (still legally married/ still legally in a civil partnership)" = "Divorced or separated",
                                 "Widowed/ civil partner died" = "Married or in civil partnership")

data_eu$urbrural <- factor(data_eu$urbrural)
levels(data_eu$urbrural)
data_eu$urbrural <- droplevels(data_eu$urbrural, exclude = c("No answer; SR: Not asked to everyone"))
# Recode the levels of the urbrural variable
data_eu$urbrural <- recode_factor(data_eu$urbrural,
                                  "A big city; ZA: Urban, Formal" = "Urban",
                                  "A country village; ZA: Trad. Auth Areas" = "Rural",
                                  "A farm or home in the country; ZA: Farms" = "Rural",
                                  "A town or a small city" = "Rural",
                                  "The suburbs or outskirts of a big city; ZA: Urban, Informal" = "Urban")


# plot the cloud of variables
var.cloud <- ggcloud_variables(mca.spe, 
                               axes = c(1,2), 
                               points = "best")

# add a supplementary variable (sex)
ggadd_supvar(var.cloud, 
             mca.spe, 
             data_eu$sex,
             axes = c(1,2))

# add a supplementary variable (location)
ggadd_supvar(var.cloud, 
             mca.spe, 
             data_eu$urbrural,
             axes = c(1,2))

# add a supplementary variable (marital)
ggadd_supvar(var.cloud, 
             mca.spe, 
             data_eu$marital,
             axes = c(1,2))

# Adding lines to link ordered supplementary variables
ggadd_supvar(var.cloud,
               mca.spe,
               data_eu$age,
               segment = T)

ggadd_supvar(var.cloud,
             mca.spe,
             data_eu$degree,
             segment = T)

## Plotting ellipses
# There are two different functions depending on whether you wish to plot an 
# indicator/concentration ellipse, or a confidence ellipse.

# confidence ellipse
ind.cloud <- ggcloud_indiv(mca.spe,
                           axes = c(1,2),
                           col = "grey") 

ggadd_ellipses(ind.cloud,
               mca.spe,
               data_eu$sex,
               axes = c(1,2),
               level = 0.05, # for 95 confidence level
               legend = "none",
               size = 1,
               points = F # it is very difficult to see ellipses with coloured points
               ) +
  ggtitle("Confidence ellipse for sex, 1st factorial plane")


ggadd_ellipses(ind.cloud,
               mca.spe,
               data_eu$urbrural,
               axes = c(1,2),
               level = 0.05, # for 95 confidence level
               legend = "none",
               size = 1,
               points = F # it is very difficult to see ellipses with coloured points
) +
  ggtitle("Confidence ellipse for location, 1st factorial plane")


ggadd_ellipses(ind.cloud,
               mca.spe,
               data_eu$degree,
               axes = c(1,2),
               sel = c(1,4), # no formal and higher degree categories selected
               level = 0.05, # for 95 confidence level
               legend = "none",
               size = 1,
               points = F # it is very difficult to see ellipses with coloured points
) +
  ggtitle("Confidence ellipse for education, 1st factorial plane")


ggadd_ellipses(ind.cloud,
               mca.spe,
               data_eu$marital,
               axes = c(1,2),
               sel = c(1,3), # no formal and higher degree categories selected
               level = 0.05, # for 95 confidence level
               legend = "none",
               size = 1,
               points = F # it is very difficult to see ellipses with coloured points
) +
  ggtitle("Confidence ellipse for marital status, 1st factorial plane")


ggadd_ellipses(ind.cloud,
               mca.spe,
               data_eu$age,
               axes = c(1,2),
               sel = c(1,4), # no formal and higher degree categories selected
               level = 0.05, # for 95 confidence level
               legend = "none",
               size = 1,
               points = F # it is very difficult to see ellipses with coloured points
) +
  ggtitle("Confidence ellipse for age, 1st factorial plane")
