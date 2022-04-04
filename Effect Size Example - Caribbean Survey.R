#==============================================
# EFFECT SIZE EXAMPLE
# Caribbean Meadow Survey Data
#==============================================

# DATA SOURCE: Johnson et al. (2019) Journal of Ecology
# Manuscript: https://besjournals-onlinelibrary-wiley-com.ezproxy.library.wisc.edu/doi/full/10.1111/1365-2745.13306
# Data Source & Metadata as EDI package 422.1
#---------------------------
# Script written by G.M. Wilkinson (March 2022) for ZOO 955 as an example of meta-analysis approaches. Analysis intended for demonstration purposes only
# DO NOT USE FOR STATISTICAL INFERENCE - consult original manuscript and analyses therein for scientific conclusions from this literature synthesis

#==============================================
# PACKAGES
if (!require(tidyverse)) install.packages('tidyverse')
library(tidyverse)
if (!require(effectsize)) install.packages('effectsize')
library(effectsize)
if (!require(lme4)) install.packages('lme4')
library(lme4)
#==============================================
# DATA DOWLOAD

# Package ID: edi.422.1 Cataloging System:https://pasta.edirepository.org.
# Contact:  Robert Johnson -  University of Florida  - johnson.robert@ufl.edu
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl2  <- "https://pasta.lternet.edu/package/data/eml/edi/422/1/d035311a5857aa3f958ccb26fca97869" 
infile2 <- tempfile()
try(download.file(inUrl2,infile2,method="curl"))
if (is.na(file.size(infile2))) download.file(inUrl2,infile2,method="auto")


dt2 <-read.csv(infile2, header = F, skip = 1, sep = ",", quot = '"', check.names=TRUE,
               col.names = c("Location", "Site", "Species", "Treatment", "Date", "GPP", "RE", "NEP"))
unlink(infile2)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings
if (class(dt2$Location)!="factor") dt2$Location<- as.factor(dt2$Location)
if (class(dt2$Site)!="factor") dt2$Site<- as.factor(dt2$Site)
if (class(dt2$Species)!="factor") dt2$Species<- as.factor(dt2$Species)
if (class(dt2$Treatment)!="factor") dt2$Treatment<- as.factor(dt2$Treatment)                                   
# attempting to convert dt2$Date dateTime string to R date structure (date or POSIXct)   
tmpDateFormat<-"%Y-%m-%d"
tmp2Date<-as.Date(dt2$Date,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp2Date) == length(tmp2Date[!is.na(tmp2Date)])){dt2$Date <- tmp2Date } else {print("Date conversion failed for dt2$Date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp2Date) 
if (class(dt2$GPP)=="factor") dt2$GPP <-as.numeric(levels(dt2$GPP))[as.integer(dt2$GPP) ] 
if (class(dt2$GPP)=="character") dt2$GPP <-as.numeric(dt2$GPP)
if (class(dt2$RE)=="factor") dt2$RE <-as.numeric(levels(dt2$RE))[as.integer(dt2$RE) ]     
if (class(dt2$RE)=="character") dt2$RE <-as.numeric(dt2$RE)
if (class(dt2$NEP)=="factor") dt2$NEP <-as.numeric(levels(dt2$NEP))[as.integer(dt2$NEP) ] 
if (class(dt2$NEP)=="character") dt2$NEP <-as.numeric(dt2$NEP)

# Convert Missing Values to NA for non-dates
dt2$GPP <- ifelse((trimws(as.character(dt2$GPP))==trimws("NA")),NA,dt2$GPP)               
suppressWarnings(dt2$GPP <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt2$GPP))==as.character(as.numeric("NA"))),NA,dt2$GPP))
dt2$RE <- ifelse((trimws(as.character(dt2$RE))==trimws("NA")),NA,dt2$RE)               
suppressWarnings(dt2$RE <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt2$RE))==as.character(as.numeric("NA"))),NA,dt2$RE))
dt2$NEP <- ifelse((trimws(as.character(dt2$NEP))==trimws("NA")),NA,dt2$NEP)               
suppressWarnings(dt2$NEP <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt2$NEP))==as.character(as.numeric("NA"))),NA,dt2$NEP))

#==============================================
# Filter to only the data we need    
surveyData = dt2 %>% 
  #Remove the invasive species treatment category
  filter(!(Treatment=="invasive")) %>%
  #Remove excess columns that we don't need for analysis
  select(-Species, -Date) %>%
  #rename the columns with lower case because I am a sloppy typist
  rename(location = Location,
         treatment = Treatment) %>%
  #Replace the Bahamas island name with nation
  mutate(location = as.character(location),
         location = replace(location, location=="Eleuthera, Bahamas", "Bahamas"))

#==============================================
# Hedge's g - Effect Size Estimate
# Unpaired samples, small sample size

#Bahamas ==========================
bahamas = surveyData %>%
  filter(location == "Bahamas")
#Hedges g (correction for small sample size)
bahamas_g <- hedges_g(NEP ~ treatment, data = bahamas,
                      pooled_sd = TRUE, paired = FALSE)
interpret_hedges_g(bahamas_g, rules = "cohen1988")
#Variance estimate for Hedge's g
bahamas_g$variance = 
  (  (length(bahamas$treatment=="grazed") + 
        length(bahamas$treatment=="ungrazed")) / 
       (length(bahamas$treatment=="grazed") * 
          length(bahamas$treatment=="ungrazed"))) + 
  bahamas_g$Hedges_g^2/(2*((length(bahamas$treatment=="grazed") + 
                              length(bahamas$treatment=="ungrazed"))))
#Cohens d
bahamas_d <- cohens_d(NEP ~ treatment, data = bahamas,
                      pooled_sd = TRUE, paired = FALSE)

#Bonaire ==========================
bonaire = surveyData %>%
  filter(location == "Bonaire")
bonaire_g <- hedges_g(NEP ~ treatment, data = bonaire,
                      pooled_sd = TRUE, paired = FALSE)
interpret_hedges_g(bonaire_g, rules = "cohen1988")
#Variance estimate for Hedge's g
bonaire_g$variance = 
  (  (length(bonaire$treatment=="grazed") + 
        length(bonaire$treatment=="ungrazed")) / 
       (length(bonaire$treatment=="grazed") * 
          length(bonaire$treatment=="ungrazed"))) + 
  bonaire_g$Hedges_g^2/(2*((length(bonaire$treatment=="grazed") + 
                              length(bonaire$treatment=="ungrazed"))))
#Cohens d
bonaire_d <- cohens_d(NEP ~ treatment, data = bonaire,
                      pooled_sd = TRUE, paired = FALSE)

#Florida ==========================
florida = surveyData %>%
  filter(location == "Florida")
florida_g <- hedges_g(NEP ~ treatment, data = florida,
                      pooled_sd = TRUE, paired = FALSE)
interpret_hedges_g(florida_g, rules = "cohen1988")
#Variance estimate for Hedge's g
florida_g$variance = 
  (  (length(florida$treatment=="grazed") + 
        length(florida$treatment=="ungrazed")) / 
       (length(florida$treatment=="grazed") * 
          length(florida$treatment=="ungrazed"))) + 
  florida_g$Hedges_g^2/(2*((length(florida$treatment=="grazed") + 
                                length(florida$treatment=="ungrazed"))))

#Cohens d
florida_d <- cohens_d(NEP ~ treatment, data = florida,
                      pooled_sd = TRUE, paired = FALSE)

#Little Cayman ====================
lilcayman = surveyData %>%
  filter(location == "Little Cayman")
lilcayman_g <- hedges_g(NEP ~ treatment, 
                      data = lilcayman,
                      pooled_sd = TRUE,
                      paired = FALSE)
interpret_hedges_g(lilcayman_g, rules = "cohen1988")
#Variance estimate for Hedge's g
lilcayman_g$variance = 
  (  (length(lilcayman$treatment=="grazed") + 
        length(lilcayman$treatment=="ungrazed")) / 
       (length(lilcayman$treatment=="grazed") * 
          length(lilcayman$treatment=="ungrazed"))) + 
  lilcayman_g$Hedges_g^2/(2*((length(lilcayman$treatment=="grazed") + 
                              length(lilcayman$treatment=="ungrazed"))))

#Cohens d
lilcayman_d <- cohens_d(NEP ~ treatment, 
                      data = lilcayman,
                      pooled_sd = TRUE,
                      paired = FALSE)

#St. Croix ========================
stcroix = surveyData %>%
  filter(location == "St. Croix")
stcroix_g <- hedges_g(NEP ~ treatment, 
                        data = stcroix,
                        pooled_sd = TRUE,
                        paired = FALSE)
interpret_hedges_g(stcroix_g, rules = "cohen1988")
#Variance estimate for Hedge's g
stcroix_g$variance = 
(  (length(stcroix$treatment=="grazed") + 
     length(stcroix$treatment=="ungrazed")) / 
  (length(stcroix$treatment=="grazed") * 
     length(stcroix$treatment=="ungrazed"))) + 
  stcroix_g$Hedges_g^2/(2*((length(stcroix$treatment=="grazed") + 
                              length(stcroix$treatment=="ungrazed"))))

#Cohens d
stcroix_d <- cohens_d(NEP ~ treatment, 
                        data = stcroix,
                        pooled_sd = TRUE,
                        paired = FALSE)

#Combine effect size estimates into one data frame
#Add the site name to the data frame first
bahamas_g$site = "bahamas"; bahamas_d$site = "bahamas"
bonaire_g$site = "bonaire"; bonaire_d$site = "bonaire"
florida_g$site = "florida"; florida_d$site = "florida"
lilcayman_g$site = "lilcayman"; lilcayman_d$site = "lilcayman"
stcroix_g$site = "stcroix"; stcroix_d$site = "stcroix"

hedges_df = rbind(bahamas_g, bonaire_g, florida_g, lilcayman_g, stcroix_g)
cohens_df = rbind(bahamas_d, bonaire_d, florida_d, lilcayman_d, stcroix_d)

# Plot the Effect Sizes ==================================
#HEDGES G
plot(rank(hedges_df$Hedges_g), hedges_df$Hedges_g, 
     ylim = c(-14,0), xlim = c(0.5,5.5), xaxt = "n",
     pch = 15, col = "dodgerblue3", cex = 1.5,
     xlab = "", ylab = "Effect Size (Grazed vs Ungrazed NEP0")
arrows(rank(hedges_df$Hedges_g), hedges_df$CI_low,
       x1 = rank(hedges_df$Hedges_g), y1 = hedges_df$CI_high,
       length = 0, angle = 0, lwd = 2, col = "dodgerblue3")
#COHEN'S D
points(rank(hedges_df$Hedges_g)+0.2, cohens_df$Cohens_d,
       pch = 19, col = "turquoise3", cex = 1.5)
arrows(rank(hedges_df$Hedges_g)+0.2, cohens_df$CI_low,
       x1 = rank(hedges_df$Hedges_g)+0.2, y1 = cohens_df$CI_high,
       length = 0, angle = 0, lwd = 2, col = "turquoise3")
#Arts and Crafts
axis(side = 1, at = c(1,2,3,4,5),
     label = c("Florida", "Caymans","Bonaire", "St. Croix", "Bahamas"))
legend("bottomright", legend = c("Hedge's g", "Cohen's d"), 
       col = c("dodgerblue3", "turquoise3"), pch = c(15,19), pt.cex = 2,
       inset = 0.05)

#=============================================
# FIXED EFFECTS MODEL
# Calcuate the weights for the fixed effects model
hedges_df$weights <- 1/hedges_df$variance^2

mean_effect = 
  sum(hedges_df$weights * hedges_df$Hedges_g) / sum(hedges_df$weights)

mean_var = 
  1 / sum(hedges_df$weights)

lines(c(-1,10), c(mean_effect, mean_effect), lwd = 2, lty = 3)
text(5, mean_effect-0.5 ,"Mean Effect", font = 2)
      