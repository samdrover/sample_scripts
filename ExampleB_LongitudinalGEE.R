## >FILE INFO ############################################################################

# AUTHOR:               Samantha Drover
# PURPOSE:              Run all analyses for insecticide-immune paper (by tables/figures) using GEE combining across times
# CREATED:              2019-06-19
# LAST MODIFIED:        2019-08-13


# DATA DICTIONARY:      
# [REDACTED DIRECTORY TO DATA DICTIONARY]


# HELPFUL LINKS:
    # http://ftp.uni-bayreuth.de/math/statlib/R/CRAN/doc/vignettes/Zelig/logit.gee.pdf (also saved in zotero )



## >PATHS & LIBRARIES #############################################################################

rm(list=ls()) #clears workspace

library(DescTools)
library(EnvStats) #for geometric mean and standard deviation
library(ggplot2)
library(foreign)
library(lubridate)
library(rms)
library(FSA) #for standard error
library(ghibli)
library(contrast)
library(scales)
library(stats)
library(extrafont)
library(gee)
library(xlsx)
library(doBy) #for getting predicted values from GEE
library(boot) #for booting CIs on predicted values from GEE
library(Zelig) #for getting predicted values from GEE models 
devtools::install_github("awalker89/openxlsx")

Sys.setenv(TZ = "Africa/Johannesburg") #setting timezone so that standard setting will not apply daylight savings time to date/time values
base::Sys.timezone()



#extrafont::font_import() #for importing fonts - takes a while
extrafont::loadfonts(device="win")
extrafont::fonttable()

## Reading in the updated dataframe with the correct pyrethroid variables, as we as additional vars needed ##
df = read.csv("[REDACTED PATH TO INSECTICIDE EXPOSURE/IMMUNE ANTIBODY TITRE DATA]")
labels = read.csv("[REDACTED PATH TO LABELS FOR  ABOVE DATA FRAME]")


## Reading in file with LOD variable ##
lod = foreign::read.dta("[REDACTED PATH TO LIMIT OF DETECTION (LOD) DATA]")
colnames(lod) = tolower(colnames(lod))


## Reading in the 5-year antibody levels ##
yr5 = read.xlsx(file = "[REDACTED DIRECTORY TO ANTIBODY TITRES AT 5 YEARS]", 
                sheetName ="Results")
colnames(yr5) = tolower(colnames(yr5))
colnames(yr5) = c("hsn", paste(colnames(yr5)[2:4], "_5yr", sep=""))


## Reading in the SCHEDULED date data for the 5-year blood collection --THESE ARE THE DATES TO USE ##
blooddate = read.xlsx(file = "[REDACTED DIRECTORY TO DATE OF BLOOD COLLECTION DATA]",
                      sheetIndex = 1)
colnames(blooddate) = tolower(colnames(blooddate))



## Reading in the date data for the 5-year visit ##
visitdate = read.xlsx(file = "[REDACTED DIRECTORY TO VISIT DATE DATA]",
                      sheetIndex = 1)
colnames(visitdate) = tolower(colnames(visitdate))


## Setting working directory for output files ##
setwd("[REDACTED WORKING DIRECTORY]")


## >SET DATA ################################################################################

## There is a double ID somewhere -- making a function to detect where it is ##
doubleID = function(data) {
  
  t = as.data.frame(table(data$hsn))
  id = t$Var1[t$Freq>=2]
  return(id)
  
}

doubleID(blooddate) #1105 is doubled here
doubleID(visitdate) #no
doubleID(df) #no
doubleID(lod) #no
doubleID(visitdate) #no
doubleID(yr5) #no

blooddate[blooddate$hsn == 1105,] #one of the samples is apparently from the mother
table(blooddate$fmem) #only one such sample 
blooddate = blooddate[blooddate$fmem != "MOTHER",] #deleting that one sample from mother (since the 5 year sample we are interested  in is from child)


## Merge LOD dataset with levels ##
d2=df
str(d2$hsn) #integer
str(lod$hsn) #character
lod$hsn = as.integer(lod$hsn) #making in to integer
d2 = d2[order(d2$hsn),]
lod = lod[order(lod$hsn),]

d3 = merge(d2, lod, by="hsn", all=T)
doubleID(d3) #no



## Merge 5 year tetanus values in ##
str(d3$hsn) #integer
str(yr5$hsn) #numeric
yr5$hsn = as.integer(yr5$hsn)
d3 = d3[order(d3$hsn),]
yr5 = yr5[order(yr5$hsn),]

d4 = merge(d3, yr5, by="hsn", all=T) #keeping all observations, though yr 5 has only approx 550 observations atm
doubleID(d4)

## Merge blooddate in (SCHEDULED) ##
blooddate = blooddate[, colnames(blooddate) %in% c("hsn", "datecol")] #restricting to only the columns we need
colnames(blooddate) = c("hsn", "datecol_5") #renaming because main dataset already has datecol variable
str(blooddate$hsn) #factor
str(d4$hsn) #integer
blooddate$hsn = as.integer(as.character(blooddate$hsn))
blooddate = blooddate[order(blooddate$hsn),]
d4 = d4[order(d4$hsn),]

d5 = merge(d4, blooddate, by="hsn", all=T)
doubleID(d5) #no

## Merge Visit Date in ##
str(visitdate$hsn) #factor
str(d5$hsn) #integer
visitdate$hsn = as.integer(as.character(visitdate$hsn))
visitdate = visitdate[order(visitdate$hsn),]
d5 = d5[order(d5$hsn),]

d6 = merge(d5, visitdate, by="hsn", all=T)
doubleID(d6) #no

## Setting Final Dataframe ##
d = d6

d = d[d$hsn < 8000,] #exluding the practice IDs etc

#rm(d2, d3, d4, d5,d6, blooddate, visitdate, df, lod, yr5)



## >VARIABLE RECODING/CREATION #####################################################################################


## >>Reformatting Date Variables ------------------------------------------------------------------

## Looking at 5-Year Dates Before Reformatting ##
temp = d[,colnames(d) %in% c("hsn", "datecol_5", "v5date")] #want to have a look and see if my date formatting is inducing differences for some reason
rm(temp)

## Reformatting as Date Variables ##
d$meas1date = as.POSIXct(d$meas1date, format = "%Y-%m-%d") #reformatting as dates so can easily manipulate in next steps
d$meas2date = as.POSIXct(d$meas2date, format = "%Y-%m-%d")
d$tet3date = as.POSIXct(d$tet3date, format = "%Y-%m-%d")
d$tet4date = as.POSIXct(d$tet4date, format = "%Y-%m-%d")
d$datecol_3_5yr = as.POSIXct(d$datecol_3_5yr, format = "%Y-%m-%d")
d$datecol_5 = as.POSIXct(d$datecol_5, format = "%Y-%m-%d")
d$v5date = as.POSIXct(d$v5date, format = "%Y-%m-%d")


## Checking to see if time since lines up (as check of date recode) ##
temp = lubridate::interval(start=d$meas2date, end=d$datecol_3_5yr) %/% months(1) #months between these dates
plot(d$time_meas2, temp) #lines up essentially
rm(temp)

temp = lubridate::interval(start=d$tet4date, end=d$datecol_3_5yr) %/% months(1) #months between these dates
plot(d$time_dtp4, temp) #lines up essentially
rm(temp)


## >>Years Between Visits ------------------------------------------------------------------------------


## First seeing if the two different date variables align ##
temp = d[, colnames(d) %in% c("hsn","datecol_5", "v5date")]
temp[temp$datecol_5 != temp$v5date,] #many differences -- lets see how many days difference are between the dates
temp$difference = difftime( temp$v5date,temp$datecol_5, units="days") 
summary(as.numeric(temp$difference)) #quite a lot of difference
temp$compare = (d$datecol_5 == d$v5date)
table(temp[as.numeric(temp$hsn)<4000,]$compare)

ggplot() +
  geom_freqpoly(data=temp, aes(temp$datecol_5), colour="red") +
  geom_freqpoly(data=temp, aes(temp$v5date), colour="blue") #making some quick figures looking at the date distribution

ggplot() +
  geom_point(data=temp, aes(x=datecol_5, y=v5date))


rm(temp)


summary(temp$datecol_5)
summary(temp$v5date)


## Creating Years Between Visits Variable ##
d$v3.5_to_v5_yrs = lubridate::time_length(difftime(d$datecol_5, d$datecol_3_5yr), "years")
summary(d$v3.5_to_v5_yrs) #makes sense: we would expect 1.5 years between the 5 and 3.5 year visit on average
hist(d$v3.5_to_v5_yrs) #looks right



## >>Time Since Vaccine (5 Year) ---------------------------------------------------------------------------


## Measles ##
time_meas1_5= lubridate::interval(start=d$meas1date, end=d$datecol_5) %/% months(1)
time_meas2_5 = lubridate::interval(start=d$meas2date, end=d$datecol_5) %/% months(1)
summary(time_meas1_5)
summary(time_meas2_5)

time_meas2_5[is.na(time_meas2_5)] = time_meas1_5[is.na(time_meas2_5)] #replace missing values in time2 with values from time1

d$time_meas2_5 = time_meas2_5 

rm(time_meas1_5, time_meas2_5)


## Tetanus/HIB ##
time_dtp3_5 = lubridate::interval(start=d$tet3date, end=d$datecol_5) %/% months(1)
time_dtp4_5 = lubridate::interval(start=d$tet4date, end=d$datecol_5) %/% months(1)
summary(time_dtp3_5)
summary(time_dtp4_5)

time_dtp4_5[is.na(time_dtp4_5)] = time_dtp3_5[is.na(time_dtp4_5)] #replace missing values in time4 with time3

d$time_dtp4_5 = time_dtp4_5

rm(time_dtp3_5, time_dtp4_5)




## >>Time since last vaccine Categories (for T1) ---------------------------------------------------------------

hist(d$time_meas2)
quantile(d$time_meas2, probs=c(0.33,0.66), na.rm=T)
hist(d$time_dtp4)
quantile(d$time_dtp4, probs=c(0.33,0.66), na.rm=T)

# Creating categorical variables based on distributions (MOST scores are between 23-25 months)
d$time_meas2_cat[d$time_meas2<23 & !is.na(d$time_meas2)] = "<23"
d$time_meas2_cat[d$time_meas2>=23 & d$time_meas2<25 & !is.na(d$time_meas2)] = "23-25"
d$time_meas2_cat[d$time_meas2>25 & !is.na(d$time_meas2)] = ">25"
d$time_meas2_cat = factor(d$time_meas2_cat, levels=c("<23", "23-25", ">25"))
psych::describeBy(d$time_meas2, group=d$time_meas2_cat)

d$time_dtp4_cat[d$time_dtp4<23 & !is.na(d$time_dtp4)] = "<23"
d$time_dtp4_cat[d$time_dtp4>=23 & d$time_dtp4<25 & !is.na(d$time_dtp4)] = "23-25"
d$time_dtp4_cat[d$time_dtp4>25 & !is.na(d$time_dtp4)] = ">25"
d$time_dtp4_cat = factor(d$time_dtp4_cat, levels=c("<23", "23-25", ">25"))
psych::describeBy(d$time_dtp4, group=d$time_dtp4_cat)



## >>Not Protected Using 5 Year Values --------------------------------------------------------------------


## Measles ##
d$notprot_meas5yr = NA 
d$notprot_meas5yr[d$measles_5yr <=0.25 & !is.na(d$measles_5yr)] = 1 #not protected
d$notprot_meas5yr[d$measles_5yr > 0.25 & !is.na(d$measles_5yr)] = 0 #protected
psych::describeBy(d$measles_5yr, d$notprot_meas5yr)


## Tetanus ##
d$notprot_tet5yr = NA
d$notprot_tet5yr[d$tetanus_5yr <= 0.1 & !is.na(d$tetanus_5yr)] = 1 
d$notprot_tet5yr[d$tetanus_5yr > 0.1 & !is.na(d$tetanus_5yr)] = 0
psych::describeBy(d$tetanus_5yr, d$notprot_tet5yr)


## HIB ##
d$notprot_hib5yr = NA
d$notprot_hib5yr[d$hib_5yr <= 1 & !is.na(d$hib_5yr)] = 1
d$notprot_hib5yr[d$hib_5yr > 1 & !is.na(d$hib_5yr)] = 0
psych::describeBy(d$hib_5yr, d$notprot_hib5yr)


## >>Log of 5-Year Antibody Levels --------------------------------------------------------


d$log_measles_5yr = log10(d$measles_5yr)
d$log_tetanus_5yr = log10(d$tetanus_5yr)
d$log_hib_5yr = log10(d$hib_5yr)



## >>Change in Antibody Levels --------------------------------------------------------------

d$measles_change =  d$measles_5yr - d$mev
d$tetanus_change = d$tetanus_5yr- d$tet
d$hib_change =  d$hib_5yr - d$hib

summary(d$measles_change)
summary(d$tetanus_change)
summary(d$hib_change)


## >>Change per Year in Antibody Levels -----------------------------------------------------

d$measles_change_peryr = d$measles_change/d$v3.5_to_v5_yrs
d$tetanus_change_peryr = d$tetanus_change/d$v3.5_to_v5_yrs
d$hib_change_peryr = d$hib_change/d$v3.5_to_v5_yrs

summary(d$measles_change_peryr)
summary(d$tetanus_change_peryr)
summary(d$hib_change_peryr)

## >>Not Protected Including those Diagnosed ------------------------------------------------

table(d$notprot_meas, d$bqbmeasles, useNA="always") #looking at crosstabs of those protected and those who were diagnosed
table(d$notprot_tet, d$bqbtetanus, useNA="always")

# Creating a variable with value of 1 if the antibodies are too low OR child diagnosed with measles 
d$meas_low_or_diag[(d$notprot_meas==1 | d$bqbmeasles==1) ] = 1
d$meas_low_or_diag[d$notprot_meas==0 & d$bqbmeasles==0  & !is.na(d$notprot_meas) & !is.na(d$bqbmeasles)] = 0
addmargins(table(d$bqbmeasles, d$meas_low_or_diag, useNA="always"))
addmargins(table(d$notprot_meas, d$meas_low_or_diag, useNA="always")) #NOTE: if missing bqbmeasles, but has low antibodies, catgorized as 1 still (and vice versa)

# Creating a variable with value of 1 if tetanus antibodies too low OR child diagnosed with tetanus
d$tet_low_or_diag[(d$notprot_tet==1 | d$bqbtetanus==1)] = 1
d$tet_low_or_diag[d$notprot_tet==0 & d$bqbtetanus==0 & !is.na(d$notprot_tet) & !is.na(d$bqbtetanus)] = 0
addmargins(table(d$bqbtetanus, d$tet_low_or_diag, useNA="always"))
addmargins(table(d$notprot_tet, d$tet_low_or_diag, useNA="always"))


## >>Breastfeeding -----------------------------------------------------------------------

summary(d$bf_mo_nc) #this is the updated breastfeeding variable, not censored 
quantile(d$bf_mo_nc, probs=c(0.33,0.67), na.rm=T)
histogram(d$bf_mo_nc)

d$bf_cat_nc = NA
d$bf_cat_nc[d$bf_mo_nc<=12 & !is.na(d$bf_mo_nc)] = "<=12"
d$bf_cat_nc[d$bf_mo_nc>12 & d$bf_mo_nc<=24 & !is.na(d$bf_mo_nc)] = "12-24"
d$bf_cat_nc[d$bf_mo_nc>24 & !is.na(d$bf_mo_nc)] = ">24"

d$bf_cat_nc = factor(d$bf_cat_nc, 
                     levels=c("<=12", "12-24", ">24"),
                     labels=c("<=12", "12-24", ">24"))

psych::describeBy(x=d$bf_mo_nc, group=d$bf_cat_nc)




## >VARIABLE SETS ###############################################################################

t1_vars=c("mage_cat", "edu_cat", "poverty_food_i", "insecure_low", "parity2",
          "hiv_combined_i", "smoke_ever", "alcohol", "mal_iom_high_act",
          "csex", "bwt_cat", "ga_cat", "deliv_method_i", "bf_cat_nc",
          "time_meas2_cat", "time_dtp4_cat")

exp_list_original = c("ln_opdde_i", "ln_ppdde_i", "ln_opddt_i", "ln_ppddt_i",
             "log_cisdbca_sg", "log_cisdcca_sg", "log_transdcca_sg", "log_pba_sg")

exp_list = c( "ln_ppddt_i","ln_ppdde_i",
                      "log_cisdbca_sg", "log_cisdcca_sg", "log_transdcca_sg", "log_pba_sg")

exp_list_notlog = c("n_opdde_i", "n_ppdde_i", "n_opddt_i", "n_ppddt_i", "n_pcb_sum_i",
                    "cisdbca_sg", "cisdcca_sg", "transdcca_sg", "pba_sg")

lodloq_vars = colnames(d)[grep("lod|loq", colnames(d))]

adjset_meas = c("mage", 
                "income_pers", "parity", "hiv_combined_i", "csex",
                "rms::rcs(time_meas2,3)")

adjset_meas_lintime = c("mage", 
                        "income_pers", "parity", "hiv_combined_i", "csex",
                        "time_meas2")

adjset_tet = c("mage", 
               "income_pers", "parity", "hiv_combined_i", "csex",
               "rms::rcs(time_dtp4,3)") 

adjset_tet_lintime = c("mage", 
                       "income_pers", "parity", "hiv_combined_i", "csex",
                       "time_dtp4") 



## >FLOW CHART NUMBERS ##############################################################################


## >>Going from 5-Year -----------------------------------------------------------------------------------

nrow(d[!is.na(d$v5date) | (!is.na(d$measles_5yr) | !is.na(d$tetanus_5yr)),]) #N=641 who supposedly were at the 5 year visit -- perhaps we gained values? 
    #yes, gaining this person makes sense--Thomas added more samples      

d$in5visit = 0 #creating a flag for those who were at the 5-year visit
d$in5visit[!is.na(d$v5date) | (!is.na(d$measles_5yr) | !is.na(d$tetanus_5yr))] = 1
table(d$in5visit)

f1 = d[d$in5visit==1,] #restrict to those who participated in year 5 visit

addmargins(table(f1$measles2, f1$dtp4, useNA="always")) #looking at the breakdown of vaccine schedule
    #5 with no measles vaccine schedule
    #41 with no dtp schedule (31 + 10)
    #44 with neither vaccine schedule (42 + 2)

f2 = f1[(f1$measles2==1 & !is.na(f1$measles2)) | ((f1$dtp4==1) & !is.na(f1$dtp4)),] #N=595 restricting to those with full vaccine schedule of EITHER measles AND/OR DTaP

addmargins(table(f1$bqbmeasles, f1$bqbtetanus))
    #28 diagnosed with only measles
    #1 diagnosed with both measles and tetanus 

f3 = f2[f2$bqbmeasles !=1 | f2$bqbtetanus != 1, ] #N=594 restricting to those who are free of EITHER measles diagnosis and/or tetanus diagnosis
    
nrow(f1[(is.na(f1$measles_5yr) & is.na(f1$tetanus_5yr)) ,]) #89 missing blood measures
     #89 missing 5 year blood measures

nrow(f1[(is.na(f1$measles_5yr) & is.na(f1$tetanus_5yr)) & is.na(f1$datecol_5),]) #17 are missing collection date AND the blood measures
     #assume that 17/89 missing legitimately and 72/89 are still in South Africa


f4 = f3[!is.na(f3$measles_5yr) | !is.na(f3$tetanus_5yr),] #N=565 restricting to those with at least one of the measures (actually the same N whether or not | or & is used)
        #number before updated file: N=513

nrow(f1[!is.na(f1$tetanus_5yr) & f1$bqbtetanus !=1 & f1$dtp4==1 & !is.na(f1$dtp4),]) #N=526 for tetanus models
nrow(f1[!is.na(f1$measles_5yr) & f1$bqbmeasles !=1 & f1$measles2==1 & !is.na(f1$measles2),]) #N=538 for measles models


flag_5yr = f4$hsn #creating a flag to use later

rm(f1,f2,f3,f4) #remove temporary dataframes 



## >>Going from 3.5 Year -----------------------------------------------------------------------------------------------


nrow(d[!is.na(d$datecol_3_5yr) & (!is.na(d$mev) | !is.na(d$tet)),]) #N=642 who supposedly were at the 3.5 year visit
      #this is consistent with what I did before

d$in3visit = 0 #creating a flag for those who were at the 5-year visit
d$in3visit[!is.na(d$datecol_3_5yr) & (!is.na(d$mev) & !is.na(d$tet))] = 1
table(d$in3visit)

f1 = d[d$in3visit==1 & !is.na(d$in3visit),] #restrict to those who participated in year 3 visit
addmargins(table(d$in5visit, d$in3visit, useNA="always")) #looking at overlap between timepoints 


addmargins(table(f1$measles2, f1$dtp4, useNA="always")) #looking at the breakdown of vaccine schedule
  #5 with no measles vaccine schedule
  #43 with no dtp schedule (33+10)
  #46 with neither vaccine schedule (43+2+1)
table(d$dtp4[d$in5visit==0 & d$in3visit==1]) #there are four more children without vaccine in the year 3.5 visit that are not in the 5 year visit
table(d$dtp4[d$in5visit==1 & d$in3visit==0]) #also one child without dtp vaccine in the 5 year that is not in the 3.5 year


f2 = f1[(f1$measles2==1 & !is.na(f1$measles2)) | ((f1$dtp4==1) & !is.na(f1$dtp4)),] #N=596 restricting to those with full vaccine schedule of EITHER measles AND/OR DTaP


addmargins(table(f1$bqbmeasles, f1$bqbtetanus))
    #24 diagnosed with only measles
    #1 diagnosed with both measles and tetanus 
table(d$bqbmeasles[d$in5visit ==1 & d$in3visit==0]) #there are four diagnosed with measles who participated in visit 5 but did not participate in visit 3

f3 = f2[f2$bqbmeasles !=1 | f2$bqbtetanus != 1, ] #N=595 restricting to those who are free of EITHER measles diagnosis and/or tetanus diagnosis
    #the original number was 596? Why the reduction by 1?

nrow(f1[(is.na(f1$mev) & is.na(f1$tet)) ,]) #0 missing blood measures --> the difference between N=667 and N=642 is the 25 missing blood measures
    #94 missing 3 year blood measures


nrow(f1[!is.na(f1$tet) & f1$bqbtetanus !=1 & f1$dtp4==1 & !is.na(f1$dtp4),]) #N=552 for tetanus models
nrow(f1[!is.na(f1$mev) & f1$bqbmeasles !=1 & f1$measles2==1 & !is.na(f1$measles2),]) #N=569 for measles models


flag_3yr = f3$hsn #creating a flag to use later


not_in_3yr = flag_5yr[!(flag_5yr %in% flag_3yr)]


## >>Combining to get overall N & RESTRICTING DF-------------------------------------------------------------------------

d$include_for_gee = 0 #making a flag that is for inclusion based on available 3.5 and 5 year data ==> Note that we only need one measure at either timepoint to include
d$include_for_gee[d$hsn %in% flag_3yr] = 1
table(d$include_for_gee) #looks right
d$include_for_gee[d$hsn %in% flag_5yr] = 1
table(d$include_for_gee) #overall there are N=611 IDs that can be included 
    #originally had only 607 that could be included, but we have gained data

d_gee = d[d$include_for_gee ==1,] #restrict to only those meant to include in GEE 

nrow(d_gee[(!is.na(d$mev) | !is.na(d$measles_5yr)) &
             d$bqbmeasles !=1 &
             d$measles2==1 & !is.na(d$measles2),]) #N=582 for measles models
                #previously N=578 before additional samples

nrow(d_gee[(!is.na(d$tet) | !is.na(d$tetanus_5yr)) &
             d$bqbtetanus !=1 &
             d$dtp4==1 & !is.na(d$dtp4),]) #N=568 for measles models
                #previously N=564 before additional samples

table(d$in3visit, d$in5visit, useNA="always") #overall counts of who was in the 3.5 and 5 year visits 


## Getting a sense of who is contributing to the models ##
d$in3yrmeas = 0 #creating a series of flags to figure out how the model memberships overlap
d$in3yrmeas[d$measles2==1 & !is.na(d$measles2) & !is.na(d$mev) & d$bqbmeasles !=1] = 1
table(d$in3yrmeas) #N=569 in the measles models for 3.5 year (lines up with flowchart work above)

d$in3yrdtp = 0 
d$in3yrdtp[d$dtp4==1 & !is.na(d$dtp4) & (!is.na(d$tet | d$hib)) & d$bqbtetanus !=1] = 1 
table(d$in3yrdtp) #N=552, lines up with flowchart above

d$in5yrmeas = 0
d$in5yrmeas[d$measles2==1 & !is.na(d$measles2) & !is.na(d$measles_5yr) & d$bqbmeasles !=1]  = 1
table(d$in5yrmeas) #N=538, lines up with flowchart

d$in5yrdtp = 0
d$in5yrdtp[d$dtp4==1 & !is.na(d$dtp4) & (!is.na(d$tetanus_5yr | d$hib_5yr)) & d$bqbtetanus !=1] = 1 
table(d$in5yrdtp) #N=526, lines up with flowchart work above

d$inmodel = ifelse((d$in3yrmeas + d$in3yrdtp + d$in5yrmeas + d$in5yrdtp) > 0, 1,0)
table(d$inmodel) #now getting N=610 as opposed to N=611...

unique(d[,colnames(d) %in% c("in3yrmeas", "in3yrdtp", "in5yrmeas", "in5yrdtp")]) #see each of the unique patterns of models that make up the sample 

nrow(d[d$in3yrmeas==1 &
       d$in3yrdtp==1  &
       d$in5yrmeas==1 &
       d$in5yrdtp==1,])#now getting counts for each ofthe unique combinations
          #N=487

nrow(d[d$in3yrmeas==1 &
         d$in3yrdtp==0  &
         d$in5yrmeas==1 &
         d$in5yrdtp==0,])#N=38

nrow(d[d$in3yrmeas==1 &
         d$in3yrdtp==1  &
         d$in5yrmeas==0 &
         d$in5yrdtp==0,])#N=40

nrow(d[d$in3yrmeas==0 &
         d$in3yrdtp==0  &
         d$in5yrmeas==1 &
         d$in5yrdtp==1,])#N=13

nrow(d[d$in3yrmeas==1 &
         d$in3yrdtp==0  &
         d$in5yrmeas==0 &
         d$in5yrdtp==0,])#N=4

nrow(d[d$in3yrmeas==0 &
         d$in3yrdtp==1  &
         d$in5yrmeas==0 &
         d$in5yrdtp==0,])#N=2

nrow(d[d$in3yrmeas==0 &
         d$in3yrdtp==1  &
         d$in5yrmeas==0 &
         d$in5yrdtp==1,])#N=23

nrow(d[d$in3yrmeas==0 &
         d$in3yrdtp==0  &
         d$in5yrmeas==0 &
         d$in5yrdtp==1,])#N=3



## >REARRANGE DATAFRAME TO LONG FORMAT ###############################################################


## Segmenting Data ##
d_covs = d_gee[, colnames(d_gee) %in% c('hsn', adjset_tet_lintime, adjset_meas_lintime, 
                                "dtp4", "measles2", "bqbmeasles", "bqbtetanus",
                                "poverty_food_i", "insecure_low",
                                exp_list, t1_vars, exp_list_notlog, lodloq_vars, 
                                "in3visit", "in5visit")] #creating a dataset of the covariates to merge in with the long dataset

d3yr = d_gee[, c('hsn', 'logmev', 'logtet', 'loghib', 
                              'notprot_meas', 'notprot_tet', 'notprot_hib')]

d5yr = d_gee[, c('hsn', "log_measles_5yr", "log_tetanus_5yr", "log_hib_5yr", 
                            "notprot_meas5yr", "notprot_tet5yr", "notprot_hib5yr")]


## Adding Variables for Time ##
d3yr$time = 0
d5yr$time = 1

## Stacking Time Points ##
colnames(d3yr)
colnames(d5yr)
colnames(d5yr) = colnames(d3yr) #making colnames the same 
dlong1 = rbind(d3yr, d5yr) #stacking the timepoints


## Merging in Covariates ##
str(dlong1$hsn) #integer
str(d_covs$hsn) #integer
dlong1 = dlong1[order(dlong1$hsn),] #sorting
d_covs = d_covs[order(d_covs$hsn), ]
dlong = merge(dlong1, d_covs, by="hsn") #covariate values present for both observations from each ID 



## >RESTRICT DATAFRAME #############################################################################

# RESTRICT DATAFRAME TO THOSE WHO HAVE NOT COMPLETED VACCINE SCHEDULE (SEE NOTE FROM THOMAS)

# addmargins(table(dlong$measles2[dlong$time==0], dlong$dtp4[dlong$time==0], useNA="always"))
# #of 752, 585 have full vaccine schedule
# #2 had neither vaccine (dtp nor measles)
# #10 had the measles vaccine, but not the dtp
# #34 had the measles vaccine, but are missing information ('DIDN'T KNOW' I think) on the dtp vaccine
# #6 had the dtp vaccine, but are missing information (DIDNT KNOW - 9) on the meaesles vaccine
# #69 do not know about either vaccine
# #46 are missing data (questionnaire data missing)



dr=dlong[(dlong$measles2==1 & !is.na(dlong$measles2) & dlong$bqbmeasles !=1) |
           (dlong$dtp4==1 & !is.na(dlong$dtp4) & dlong$bqbtetanus !=1),]
  
dr_meas = dr[dr$measles2==1 & !is.na(dr$measles2),] #creating a restricted dataframe specific to measles. First keeping only those with full measles vaccine schedule
dr_meas = dr_meas[dr_meas$bqbmeasles != 1, ] #dropping anyone diagnosed with measles. Down to 582

dr_tet = dr[dr$dtp4==1 & !is.na(dr$dtp4),] #creating restricted dataframe specific to tetanus. First keeping only those w full DTaP vacc shd
dr_tet = dr_tet[dr_tet$bqbtetanus != 1, ] #dropping anyone diagnosed with tetatnus. Only lose one additional. Final N = 568


table(dr$in3visit, dr$in5visit)

summary(dlong$logmev)
summary(dlong$logtet)
summary(dlong$time_dtp4)
summary(dlong$time_meas2)

check_ids = dr[dr$hsn %in% not_in_3yr,]


## Removing extra dataframes ##
# rm(blooddate, d, d_covs, d_gee, d2, d3, d3yr, d4, d5, d5yr, d6, df, dlong1,
#    dlong, f1, f2, f3, lod, visitdate, yr5) #tidies up workspace 



## >FUNCTIONS ########################################################################


## >>Descriptive Functions -----------------------------------------------------------------------

## For Getting Mode of Factors ##
Mode = function(x){ 
  ta = table(x)
  tam = max(ta)
  if (all(ta == tam))
    mod = NA
  else
    if(is.numeric(x))
      mod = as.numeric(names(ta)[ta == tam])
  else
    mod = names(ta)[ta == tam]
  return(mod)
}


## Table 1: Categorical ##
catTable1 = function(df, varlist) { 
  
  t1 = list()
  for (i in seq_along(varlist)) {
    
    levels = as.data.frame(table(df[,colnames(df)==varlist[i]])) #Getting frequencies
    sum=sum(levels$Freq)
    levels$Pct = round(levels$Freq/sum*100,1) #Getting percentages
    
    levels$Var1 = as.character(levels$Var1) #Converting all to characters to make merge easier
    levels$Freq = as.character(levels$Freq)
    levels$Pct = as.character(levels$Pct)
    
    name = as.data.frame(cbind("Var1" = varlist[i], "Freq" = NA, "Pct" = NA)) #getting name of the variable to put on top
    name$Var1 = as.character(name$Var1)
    name$Freq = as.character(name$Freq)
    name$Pct = as.character(name$Pct)
    
    t1[i] = list(
      rbind(name,levels)
    )
  }
  
  out = as.data.frame(do.call("rbind", t1))
  
  return(out)
  
}


## Table 2 ##
makeT2 = function(df, varlist) {
  
  t2 = list()
  for(i in seq_along(varlist)) {
    
    gm = EnvStats::geoMean(df[, colnames(df)==varlist[i]], na.rm=T)
    gsd = EnvStats::geoSD(df[, colnames(df)==varlist[i]], na.rm=T)
    gm_sd = paste(round(gm,2), " (", round(gsd,2), ")", sep="")
    
    q = rbind(quantile(df[, colnames(df)==varlist[i]], probs=c(0,0.05, 0.25, 0.50, 0.75, 0.95,1), na.rm=T))
    
    nobs = sum(!is.na(df[,colnames(df)==varlist[i]]))
    
    compound = varlist[i]
    
    t2[i] = list(
      cbind(compound, nobs, gm_sd, round(q,2))
    )
  }
  
  out = as.data.frame(do.call("rbind", t2))
  return(out)
  
}


## >>Effects at Single Timepoints ----------------------------------------------------------------

## to get point RR estimates ##
getEstimatesME = function(df, outcome, exposures, adjustment) {
  
  out=list()
  for(i in seq_along(exposures)) {
    m2 = glm(data=df,
             as.formula(paste(outcome, "~", exposures[i], "+", paste(adjustment, collapse="+"))),
             family=poisson(link="log")) #poisson reg model using an offset 
    
    m2_cov = sandwich::vcovHC(m2, type="HC0") #Robust covariance matrix estimators a la White for panel models
    m2_se = sqrt(diag(m2_cov))
    
    cs = coef(m2) #getting the coefficients
    estimate = unname(cs[grep(exposures[i], names(cs))])
    
    robust_se = unname(m2_se[grep(exposures[i], names(m2_se))]) #getting the robust standard error
    
    rr = exp(estimate) #getting the risk ratio estimates and confidence intervals using the robust SE
    lo_rr = exp(estimate -1.96*robust_se)
    hi_rr = exp(estimate + 1.96*robust_se)
    pval = 2*pnorm(abs(estimate/robust_se), lower.tail=F)
    
    nobs = nobs(m2)
    predictor = exposures[i]
    dv = outcome
    
    out[i] = list(
      as.data.frame(cbind(
        predictor, dv, nobs,
        rr, lo_rr, hi_rr, pval, robust_se
      ))
    )
    rm(predictor, dv, nobs, rr, lo_rr, hi_rr, pval, robust_se )
    
  }
  
  output = do.call("rbind", out)
  
  output$rr = as.numeric(as.character(output$rr)) #making sure all numeric variables are coded as such
  output$lo_rr = as.numeric(as.character(output$lo_rr))
  output$hi_rr = as.numeric(as.character(output$hi_rr))
  output$pval = as.numeric(as.character(output$pval))
  output$robust_se = as.numeric(as.character(output$robust_se))
  
  
  return(output)
  
}


## to clean point estimates for putting in table ##
cleanEstimatesME = function(input) {
  
  input$rr_ci = paste(round(input$rr,2), " (", round(input$lo_rr,2), ", ", round(input$hi_rr, 2), ")", sep="")
  
  output = input[,c("predictor", "nobs", "rr_ci")]
  return(output)
  
}


## to get point beta (se) estimates ##
getLinEstimatesME = function(df, outcome, exposures, adjustment) {
  
  out=list()
  for(i in seq_along(exposures)) {
    
    m2 = glm(data=df,
             as.formula(paste(outcome, "~", exposures[i], "+", paste(adjustment, collapse="+"))),
             family=gaussian) #note that when there are nonlinear terms in the model (eg time), the coefficients are not directly interpretable themselves 
    
    
    s = summary(m2)
    coef = as.data.frame(s$coefficients)
    beta = coef$Estimate[rownames(coef)==exposures[i]] #exponentiated so that this should now be the increase in mev per 10-fold increase in the exposure
    se = coef$`Std. Error`[rownames(coef)==exposures[i]]
    
    beta_ci = paste(round(beta,2), " (", round((beta-(1.96*se)),2), ", ", round((beta+(1.96+se)),2), ")", sep="")
    
    nobs_lin = nobs(m2)
    predictor_lin = exposures[i]
    dv_lin = outcome
    
    out[i] = list(
      as.data.frame(cbind(
        predictor_lin, dv_lin, nobs_lin,
        beta_ci
      ))
    )
    rm(predictor_lin, dv_lin, nobs_lin, beta_ci)
    
  }
  
  output = do.call("rbind", out)
  
  return(output)
  
}


## to test interaction ## 
testIntx = function(df, outcome, exposures, adjustment, modifier) {
  
  out=list()
  for(i in seq_along(exposures)) {
    
    f = as.formula(paste(outcome, "~", exposures[i], "+", exposures[i], "*", modifier, "+", paste(adjustment, collapse="+")))
    
    m2 = glm(data=df, f, family = poisson(link="log")) #poisson regression model 
    
    m2_cov = sandwich::vcovHC(m2, type="HC0") #Robust covariance matrix estimators a la White for panel models
    m2_se = sqrt(diag(m2_cov))
    
    cs = coef(m2) #getting the coefficients
    estimate = cs[grep(exposures[i], names(cs))] 
    estimate = unname(estimate[grep(modifier, names(estimate))]) #gets just the coefficient for the interaction term
    
    robust_se = m2_se[grep(exposures[i], names(m2_se))]
    robust_se = unname(robust_se[grep(modifier, names(robust_se))]) #getting the robust standard error for the interaction term
    
    pval = 2*pnorm(abs(estimate/robust_se), lower.tail=F)
    
    nobs = nobs(m2)
    predictor = exposures[i]
    dv = outcome
    
    out[i] = list(
      as.data.frame(cbind(
        predictor, dv, nobs, estimate, robust_se, pval 
      ))
    )
    rm(predictor, dv, nobs, pval, robust_se )
    
  }
  
  output = do.call("rbind", out)
  
  output$pval = as.numeric(as.character(output$pval))
  output$robust_se = as.numeric(as.character(output$robust_se))
  output$estimate = as.numeric(as.character(output$estimate))
  output$siglabel = ""
  output$siglabel[output$pval<.05] = "*"
  
  return(output)
  
}


## to get stratified point RR estimates ##
getStratEstimatesME = function(df, outcome, exposures, adjustment, l1, l2, stratby) {
  
  s = df[, colnames(df) == stratby]
  
  df1 = df[s==l1,]
  df2 = df[s==l2,]
  
  out=list()
  
  for(i in seq_along(exposures)) {
    
    ## For first level of modifier ##
    m1 = glm(data=df1,
             as.formula(paste(outcome, "~", exposures[i], "+", paste(adjustment, collapse="+"))),
             family=poisson(link="log")) #poisson reg model using an offset 
    
    m1_cov = sandwich::vcovHC(m1, type="HC0") #Robust covariance matrix estimators a la White for panel models
    m1_se = sqrt(diag(m1_cov))
    
    cs = coef(m1) #getting the coefficients
    estimate = unname(cs[grep(exposures[i], names(cs))])
    
    robust_se = unname(m1_se[grep(exposures[i], names(m1_se))]) #getting the robust standard error
    
    rr = exp(estimate) #getting the risk ratio estimates and confidence intervals using the robust SE
    lo_rr = exp(estimate -1.96*robust_se)
    hi_rr = exp(estimate + 1.96*robust_se)
    
    predictor = exposures[i]
    dv = outcome
    level = l1
    
    out1 = as.data.frame(cbind(predictor, dv, level, rr, lo_rr, hi_rr))
    
    rm(predictor, dv, level, rr, lo_rr, hi_rr)
    
    ## For second level of modifier ##
    m2 = glm(data=df2,
             as.formula(paste(outcome, "~", exposures[i], "+", paste(adjustment, collapse="+"))),
             family=poisson(link="log")) #poisson reg model using an offset 
    
    m2_cov = sandwich::vcovHC(m2, type="HC0") #Robust covariance matrix estimators a la White for panel models
    m2_se = sqrt(diag(m2_cov))
    
    cs = coef(m2) #getting the coefficients
    estimate = unname(cs[grep(exposures[i], names(cs))])
    
    robust_se = unname(m2_se[grep(exposures[i], names(m2_se))]) #getting the robust standard error
    
    rr = exp(estimate) #getting the risk ratio estimates and confidence intervals using the robust SE
    lo_rr = exp(estimate -1.96*robust_se)
    hi_rr = exp(estimate + 1.96*robust_se)
    
    predictor = exposures[i]
    dv = outcome
    level = l2
    
    out2 = as.data.frame(cbind(predictor, dv, level, rr, lo_rr, hi_rr))
    
    
    
    ## Binding together both levels ##
    out[i] = list(
      as.data.frame(rbind(out1, out2))
    )
    
    rm(predictor, dv,  rr, lo_rr, hi_rr, robust_se )
    
  }
  
  output = do.call("rbind", out)
  
  output$rr = as.numeric(as.character(output$rr)) #making sure all numeric variables are coded as such
  output$lo_rr = as.numeric(as.character(output$lo_rr))
  output$hi_rr = as.numeric(as.character(output$hi_rr))
  
  return(output)
  
}


## Test Interaction on Continuous Outcomes ##
runMEIntxCont = function(outcome, exposures, adjustment, modifier, df) {
  
  out=list()
  for (i in seq_along(exposures)) {
    
    f_intx =  as.formula(paste(outcome, "~", exposures[i],"+", modifier, "+", modifier, "*", exposures[i], "+", 
                               paste(adjustment, collapse="+")))
    
    m_intx = glm(data=df,
                 f_intx,
                 family=gaussian)
    
    coef_intx = as.data.frame(coefficients(summary(m_intx))) #extracting relevent coefficients 
    coef_intx = coef_intx[grep(exposures[i], rownames(coef_intx)),] #restricting to only estimates pertaining to the exposure
    coef_intx = coef_intx[grep(modifier, rownames(coef_intx)), ] #now restricting to only the estimate of the interaction between exposure and modifier
    
    est_robust95ci = paste(round(exp(coef_intx$Estimate),3), " [", 
                           round(exp(coef_intx$Estimate-(1.96*coef_intx$`Std. Error`)),3),
                           ", ",
                           round(exp(coef_intx$Estimate+(1.96*coef_intx$`Std. Error`)),3),
                           "]", sep="")
    
    pval = coef_intx$`Pr(>|t|)`
    
    nobs = nobs(m_intx)
    predictor = exposures[i]
    dv = outcome
    
    out[i] = list(
      as.data.frame(cbind(
        predictor, dv, nobs, est_robust95ci, pval 
      ))
    )
    rm(predictor, dv, nobs, pval, est_robust95ci)
    
    
    
  }
  
  output = do.call("rbind", out)
  
  output$pval = as.numeric(as.character(output$pval))
  output$siglabel = ""
  output$siglabel[output$pval<.10] = "*"
  
  return(output) 
  
}



## Getting Stratified Estimates for Continuous Outcome @ Single Timpepoint ##
getMEStratEstsCont = function(outcome, exposures, adjustment, l1, l2, stratby, df) {
  
  s = df[, colnames(df) == stratby]
  
  df1 = df[s==l1 & !is.na(s),]
  df2 = df[s==l2 & !is.na(s),]
  
  
  out=list()
  for (i in seq_along(exposures)) {
    
 
    ## For first level of stratifying variable ##
    m_intx = glm(data=df1,
             as.formula(paste(outcome, "~", exposures[i], "+", paste(adjustment, collapse="+"))),
             family=gaussian) #note that when there are nonlinear terms in the model (eg time), the coefficients are not directly interpretable themselves 
    
    coef_intx = as.data.frame(coefficients(summary(m_intx))) #extracting relevent coefficients 
    coef_intx = coef_intx[rownames(coef_intx) ==exposures[i],]
    
    rr = 100*(exp(coef_intx$Estimate)-1) #getting the percent (?) beta and robust 95% CI
    lo_rr = 100*(exp(coef_intx$Estimate - (1.96*coef_intx$`Std. Error`))-1)
    hi_rr = 100*(exp(coef_intx$Estimate + (1.96*coef_intx$`Std. Error`))-1)
    
    predictor = exposures[i]
    dv = outcome
    level = l1
    
    out1 = as.data.frame(cbind(predictor, dv, level, rr, lo_rr, hi_rr))
    
    rm(predictor, dv, level, rr, lo_rr, hi_rr, coef_intx, m_intx)
    
    
    ## For second level of stratifying variable ##
    m_intx = glm(data=df2,
                 as.formula(paste(outcome, "~", exposures[i], "+", paste(adjustment, collapse="+"))),
                 family=gaussian) 
    
    coef_intx = as.data.frame(coefficients(summary(m_intx))) #extracting relevent coefficients 
    coef_intx = coef_intx[rownames(coef_intx) ==exposures[i],]
    
    rr = 100*(exp(coef_intx$Estimate)-1) #getting the beta and robust 95% CI
    lo_rr = 100*(exp(coef_intx$Estimate - (1.96*coef_intx$`Std. Error`))-1)
    hi_rr = 100*(exp(coef_intx$Estimate + (1.96*coef_intx$`Std. Error`))-1)
    
    predictor = exposures[i]
    dv = outcome
    level = l2
    
    out2 = as.data.frame(cbind(predictor, dv, level, rr, lo_rr, hi_rr))
    
    rm(predictor, dv, level, rr, lo_rr, hi_rr)
    
    
    
    ## Binding together both levels ##
    out[i] = list(
      as.data.frame(rbind(out1, out2))
    )
    
    
    
    
  }
  output = do.call("rbind", out)
  
  output$rr = as.numeric(as.character(output$rr)) #making sure all numeric variables are coded as such
  output$lo_rr = as.numeric(as.character(output$lo_rr))
  output$hi_rr = as.numeric(as.character(output$hi_rr))
  
  colnames(output) = c("predictor", "dv", "level", "beta", "lo_ci", "hi_ci")
  
  return(output)
}


## To Plot Stratified Estimates ##
plotStratEstimates = function(df,sigdf, savename) {
  plot = 
    ggplot() + 
    geom_pointrange(data=df,
                    aes(x=predictor, y=rr, ymin=lo_rr, ymax=hi_rr, fill=level), pch=21,
                    position=position_dodge(width=0.6), size=2) +
    geom_text(data=sigdf, aes(x=predictor, y=2, label=siglabel), size=24) + 
    facet_grid(dv~.) + 
    scale_y_continuous(trans="log", breaks=c(0.5, 1, 2, 4)) +
    #coord_cartesian(ylim=c(0.25,6)) +
    scale_fill_manual(values=c("black", "white")) + 
    scale_x_discrete(labels = c(#expression(paste(italic("o,p'"), "-DDE")),
                                expression(paste(italic("p,p'"), "-DDT")),
                                expression(paste(italic("p,p'"), "-DDE")),
                                #expression(paste(italic("o,p'"), "-DDT")),
                                expression(paste(italic("cis"), "-DBCA")),
                                expression(paste(italic("cis"), "-DCCA")),
                                expression(paste(italic("trans"), "-DCCA")),
                                "3-PBA")) +
    geom_hline(yintercept=1, linetype='dashed') +
    labs(
      y="Relative Risk per 10-fold increase\n",
      x="\nExposure") +
    theme_light() +
    theme(
      text = element_text(family = "Century Gothic"),
      legend.title = element_blank(),
      legend.text = element_text(size=20),
      legend.key.height = unit(3, "cm"),
      legend.key.width = unit(3, "cm"),
      legend.position = "top",
      axis.text = element_text(size=20),
      axis.text.x = element_text(angle = 45, size=20, hjust=1),
      axis.title = element_text(size=25),
      plot.caption = element_text(size=20, hjust=0),
      strip.background = element_blank(),
      strip.text = element_text(size=20, colour="black")
    ) 
  ggsave(plot, file=paste(savename,".png", sep=""), height=15, width =12)
  
}


## >>GEE Effects -------------------------------------------------------------------------------


## Creating a Prediction Dataframe for Tetanus ##
createPredictDf = function(expvect, expname, mage=26.3, income_pers=522, parity=1.09,
                           hiv_combined_i="HIV Negative", csex="Boy", time_dtp4=23, time=0) { #defaults are the mean/mode in the dr_tet dataset
  
  minexp = min(expvect, na.rm=T)
  maxexp = max(expvect, na.rm=T )
  increment = (maxexp - minexp)/1000 
  
  range_exp = seq(from=minexp, to=maxexp, by=increment) #creating a vector that ranges evenly across range of exposure
  
  preddf = data.frame(cbind(range_exp,mage, income_pers, parity, hiv_combined_i, csex, time_dtp4, time )) #dataframe with all values for prediction
  
  preddf$range_exp = as.numeric(as.character(preddf$range_exp))
  preddf$income_pers = as.numeric(as.character(preddf$income_pers))
  preddf$parity = as.integer(as.character(preddf$parity))
  preddf$time_dtp4 = as.integer(as.character(preddf$time_dtp4))
  preddf$time = as.numeric(as.character(preddf$time))
  preddf$mage = as.numeric(as.character(preddf$mage))
  
  
  colnames(preddf) = c(expname, colnames(preddf)[2:length(colnames(preddf))]) 
  
  return(preddf)
  
}


## Get Linear Main Effects with GEE (no interaction with time) ##
runGEE = function(outcome, exposures, adjustment, df) {
  
  iloop=list()
  for (i in seq_along(exposures)) {
    
    f_intx =  as.formula(paste(outcome, "~", exposures[i], "+", "time", "+", paste(adjustment, collapse="+")))
    
    m_intx = gee(data=df, formula=f_intx, id=hsn, corstr="exchangeable")
    
    coef_intx = as.data.frame(coefficients(summary(m_intx))) #extracting relevent coefficients 
    coef_intx = coef_intx[rownames(coef_intx) ==exposures[i],]
    coef_intx$exposure = rownames(coef_intx)
    coef_intx$outcome = outcome
    coef_intx$est_robust95ci = paste(round(coef_intx$Estimate,3), " [", 
                                   round((coef_intx$Estimate-(1.96*coef_intx$`Robust S.E.`)),3),
                                   ", ",
                                   round((coef_intx$Estimate+(1.96*coef_intx$`Robust S.E.`)),3),
                                   "]", sep="")
    coef_intx = coef_intx[, c("outcome","exposure", "est_robust95ci")]
    
    #coef_intx[4,] = rep(" ", 3) #just adding a row of blanks to make thigs easier to read
    rownames(coef_intx) = c()
    
    iloop[i] = list(coef_intx)
    
  }
  ibound = do.call("rbind", iloop)
  return(ibound)
  
}


## Get Binary (RR) Effects with GEE (no interaction with time) ##
runGEEbinary = function(outcome, exposures, adjustment, df) {
  
  iloop=list()
  for (i in seq_along(exposures)) {
    
    f_intx =  as.formula(paste(outcome, "~", exposures[i], "+", "time", "+", paste(adjustment, collapse="+")))
    
    m_intx = gee(data=df, formula=f_intx, id=hsn, corstr="exchangeable", family=poisson)
    
    coef_intx = as.data.frame(coefficients(summary(m_intx))) #extracting relevent coefficients 
    coef_intx = coef_intx[rownames(coef_intx) ==exposures[i],]
    coef_intx$exposure = rownames(coef_intx)
    coef_intx$outcome = outcome
    coef_intx$est_robust95ci = paste(round(exp(coef_intx$Estimate),3), " [", 
                                     round(exp(coef_intx$Estimate-(1.96*coef_intx$`Robust S.E.`)),3),
                                     ", ",
                                     round(exp(coef_intx$Estimate+(1.96*coef_intx$`Robust S.E.`)),3),
                                     "]", sep="")
    coef_intx = coef_intx[, c("outcome","exposure", "est_robust95ci")]
    
    #coef_intx[4,] = rep(" ", 3) #just adding a row of blanks to make thigs easier to read
    rownames(coef_intx) = c()
    
    iloop[i] = list(coef_intx)
    
  }
  ibound = do.call("rbind", iloop)
  return(ibound)
  
}


## Test Interaction on Categorical Outcomes ##
runGEEIntxBinary = function(outcome, exposures, adjustment, modifier, df) {
  
  out=list()
  for (i in seq_along(exposures)) {
    
    f_intx =  as.formula(paste(outcome, "~", exposures[i], "+", "time", "+", modifier, "+", modifier, "*", exposures[i], "+", 
                               paste(adjustment, collapse="+")))
    
    m_intx = gee(data=df, formula=f_intx, id=hsn, corstr="exchangeable", family=poisson)
    
    coef_intx = as.data.frame(coefficients(summary(m_intx))) #extracting relevent coefficients 
    coef_intx = coef_intx[grep(exposures[i], rownames(coef_intx)),] #restricting to only estimates pertaining to the exposure
    coef_intx = coef_intx[grep(modifier, rownames(coef_intx)), ] #now restricting to only the estimate of the interaction between exposure and modifier
    
    est_robust95ci = paste(round(exp(coef_intx$Estimate),3), " [", 
                           round(exp(coef_intx$Estimate-(1.96*coef_intx$`Robust S.E.`)),3),
                           ", ",
                           round(exp(coef_intx$Estimate+(1.96*coef_intx$`Robust S.E.`)),3),
                           "]", sep="")
    
    pval = 2*pnorm(abs(coef_intx$Estimate/coef_intx$`Robust S.E.`), lower.tail=F)
    
    nobs = nobs(m_intx)
    predictor = exposures[i]
    dv = outcome
    
    out[i] = list(
      as.data.frame(cbind(
        predictor, dv, nobs, est_robust95ci, pval 
      ))
    )
    rm(predictor, dv, nobs, pval, est_robust95ci)
   
    
    
  }
  
  output = do.call("rbind", out)
  
  output$pval = as.numeric(as.character(output$pval))
  output$siglabel = ""
  output$siglabel[output$pval<.10] = "*"
  
  return(output)
  
}


## Test Interaction on Continuous Outcomes ##
runGEEIntxCont = function(outcome, exposures, adjustment, modifier, df) {
  
  out=list()
  for (i in seq_along(exposures)) {
    
    f_intx =  as.formula(paste(outcome, "~", exposures[i], "+", "time", "+", modifier, "+", modifier, "*", exposures[i], "+", 
                               paste(adjustment, collapse="+")))
    
    m_intx = gee(data=df, formula=f_intx, id=hsn, corstr="exchangeable")
    
    coef_intx = as.data.frame(coefficients(summary(m_intx))) #extracting relevent coefficients 
    coef_intx = coef_intx[grep(exposures[i], rownames(coef_intx)),] #restricting to only estimates pertaining to the exposure
    coef_intx = coef_intx[grep(modifier, rownames(coef_intx)), ] #now restricting to only the estimate of the interaction between exposure and modifier
    
    est_robust95ci = paste(round(exp(coef_intx$Estimate),3), " [", 
                           round(exp(coef_intx$Estimate-(1.96*coef_intx$`Robust S.E.`)),3),
                           ", ",
                           round(exp(coef_intx$Estimate+(1.96*coef_intx$`Robust S.E.`)),3),
                           "]", sep="")
    
    pval = 2*pnorm(abs(coef_intx$Estimate/coef_intx$`Robust S.E.`), lower.tail=F)
    
    nobs = nobs(m_intx)
    predictor = exposures[i]
    dv = outcome
    
    out[i] = list(
      as.data.frame(cbind(
        predictor, dv, nobs, est_robust95ci, pval 
      ))
    )
    rm(predictor, dv, nobs, pval, est_robust95ci)
    
    
    
  }
  
  output = do.call("rbind", out)
  
  output$pval = as.numeric(as.character(output$pval))
  output$siglabel = ""
  output$siglabel[output$pval<.10] = "*"
  
  return(output) 
  
}


## Getting Overall Estimates for Categorical Outcome (for plotting) ##
getGEEOverallEsts = function(outcome, exposures, adjustment,  df) {
  
  df1=df

  out=list()
  for (i in seq_along(exposures)) {
    
    f_intx =  as.formula(paste(outcome, "~", exposures[i], "+", "time", "+", paste(adjustment, collapse="+")))
    
    
    ## For first level of stratifying variable ##
    m_intx = gee(data=df1, formula=f_intx, id=hsn, corstr="exchangeable", family=poisson)
    
    coef_intx = as.data.frame(coefficients(summary(m_intx))) #extracting relevent coefficients 
    coef_intx = coef_intx[rownames(coef_intx) ==exposures[i],]
    
    rr = exp(coef_intx$Estimate) #getting the risk ratio and robust 95% CI
    lo_rr = exp(coef_intx$Estimate - (1.96*coef_intx$`Robust S.E.`))
    hi_rr = exp(coef_intx$Estimate + (1.96*coef_intx$`Robust S.E.`))
    
    predictor = exposures[i]
    dv = outcome

    out1 = as.data.frame(cbind(predictor, dv, rr, lo_rr, hi_rr))
    
    
    ## Binding together both levels ##
    out[i] = list(
      as.data.frame(rbind(out1))
    )
    
  }
  output = do.call("rbind", out)
  
  output$rr = as.numeric(as.character(output$rr)) #making sure all numeric variables are coded as such
  output$lo_rr = as.numeric(as.character(output$lo_rr))
  output$hi_rr = as.numeric(as.character(output$hi_rr))
  
  return(output)
  
  
}


## Getting Overall Estimates for Categorical Outcome (for plotting) ##
getGEEOverallEstsCont = function(outcome, exposures, adjustment,  df) {
  
  df1=df
  
  out=list()
  for (i in seq_along(exposures)) {
    
    f_intx =  as.formula(paste(outcome, "~", exposures[i], "+", "time", "+", paste(adjustment, collapse="+")))
    
    
    ## For first level of stratifying variable ##
    m_intx = gee(data=df1, formula=f_intx, id=hsn, corstr="exchangeable")
    
    coef_intx = as.data.frame(coefficients(summary(m_intx))) #extracting relevent coefficients 
    coef_intx = coef_intx[rownames(coef_intx) ==exposures[i],]
    
    rr = 100*((10^(coef_intx$Estimate))-1) #getting the percents for the beta? and robust 95% CI
    lo_rr = 100*((10^(coef_intx$Estimate - (1.96*coef_intx$`Robust S.E.`)))-1) 
    hi_rr = 100*((10^(coef_intx$Estimate + (1.96*coef_intx$`Robust S.E.`)))-1)
    
    predictor = exposures[i]
    dv = outcome
    
    out1 = as.data.frame(cbind(predictor, dv, rr, lo_rr, hi_rr))
    
    
    ## Binding together both levels ##
    out[i] = list(
      as.data.frame(rbind(out1))
    )
    
  }
  output = do.call("rbind", out)
  
  output$rr = as.numeric(as.character(output$rr)) #making sure all numeric variables are coded as such
  output$lo_rr = as.numeric(as.character(output$lo_rr))
  output$hi_rr = as.numeric(as.character(output$hi_rr))
  
  return(output)
  
  
}


## Getting Stratified Estimates for Categorical Outcome ##
getGEEStratEsts = function(outcome, exposures, adjustment, l1, l2, stratby, df) {
  
  s = df[, colnames(df) == stratby]
  
  df1 = df[s==l1 & !is.na(s),]
  df2 = df[s==l2 & !is.na(s),]
  
  
  out=list()
  for (i in seq_along(exposures)) {
    
    f_intx =  as.formula(paste(outcome, "~", exposures[i], "+", "time", "+", paste(adjustment, collapse="+")))
    
    
    ## For first level of stratifying variable ##
    m_intx = gee(data=df1, formula=f_intx, id=hsn, corstr="exchangeable", family=poisson)
    
    coef_intx = as.data.frame(coefficients(summary(m_intx))) #extracting relevent coefficients 
    coef_intx = coef_intx[rownames(coef_intx) ==exposures[i],]
    
    rr = exp(coef_intx$Estimate) #getting the risk ratio and robust 95% CI
    lo_rr = exp(coef_intx$Estimate - (1.96*coef_intx$`Robust S.E.`))
    hi_rr = exp(coef_intx$Estimate + (1.96*coef_intx$`Robust S.E.`))
    
    predictor = exposures[i]
    dv = outcome
    level = l1
    
    out1 = as.data.frame(cbind(predictor, dv, level, rr, lo_rr, hi_rr))
    
    rm(predictor, dv, level, rr, lo_rr, hi_rr, coef_intx, m_intx)
    
    
    ## For second level of stratifying variable ##
    m_intx = gee(data=df2, formula=f_intx, id=hsn, corstr="exchangeable", family=poisson)
    
    coef_intx = as.data.frame(coefficients(summary(m_intx))) #extracting relevent coefficients 
    coef_intx = coef_intx[rownames(coef_intx) ==exposures[i],]
    
    rr = exp(coef_intx$Estimate) #getting the risk ratio and robust 95% CI
    lo_rr = exp(coef_intx$Estimate - (1.96*coef_intx$`Robust S.E.`))
    hi_rr = exp(coef_intx$Estimate + (1.96*coef_intx$`Robust S.E.`))
    
    predictor = exposures[i]
    dv = outcome
    level = l2
    
    out2 = as.data.frame(cbind(predictor, dv, level, rr, lo_rr, hi_rr))
    
    rm(predictor, dv, level, rr, lo_rr, hi_rr)
    
    
    
    ## Binding together both levels ##
    out[i] = list(
      as.data.frame(rbind(out1, out2))
    )
    
    
    
    
  }
  output = do.call("rbind", out)
  
  output$rr = as.numeric(as.character(output$rr)) #making sure all numeric variables are coded as such
  output$lo_rr = as.numeric(as.character(output$lo_rr))
  output$hi_rr = as.numeric(as.character(output$hi_rr))
  
  return(output)
  
  
}




## Getting Stratified Estimates for Continuous Outcome ##
getGEEStratEstsCont = function(outcome, exposures, adjustment, l1, l2, stratby, df) {
  
  s = df[, colnames(df) == stratby]
  
  df1 = df[s==l1 & !is.na(s),]
  df2 = df[s==l2 & !is.na(s),]
  
  
  out=list()
  for (i in seq_along(exposures)) {
    
    f_intx =  as.formula(paste(outcome, "~", exposures[i], "+", "time", "+", paste(adjustment, collapse="+")))
    
    
    ## For first level of stratifying variable ##
    m_intx = gee(data=df1, formula=f_intx, id=hsn, corstr="exchangeable")
    
    coef_intx = as.data.frame(coefficients(summary(m_intx))) #extracting relevent coefficients 
    coef_intx = coef_intx[rownames(coef_intx) ==exposures[i],]
    
    rr = 100*(exp(coef_intx$Estimate)-1) #getting the percent (?) beta and robust 95% CI
    lo_rr = 100*(exp(coef_intx$Estimate - (1.96*coef_intx$`Robust S.E.`))-1)
    hi_rr = 100*(exp(coef_intx$Estimate + (1.96*coef_intx$`Robust S.E.`))-1)
    
    predictor = exposures[i]
    dv = outcome
    level = l1
    
    out1 = as.data.frame(cbind(predictor, dv, level, rr, lo_rr, hi_rr))
    
    rm(predictor, dv, level, rr, lo_rr, hi_rr, coef_intx, m_intx)
    
    
    ## For second level of stratifying variable ##
    m_intx = gee(data=df2, formula=f_intx, id=hsn, corstr="exchangeable")
    
    coef_intx = as.data.frame(coefficients(summary(m_intx))) #extracting relevent coefficients 
    coef_intx = coef_intx[rownames(coef_intx) ==exposures[i],]
    
    rr = 100*(exp(coef_intx$Estimate)-1) #getting the beta and robust 95% CI
    lo_rr = 100*(exp(coef_intx$Estimate - (1.96*coef_intx$`Robust S.E.`))-1)
    hi_rr = 100*(exp(coef_intx$Estimate + (1.96*coef_intx$`Robust S.E.`))-1)
    
    predictor = exposures[i]
    dv = outcome
    level = l2
    
    out2 = as.data.frame(cbind(predictor, dv, level, rr, lo_rr, hi_rr))
    
    rm(predictor, dv, level, rr, lo_rr, hi_rr)
    
    
    
    ## Binding together both levels ##
    out[i] = list(
      as.data.frame(rbind(out1, out2))
    )
    
    
    
    
  }
  output = do.call("rbind", out)
  
  output$rr = as.numeric(as.character(output$rr)) #making sure all numeric variables are coded as such
  output$lo_rr = as.numeric(as.character(output$lo_rr))
  output$hi_rr = as.numeric(as.character(output$hi_rr))
  
  colnames(output) = c("predictor", "dv", "level", "beta", "lo_ci", "hi_ci")
  
  return(output)
  
  
}


## For Running GEE and Outputting Predicted Vector ##
runGEEAndPredict = function(data, indices, f, preddf) {
  
  dt = data[indices,] #indices for the boot function
  
  zmodel = Zelig::zelig(f, model = "normal.gee", id="hsn", data=dt, corstr="exchangeable") #using zelig to run GEE so can get the predicted more easily
  
  zpred = c(do.call("cbind", Zelig::predict(zmodel, preddf)))
  
  return(zpred)
  
}


## For Running GEE With Binary Outcome and Outputting Predicted Vector, Comparing to 1st Quartile ##
runGEEAndPredictBinary = function(data, indices, f, preddf, q25) {
  
  dt = data[indices,] #indices for the boot function
  
  zmodel = Zelig::zelig(f, model = "poisson.gee", id="hsn", data=dt, corstr="exchangeable") #using zelig to run GEE so can get the predicted more easily
  
  zpred = c(do.call("cbind", Zelig::predict(zmodel, preddf))) #getting predicted log(odds) at each value of the exposure

  zpred2 = cbind(preddf, zpred)
  zpred2$ref = zpred2$zpred[DescTools::Closest(zpred2[,1], q25, which=T)] #making the referent value the predicted value at the 25th percentile cutoff in the exposure
  
  zpred2$logor = zpred2$zpred - zpred2$ref #computing the log(odds ratio) comparing the predicted values to the value at the 25th percentile cutoff of the exposure
  
  out = c(exp(zpred2$logor))
  
  for (i in 1:length(indices)) {
    print(paste(i, "/", length(indices)))
    Sys.sleep(0.01)
    flush.console()
  }

  return(out)
  
}


## For Running GEE With Binary Outcome (WITH LOGIT) and Outputting Predicted Vector, Comparing to 1st Quartile ##
runGEEAndPredictBinaryLogit = function(data, indices, f, preddf, q25) {
  
  dt = data[indices,] #indices for the boot function
  
  zmodel = Zelig::zelig(f, model = "logit.gee", id="hsn", data=dt, corstr="exchangeable") #using zelig to run GEE so can get the predicted more easily
  
  zpred = c(do.call("cbind", Zelig::predict(zmodel, preddf))) #getting predicted log(odds) at each value of the exposure
  
  zpred2 = cbind(preddf, zpred)
  zpred2$ref = zpred2$zpred[DescTools::Closest(zpred2[,1], q25, which=T)] #making the referent value the predicted value at the 25th percentile cutoff in the exposure
  
  zpred2$logor = zpred2$zpred - zpred2$ref #computing the log(odds ratio) comparing the predicted values to the value at the 25th percentile cutoff of the exposure
  
  out = c(exp(zpred2$logor))
  
  # for (i in 1:length(indices)) {
  #   print(paste(i, "/", length(indices)))
  #   Sys.sleep(0.01)
  #   flush.console()
  # }
  # 
  return(out)
  
}


## Function to Take Predicted Values and Boot CI (and append exposure values) ##
bootCIsGEE = function(results, range_exp) {
  
  pred_values = list()
  for (i in 1:ncol(results$t)) {
    
    r = boot::boot.ci(results, type="basic",index=i)
    estq = r[[2]]
    cis = r[[4]]
    colnames(cis) = NULL
    lcl = cis[1,4]
    ucl = cis[1,5]
    rdf = data.frame(estq, lcl, ucl)
    pred_values[i] = list(rdf)
  }
  pred_values_df = do.call("rbind", pred_values)
  pred_values_df$exp = range_exp
  
  return(pred_values_df)
  
}


## Function to Take Predicted Values and Boot CI (and append exposure values) for the Binary Comparison ##
bootCIsGEEBinary = function(results, range_exp, q25) {
  
  breakcol = DescTools::Closest(range_exp, q25, which=T) #basically the CI will be (1,1) at the referent point, so taking that in to consideration
  
  
  pred_values1 = list()
  for (i in 1:(breakcol-1)) {
    
    r = boot::boot.ci(results, type="basic",index=i)
    estq = r[[2]]
    cis = r[[4]]
    colnames(cis) = NULL
    lcl = cis[1,4]
    ucl = cis[1,5]
    rdf = data.frame(estq, lcl, ucl)
    pred_values1[i] = list(rdf)
  }
  pred_values_df1 = do.call("rbind", pred_values1)
  
  
  pred_values2 = list()
  for (i in (breakcol+1):ncol(results$t)) {
    
    r = boot::boot.ci(results, type="basic",index=i)
    estq = r[[2]]
    cis = r[[4]]
    colnames(cis) = NULL
    lcl = cis[1,4]
    ucl = cis[1,5]
    rdf = data.frame(estq, lcl, ucl)
    pred_values2[i] = list(rdf)
  }
  pred_values_df2 = do.call("rbind", pred_values2)

  pred_values_df3 = data.frame(estq = 1, lcl=1, ucl=1)
  
  pred_values_df = rbind(pred_values_df1, pred_values_df3, pred_values_df2)
  pred_values_df$exp = range_exp
  
  return(pred_values_df)
  
}


## Function to Save Two Stratified Plots for CONTINUOUS outcomes (First strata, then both) ##
produceStratPlot = function(d1, d2,d3, name, lims=c(-35,35), cuts=c(-20,0,20), ast = 25) {
  
  plot1=
    ggplot() + 
    geom_pointrange(data=d1,
                    aes(x=predictor, y=beta, ymin=lo_ci, ymax=hi_ci, fill=level), pch=21,
                    position=position_dodge(width=0.6), size=1.5) +
    #geom_text(data=d2, aes(x=predictor, y=2, label=siglabel), size=20) + 
    facet_grid(dv~.) + 
    scale_y_continuous(breaks=cuts) +
    coord_cartesian(ylim=lims) +
    scale_fill_manual(values=c("black", "white")) + 
    scale_x_discrete(labels = c(expression(paste(italic("p,p'"), "-DDT")),
                                expression(paste(italic("p,p'"), "-DDE")),
                                expression(paste(italic("cis"), "-DBCA")),
                                expression(paste(italic("cis"), "-DCCA")),
                                expression(paste(italic("trans"), "-DCCA")),
                                "3-PBA")) +
    geom_hline(yintercept=0, linetype='dashed') +
    labs(
      y="Percent Difference in Antibody Titers\n",
      x="")  +
    theme_light() +
    theme(
      text = element_text(family = "Century Gothic"),
      legend.title = element_blank(),
      legend.text = element_text(size=18),
      legend.key.height = unit(3, "cm"),
      legend.key.width = unit(3, "cm"),
      
      legend.position = "top",
      legend.justification = c(0,0),
      legend.direction = "horizontal",
      legend.background = element_rect(fill="transparent"),
      legend.key = element_rect(fill="transparent"),
      
      axis.text = element_text(size=20),
      axis.text.x = element_text(size=20),
      axis.title = element_text(size=25),
      plot.caption = element_text(size=20, hjust=0),
      strip.background = element_blank(),
      strip.text = element_text(size=20, colour="black")
    ) 
  ggsave(plot1, file=paste(name, "_1.png", sep=""), height=15, width =12)
  
  plot2=
    ggplot() + 
    geom_pointrange(data=d3,
                    aes(x=predictor, y=beta, ymin=lo_ci, ymax=hi_ci, fill=level), pch=21,
                    position=position_dodge(width=0.6), size=1.5) +
    geom_text(data=d2, aes(x=predictor, y=ast, label=siglabel), size=20) + 
    facet_grid(dv~.) + 
    scale_y_continuous(breaks=cuts) +
    #coord_cartesian(ylim=c(0.2,4.2)) +
    coord_cartesian(ylim=lims) +
    scale_fill_manual(values=c("black", "white")) + 
    scale_x_discrete(labels = c(expression(paste(italic("p,p'"), "-DDT")),
                                expression(paste(italic("p,p'"), "-DDE")),
                                expression(paste(italic("cis"), "-DBCA")),
                                expression(paste(italic("cis"), "-DCCA")),
                                expression(paste(italic("trans"), "-DCCA")),
                                "3-PBA")) +
    geom_hline(yintercept=0, linetype='dashed') +
    labs(
      y="Percent Difference in Antibody Titers\n",
      x="")   +
    theme_light() +
    theme(
      text = element_text(family = "Century Gothic"),
      legend.title = element_blank(),
      legend.text = element_text(size=18),
      legend.key.height = unit(3, "cm"),
      legend.key.width = unit(3, "cm"),
      
      legend.position = "top",
      legend.justification = c(0,0),
      legend.direction = "horizontal",
      legend.background = element_rect(fill="transparent"),
      legend.key = element_rect(fill="transparent"),
      
      axis.text = element_text(size=20),
      axis.text.x = element_text(size=20),
      axis.title = element_text(size=25),
      plot.caption = element_text(size=20, hjust=0),
      strip.background = element_blank(),
      strip.text = element_text(size=20, colour="black")
    ) 
  ggsave(plot2, file=paste(name, "_2.png", sep=""), height=15, width =12)
  
}


## Function to Save Three Overall Plots for CONTINUOUS (First strata, then two, then all) -- stratified by immune outcome ##
produceStratPlot3 = function(d1,d3,d4, name, cuts=c(-0.3,0.3)) {
  
  plot1=
    ggplot() + 
    geom_pointrange(data=d1,
                    aes(x=predictor, y=rr, ymin=lo_rr, ymax=hi_rr, fill=dv), pch=21,
                    position=position_dodge(width=0.6), size=1.5) +
    #geom_text(data=d2, aes(x=predictor, y=2, label=siglabel), size=20) + 
    #scale_y_continuous(breaks=c(0.25, 0.5, 1, 2, 4)) +
    coord_cartesian(ylim=cuts) +
    scale_fill_manual(values=c(ghibli_palettes$MononokeMedium[c(3)])) + 
    scale_x_discrete(labels = c(expression(paste(italic("p,p'"), "-DDT")),
                                expression(paste(italic("p,p'"), "-DDE")),
                                expression(paste(italic("cis"), "-DBCA")),
                                expression(paste(italic("cis"), "-DCCA")),
                                expression(paste(italic("trans"), "-DCCA")),
                                "3-PBA")) +
    geom_hline(yintercept=0, linetype='dashed') +
    labs(
      y="Percent Difference in Antibody Titers\n",
      x="")  +
    theme_light() +
    theme(
      text = element_text(family = "Century Gothic"),
      legend.title = element_blank(),
      legend.text = element_text(size=18),
      legend.key.height = unit(3, "cm"),
      legend.key.width = unit(3, "cm"),
      
      legend.position = c(0,0.75),
      legend.justification = c(0,0),
      legend.direction = "horizontal",
      legend.background = element_rect(fill="transparent"),
      legend.key = element_rect(fill="transparent"),
      
      axis.text = element_text(size=20),
      axis.text.x = element_text(size=20),
      axis.title = element_text(size=25),
      plot.caption = element_text(size=20, hjust=0),
      strip.background = element_blank(),
      strip.text = element_text(size=20, colour="black")
    ) 
  ggsave(plot1, file=paste(name, "_1.png", sep=""), height=8, width =12)
  
  plot2=
    ggplot() + 
    geom_pointrange(data=d3,
                    aes(x=predictor, y=rr, ymin=lo_rr, ymax=hi_rr, fill=dv), pch=21,
                    position=position_dodge(width=0.6), size=1.5) +
    #geom_text(data=d2, aes(x=predictor, y=2, label=siglabel), size=20) + 
    #scale_y_continuous(breaks=c(0.25, 0.5, 1, 2, 4)) +
    coord_cartesian(ylim=cuts) +
    scale_fill_manual(values=c(ghibli_palettes$MononokeMedium[c(3,5)])) + 
    scale_x_discrete(labels = c(expression(paste(italic("p,p'"), "-DDT")),
                                expression(paste(italic("p,p'"), "-DDE")),
                                expression(paste(italic("cis"), "-DBCA")),
                                expression(paste(italic("cis"), "-DCCA")),
                                expression(paste(italic("trans"), "-DCCA")),
                                "3-PBA")) +
    geom_hline(yintercept=0, linetype='dashed') +
    labs(
      y="Percent Difference in Antibody Titers\n",
      x="")   +
    theme_light() +
    theme(
      text = element_text(family = "Century Gothic"),
      legend.title = element_blank(),
      legend.text = element_text(size=18),
      legend.key.height = unit(3, "cm"),
      legend.key.width = unit(3, "cm"),
      
      
      legend.position = c(0,0.75),
      legend.justification = c(0,0),
      legend.direction = "horizontal",
      legend.background = element_rect(fill="transparent"),
      legend.key = element_rect(fill="transparent"),
      
      axis.text = element_text(size=20),
      axis.text.x = element_text(size=20),
      axis.title = element_text(size=25),
      plot.caption = element_text(size=20, hjust=0),
      strip.background = element_blank(),
      strip.text = element_text(size=20, colour="black")
    ) 
  ggsave(plot2, file=paste(name, "_2.png", sep=""), height=8, width =12)
  
  
  plot3=
    ggplot() + 
    geom_pointrange(data=d4,
                    aes(x=predictor, y=rr, ymin=lo_rr, ymax=hi_rr, fill=dv), pch=21,
                    position=position_dodge(width=0.6), size=1.5) +
   # geom_text(data=d2, aes(x=predictor, y=2, label=siglabel), size=20) + 
    #scale_y_continuous(breaks=c(0.25, 0.5, 1, 2, 4)) +
    coord_cartesian(ylim=cuts) +
    scale_fill_manual(values=c(ghibli_palettes$MononokeMedium[c(3,5,7)])) + 
    scale_x_discrete(labels = c(expression(paste(italic("p,p'"), "-DDT")),
                                expression(paste(italic("p,p'"), "-DDE")),
                                expression(paste(italic("cis"), "-DBCA")),
                                expression(paste(italic("cis"), "-DCCA")),
                                expression(paste(italic("trans"), "-DCCA")),
                                "3-PBA")) +
    geom_hline(yintercept=0, linetype='dashed') +
    labs(
      y="Percent Difference in Antibody Titers\n",
      x="")  +
    theme_light() +
    theme(
      text = element_text(family = "Century Gothic"),
      legend.title = element_blank(),
      legend.text = element_text(size=18),
      legend.key.height = unit(3, "cm"),
      legend.key.width = unit(3, "cm"),
      
      legend.position = c(0,0.75),
      legend.justification = c(0,0),
      legend.direction = "horizontal",
      legend.background = element_rect(fill="transparent"),
      legend.key = element_rect(fill="transparent"),
      
      axis.text = element_text(size=20),
      axis.text.x = element_text(size=20),
      axis.title = element_text(size=25),
      plot.caption = element_text(size=20, hjust=0),
      strip.background = element_blank(),
      strip.text = element_text(size=20, colour="black")
    ) 
  ggsave(plot3, file=paste(name, "_3.png", sep=""), height=8, width =12)
  
}




## >TABLE 1: DESCRIPTIVES ##########################################################################

dr_t1 = dr[,colnames(dr) %in% c("hsn", t1_vars)] #creating dataframe with only one row per HSN for the table 1 variables
dr_t1 = unique(dr_t1)

t1 = catTable1(df=dr_t1, varlist=t1_vars)
write.xlsx(t1, file="analyses_2019-10-30.xlsx", sheetName="t1", append=T)


## >TABLE 2: EXPOSURE DESCRIPTIVES #################################################################

dr_t2 = dr[,colnames(dr) %in% c("hsn", exp_list_notlog,lodloq_vars)]
dr_t2 = unique(dr_t2)

t2 = makeT2(df=dr_t2, varlist=exp_list_notlog)
write.xlsx(t2, file="analyses_2019-10-30.xlsx", sheetName="t2", append=T)

## For getting detection/quantification frequencies ##
write.xlsx(summary(dr_t2[,colnames(dr_t2) %in% lodloq_vars]), file="analyses_2019-10-30.xlsx",  sheetName="LODs", append=T)

## For seeing if the min is below the detection/quantification freq ## 
dr$opdde_lod3[dr$n_opdde_i == min(dr$n_opdde_i, na.rm=T)]
dr$ppdde_lod3[dr$n_ppdde_i == min(dr$n_ppdde_i, na.rm=T)]
dr$opddt_lod3[dr$n_opddt_i == min(dr$n_opddt_i, na.rm=T)]
dr$ppddt_lod3[dr$n_ppddt_i == min(dr$n_ppddt_i, na.rm=T)]
dr$cisdbca_lod[dr$cisdcca_sg == min(dr$cisdcca_sg, na.rm=T)]
dr$cisdbca_loq[dr$cisdcca_sg == min(dr$cisdcca_sg, na.rm=T)]
dr$cisdcca_lod[dr$cisdcca_sg == min(dr$cisdcca_sg, na.rm=T)]
dr$cisdcca_loq[dr$cisdcca_sg == min(dr$cisdcca_sg, na.rm=T)]
dr$transdcca_lod[dr$transdcca_sg == min(dr$transdcca_sg, na.rm=T)]
dr$transdcca_loq[dr$transdcca_sg == min(dr$transdcca_sg, na.rm=T)]
dr$pba_lod[dr$pba_sg == min(dr$pba_sg, na.rm=T)]
dr$pba_loq[dr$pba_sg == min(dr$pba_sg, na.rm=T)]



## >MAIN EFFECTS @ EACH TIME ##############################################################################


## >>At 3.5 Years -----------------------------------------------------------------------------------------

## Getting the Dataframes for 3.5 Year Timepoint ##
dr_meas_3yr = dr_meas[dr_meas$time==0 & !is.na(dr_meas$time) & !is.na(dr_meas$logmev),]
doubleID(dr_meas_3yr) #no doubled ID

dr_tet_3yr = dr_tet[dr_tet$time==0 & !is.na(dr_tet$time) & (!is.na(dr_tet$logtet) | !is.na(dr_tet$loghib)),]
doubleID(dr_tet_3yr) #no doubled ID


## Measles ##
me_meas_point_3 = getEstimatesME(df=dr_meas_3yr, outcome = "notprot_meas", exposures=exp_list, adjustment=adjset_meas)
me_meas_point_3 = cleanEstimatesME(me_meas_point_3)

me_meas_lin_3 = getLinEstimatesME(df=dr_meas_3yr, outcome="logmev", exposures=exp_list, adjustment=adjset_meas)

me_meas_3 = cbind(me_meas_point_3, me_meas_lin_3)


## Tetatanus ##
me_tet_point_3 = getEstimatesME(df=dr_tet_3yr, outcome = "notprot_tet", exposures=exp_list, adjustment=adjset_tet)
me_tet_point_3 = cleanEstimatesME(me_tet_point_3)

me_tet_lin_3 = getLinEstimatesME(df=dr_tet_3yr, outcome="logtet", exposures=exp_list, adjustment=adjset_tet)

me_tet_3 = cbind(me_tet_point_3, me_tet_lin_3)


## HIB ##
me_hib_point_3 = getEstimatesME(df=dr_tet_3yr, outcome = "notprot_hib", exposures=exp_list, adjustment=adjset_tet)
me_hib_point_3 = cleanEstimatesME(me_hib_point_3)

me_hib_lin_3 = getLinEstimatesME(df=dr_tet_3yr, outcome="loghib", exposures=exp_list, adjustment=adjset_tet)

me_hib_3 = cbind(me_hib_point_3, me_hib_lin_3)


## Binding all together ##
mes_3 = rbind(me_meas_3, me_tet_3, me_hib_3)
write.xlsx(mes_3, file="analyses_2019-10-30.xlsx",  sheetName="MEs3yr", append=T)



## >>At 5 Years ----------------------------------------------------------------------------

## Getting the Dataframes for 3.5 Year Timepoint ##
dr_meas_5yr = dr_meas[dr_meas$time==1 & !is.na(dr_meas$logmev),] #hmm why is this number the same as the 3.5 year?? Need to think about this
doubleID(dr_meas_5yr) #no doubled ID

dr_tet_5yr = dr_tet[dr_tet$time==1 & (!is.na(dr_tet$logtet) | !is.na(dr_tet$loghib)),]
doubleID(dr_tet_5yr) #no doubled ID


## Measles ##
me_meas_point_5 = getEstimatesME(df=dr_meas_5yr, outcome = "notprot_meas", exposures=exp_list, adjustment=adjset_meas)
me_meas_point_5 = cleanEstimatesME(me_meas_point_5)

me_meas_lin_5 = getLinEstimatesME(df=dr_meas_5yr, outcome="logmev", exposures=exp_list, adjustment=adjset_meas)

me_meas_5 = cbind(me_meas_point_5, me_meas_lin_5)


## Tetatanus ##
me_tet_point_5 = getEstimatesME(df=dr_tet_5yr, outcome = "notprot_tet", exposures=exp_list, adjustment=adjset_tet)
me_tet_point_5 = cleanEstimatesME(me_tet_point_5)

me_tet_lin_5 = getLinEstimatesME(df=dr_tet_5yr, outcome="logtet", exposures=exp_list, adjustment=adjset_tet)

me_tet_5 = cbind(me_tet_point_5, me_tet_lin_5)


## HIB ##
me_hib_point_5 = getEstimatesME(df=dr_tet_5yr, outcome = "notprot_hib", exposures=exp_list, adjustment=adjset_tet)
me_hib_point_5 = cleanEstimatesME(me_hib_point_5)

me_hib_lin_5 = getLinEstimatesME(df=dr_tet_5yr, outcome="loghib", exposures=exp_list, adjustment=adjset_tet)

me_hib_5 = cbind(me_hib_point_5, me_hib_lin_5)


## Binding all together ##
mes_5 = rbind(me_meas_5, me_tet_5, me_hib_5)
write.xlsx(mes_5, file="analyses_2019-10-30.xlsx",  sheetName="MEs5yr", append=T)


## >CONTINUOUS STRATIFIED EFFECTS (EMM) @ EACH TIME POINT ##################################################

## >>Child Sex in 3.5 Year Data --------------------------------------------------------------------------

## Running Interaction Tests ##
meas_csex = runMEIntxCont(outcome="logmev", exposures=exp_list, adjustment=adjset_meas[!(adjset_meas=="csex")], modifier="csex", df=dr_meas_3yr)
tet_csex = runMEIntxCont(outcome="logtet", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="csex")], modifier="csex", df=dr_tet_3yr)
hib_csex = runMEIntxCont(outcome="loghib", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="csex")], modifier="csex", df=dr_tet_3yr)

intxtests_csex_cont = rbind(meas_csex, tet_csex, hib_csex)
rm(meas_csex, tet_csex, hib_csex)

intxtests_csex_cont$dv = factor(intxtests_csex_cont$dv,
                                levels = c("logmev", "logtet", "loghib"),
                                labels = c("Measles", "Tetanus", "HIB"))

## Getting Stratified Estimates ##
meas_csex = getMEStratEstsCont(outcome="logmev", exposures=exp_list, adjustment=adjset_meas[!(adjset_meas=="csex")], 
                                stratby="csex", l1="Boy", l2="Girl", df=dr_meas_3yr)

tet_csex = getMEStratEstsCont(outcome="logtet", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="csex")], 
                               stratby="csex", l1="Boy", l2="Girl",  df=dr_tet_3yr)

hib_csex = getMEStratEstsCont(outcome="loghib", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="csex")], 
                               stratby="csex", l1="Boy", l2="Girl",  df=dr_tet_3yr)

mod_csex_cont = rbind(meas_csex, tet_csex, hib_csex)

mod_csex_cont$dv = factor(mod_csex_cont$dv,
                          levels = c("logmev", "logtet", "loghib"),
                          labels = c("Measles", "Tetanus", "HIB"))

mod_csex_cont$level = factor(mod_csex_cont$level, 
                             levels = c("Boy", "Girl"),
                             labels = c("Boy ", "Girl")) #again, stupid way to make sure spacing works out, but do not have time to do it better right now

## Producing strata plots ##
produceStratPlot(d1=mod_csex_cont[mod_csex_cont$level  == "Boy ",], 
                 d2=intxtests_csex_cont, d3=mod_csex_cont, name="strat_csex_3yr", lims=c(-30, 35), ast=20, cuts=c(-20,0,20))




## >>Poverty in 3.5 Year Data --------------------------------------------------------------------------

## Running Interaction Tests ##
meas_pov = runMEIntxCont(outcome="logmev", exposures=exp_list, 
                          adjustment=adjset_meas, modifier="poverty_food_i", df=dr_meas_3yr)
tet_pov = runMEIntxCont(outcome="logtet", exposures=exp_list, 
                         adjustment=adjset_tet, modifier="poverty_food_i", df=dr_tet_3yr)
hib_pov = runMEIntxCont(outcome="loghib", exposures=exp_list, 
                         adjustment=adjset_tet, modifier="poverty_food_i", df=dr_tet_3yr)

intxtests_pov_cont = rbind(meas_pov, tet_pov, hib_pov)
rm(meas_pov, tet_pov, hib_pov)

intxtests_pov_cont$dv = factor(intxtests_pov_cont$dv,
                                levels = c("logmev", "logtet", "loghib"),
                                labels = c("Measles", "Tetanus", "HIB"))

## Getting Stratified Estimates ##
meas_pov = getMEStratEstsCont(outcome="logmev", exposures=exp_list, adjustment=adjset_meas, 
                               stratby="poverty_food_i", l1="Above poverty level", l2="At or below poverty level", df=dr_meas_3yr)

tet_pov = getMEStratEstsCont(outcome="logtet", exposures=exp_list, adjustment=adjset_tet, 
                              stratby="poverty_food_i", l1="Above poverty level", l2="At or below poverty level",  df=dr_tet_3yr)

hib_pov = getMEStratEstsCont(outcome="loghib", exposures=exp_list, adjustment=adjset_tet, 
                              stratby="poverty_food_i", l1="Above poverty level", l2="At or below poverty level",  df=dr_tet_3yr)

mod_pov_cont = rbind(meas_pov, tet_pov, hib_pov)

mod_pov_cont$dv = factor(mod_pov_cont$dv,
                          levels = c("logmev", "logtet", "loghib"),
                          labels = c("Measles", "Tetanus", "HIB"))

mod_pov_cont$level = factor(mod_pov_cont$level, 
                             levels = c("Above poverty level", "At or below poverty level"),
                             labels = c("Above poverty level", "At or below poverty level"))#again, stupid way to make sure spacing works out, but do not have time to do it better right now

## Producing strata plots ##
produceStratPlot(d1=mod_pov_cont[mod_pov_cont$level  == "Above poverty level",], 
                 d2=intxtests_pov_cont, d3=mod_pov_cont, name="strat_pov_3yr", lims=c(-30, 35), ast=20, cuts=c(-20,0,20))





## >>Child Sex in 5 Year Data --------------------------------------------------------------------------

## Running Interaction Tests ##
meas_csex = runMEIntxCont(outcome="logmev", exposures=exp_list, adjustment=adjset_meas[!(adjset_meas=="csex")], modifier="csex", df=dr_meas_5yr)
tet_csex = runMEIntxCont(outcome="logtet", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="csex")], modifier="csex", df=dr_tet_5yr)
hib_csex = runMEIntxCont(outcome="loghib", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="csex")], modifier="csex", df=dr_tet_5yr)

intxtests_csex_cont = rbind(meas_csex, tet_csex, hib_csex)
rm(meas_csex, tet_csex, hib_csex)

intxtests_csex_cont$dv = factor(intxtests_csex_cont$dv,
                                levels = c("logmev", "logtet", "loghib"),
                                labels = c("Measles", "Tetanus", "HIB"))

## Getting Stratified Estimates ##
meas_csex = getMEStratEstsCont(outcome="logmev", exposures=exp_list, adjustment=adjset_meas[!(adjset_meas=="csex")], 
                               stratby="csex", l1="Boy", l2="Girl", df=dr_meas_5yr)

tet_csex = getMEStratEstsCont(outcome="logtet", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="csex")], 
                              stratby="csex", l1="Boy", l2="Girl",  df=dr_tet_5yr)

hib_csex = getMEStratEstsCont(outcome="loghib", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="csex")], 
                              stratby="csex", l1="Boy", l2="Girl",  df=dr_tet_5yr)

mod_csex_cont = rbind(meas_csex, tet_csex, hib_csex)

mod_csex_cont$dv = factor(mod_csex_cont$dv,
                          levels = c("logmev", "logtet", "loghib"),
                          labels = c("Measles", "Tetanus", "HIB"))

mod_csex_cont$level = factor(mod_csex_cont$level, 
                             levels = c("Boy", "Girl"),
                             labels = c("Boy ", "Girl")) #again, stupid way to make sure spacing works out, but do not have time to do it better right now

## Producing strata plots ##
produceStratPlot(d1=mod_csex_cont[mod_csex_cont$level  == "Boy ",], 
                 d2=intxtests_csex_cont, d3=mod_csex_cont, name="strat_csex_5yr", lims=c(-30, 54), ast=20, cuts=c(-20,0,20,40))




## >>Poverty in 5 Year Data --------------------------------------------------------------------------

## Running Interaction Tests ##
meas_pov = runMEIntxCont(outcome="logmev", exposures=exp_list, 
                         adjustment=adjset_meas, modifier="poverty_food_i", df=dr_meas_5yr)
tet_pov = runMEIntxCont(outcome="logtet", exposures=exp_list, 
                        adjustment=adjset_tet, modifier="poverty_food_i", df=dr_tet_5yr)
hib_pov = runMEIntxCont(outcome="loghib", exposures=exp_list, 
                        adjustment=adjset_tet, modifier="poverty_food_i", df=dr_tet_5yr)

intxtests_pov_cont = rbind(meas_pov, tet_pov, hib_pov)
rm(meas_pov, tet_pov, hib_pov)

intxtests_pov_cont$dv = factor(intxtests_pov_cont$dv,
                               levels = c("logmev", "logtet", "loghib"),
                               labels = c("Measles", "Tetanus", "HIB"))

## Getting Stratified Estimates ##
meas_pov = getMEStratEstsCont(outcome="logmev", exposures=exp_list, adjustment=adjset_meas, 
                              stratby="poverty_food_i", l1="Above poverty level", l2="At or below poverty level", df=dr_meas_5yr)

tet_pov = getMEStratEstsCont(outcome="logtet", exposures=exp_list, adjustment=adjset_tet, 
                             stratby="poverty_food_i", l1="Above poverty level", l2="At or below poverty level",  df=dr_tet_5yr)

hib_pov = getMEStratEstsCont(outcome="loghib", exposures=exp_list, adjustment=adjset_tet, 
                             stratby="poverty_food_i", l1="Above poverty level", l2="At or below poverty level",  df=dr_tet_5yr)

mod_pov_cont = rbind(meas_pov, tet_pov, hib_pov)

mod_pov_cont$dv = factor(mod_pov_cont$dv,
                         levels = c("logmev", "logtet", "loghib"),
                         labels = c("Measles", "Tetanus", "HIB"))

mod_pov_cont$level = factor(mod_pov_cont$level, 
                            levels = c("Above poverty level", "At or below poverty level"),
                            labels = c("Above poverty level", "At or below poverty level"))#again, stupid way to make sure spacing works out, but do not have time to do it better right now

## Producing strata plots ##
produceStratPlot(d1=mod_pov_cont[mod_pov_cont$level  == "Above poverty level",], 
                 d2=intxtests_pov_cont, d3=mod_pov_cont, name="strat_pov_5yr", lims=c(-30, 50), ast=20, cuts=c(-20,0,20,40))




## >CATEGORICAL STRATIFIED EFFECTS (EMM) @ EACH TIME POINT ################################################################


## >>Child Sex in 3.5 Year Data ----------------------------------------------------------------------------------


## Running interaction tests ##
meas_csex = testIntx(df=dr_meas_3yr, outcome="notprot_meas", exposures=exp_list, adjustment=adjset_meas, modifier="csex")
tet_csex = testIntx(df=dr_tet_3yr, outcome="notprot_tet", exposures=exp_list, adjustment=adjset_tet, modifier="csex")
hib_csex = testIntx(df=dr_tet_3yr, outcome="notprot_hib", exposures=exp_list, adjustment=adjset_tet, modifier="csex")

intxtests_csex = rbind(meas_csex, tet_csex, hib_csex)
rm(meas_csex, tet_csex, hib_csex)

intxtests_csex$dv = factor(intxtests_csex$dv,
                           levels = c("notprot_meas", "notprot_tet", "notprot_hib"),
                           labels = c("Measles", "Tetanus", "HIB"))


## Getting stratified estimates ##
meas_csex = getStratEstimatesME(df=dr_meas_3yr, outcome="notprot_meas", exposures=exp_list, 
                                adjustment = adjset_meas[adjset_meas != "csex"], l1="Boy", l2="Girl", stratby="csex")

tet_csex = getStratEstimatesME(df=dr_tet_3yr, outcome="notprot_tet", exposures=exp_list, 
                               adjustment = adjset_tet[adjset_tet != "csex"], l1="Boy", l2="Girl", stratby="csex")

hib_csex = getStratEstimatesME(df=dr_tet_3yr, outcome="notprot_hib", exposures=exp_list, 
                               adjustment = adjset_tet[adjset_tet != "csex"], l1="Boy", l2="Girl", stratby="csex")

mod_csex = rbind(meas_csex, tet_csex, hib_csex)

mod_csex$dv = factor(mod_csex$dv,
                     levels = c("notprot_meas", "notprot_tet", "notprot_hib"),
                     labels = c("Measles", "Tetanus", "HIB"))


## Saving Stratified Plot ##
plotStratEstimates(df=mod_csex, sigdf = intxtests_csex, savename="csex_3yr")

rm(meas_csex, tet_csex, hib_csex, intxtests_csex, mod_csex)


## >>Poverty in 3.5 Year Data ---------------------------------------------------------------

## Running interaction tests ##
meas_poverty_food_i = testIntx(df=dr_meas_3yr, outcome="notprot_meas", exposures=exp_list, 
                               adjustment=c(adjset_meas, "poverty_food_i"), modifier="poverty_food_i")
tet_poverty_food_i = testIntx(df=dr_tet_3yr, outcome="notprot_tet", exposures=exp_list, 
                              adjustment=c(adjset_tet,"poverty_food_i"), modifier="poverty_food_i")
hib_poverty_food_i = testIntx(df=dr_tet_3yr, outcome="notprot_hib", exposures=exp_list, 
                              adjustment=c(adjset_tet, "poverty_food_i"), modifier="poverty_food_i")

intxtests_poverty_food_i = rbind(meas_poverty_food_i, tet_poverty_food_i, hib_poverty_food_i)
rm(meas_poverty_food_i, tet_poverty_food_i, hib_poverty_food_i)

intxtests_poverty_food_i$dv = factor(intxtests_poverty_food_i$dv,
                                     levels = c("notprot_meas", "notprot_tet", "notprot_hib"),
                                     labels = c("Measles", "Tetanus", "HIB"))


## Getting stratified estimates ##
meas_foodpov = getStratEstimatesME(df=dr_meas_3yr, outcome="notprot_meas", exposures=exp_list, 
                                   adjustment = adjset_meas, l1="Above poverty level", l2="At or below poverty level", stratby="poverty_food_i")

tet_foodpov = getStratEstimatesME(df=dr_tet_3yr, outcome="notprot_tet", exposures=exp_list, 
                                  adjustment = adjset_tet, l1="Above poverty level", l2="At or below poverty level", stratby="poverty_food_i")

hib_foodpov = getStratEstimatesME(df=dr_tet_3yr, outcome="notprot_hib", exposures=exp_list, 
                                  adjustment = adjset_tet, l1="Above poverty level", l2="At or below poverty level", stratby="poverty_food_i")

mod_foodpov = rbind(meas_foodpov, tet_foodpov, hib_foodpov)

mod_foodpov$dv = factor(mod_foodpov$dv,
                        levels = c("notprot_meas", "notprot_tet", "notprot_hib"),
                        labels = c("Measles", "Tetanus", "HIB"))



## Saving Stratified Plot ##
plotStratEstimates(df=mod_foodpov,sigdf=intxtests_poverty_food_i, savename="foodpov_3yr")

rm(meas_foodpov, tet_foodpov, hib_foodpov, mod_foodpov, intxtests_poverty_food_i)





## >>Child Sex in 5 Year Data ----------------------------------------------------------------------------------


## Running interaction tests ##
meas_csex = testIntx(df=dr_meas_5yr, outcome="notprot_meas", exposures=exp_list, adjustment=adjset_meas, modifier="csex")
tet_csex = testIntx(df=dr_tet_5yr, outcome="notprot_tet", exposures=exp_list, adjustment=adjset_tet, modifier="csex")
hib_csex = testIntx(df=dr_tet_5yr, outcome="notprot_hib", exposures=exp_list, adjustment=adjset_tet, modifier="csex")

intxtests_csex = rbind(meas_csex, tet_csex, hib_csex)
rm(meas_csex, tet_csex, hib_csex)

intxtests_csex$dv = factor(intxtests_csex$dv,
                           levels = c("notprot_meas", "notprot_tet", "notprot_hib"),
                           labels = c("Measles", "Tetanus", "HIB"))


## Getting stratified estimates ##
meas_csex = getStratEstimatesME(df=dr_meas_5yr, outcome="notprot_meas", exposures=exp_list, 
                                adjustment = adjset_meas[adjset_meas != "csex"], l1="Boy", l2="Girl", stratby="csex")

tet_csex = getStratEstimatesME(df=dr_tet_5yr, outcome="notprot_tet", exposures=exp_list, 
                               adjustment = adjset_tet[adjset_tet != "csex"], l1="Boy", l2="Girl", stratby="csex")

hib_csex = getStratEstimatesME(df=dr_tet_5yr, outcome="notprot_hib", exposures=exp_list, 
                               adjustment = adjset_tet[adjset_tet != "csex"], l1="Boy", l2="Girl", stratby="csex")

mod_csex = rbind(meas_csex, tet_csex, hib_csex)

mod_csex$dv = factor(mod_csex$dv,
                     levels = c("notprot_meas", "notprot_tet", "notprot_hib"),
                     labels = c("Measles", "Tetanus", "HIB"))


## Saving Stratified Plot ##
plotStratEstimates(df=mod_csex, sigdf = intxtests_csex, savename="csex_5yr")

rm(meas_csex, tet_csex, hib_csex, intxtests_csex, mod_csex)


## >>Poverty in 5 Year Data ---------------------------------------------------------------

## Running interaction tests ##
meas_poverty_food_i = testIntx(df=dr_meas_5yr, outcome="notprot_meas", exposures=exp_list, adjustment=c(adjset_meas, "poverty_food_i"), modifier="poverty_food_i")
tet_poverty_food_i = testIntx(df=dr_tet_5yr, outcome="notprot_tet", exposures=exp_list, adjustment=c(adjset_tet,"poverty_food_i"), modifier="poverty_food_i")
hib_poverty_food_i = testIntx(df=dr_tet_5yr, outcome="notprot_hib", exposures=exp_list, adjustment=c(adjset_tet, "poverty_food_i"), modifier="poverty_food_i")

intxtests_poverty_food_i = rbind(meas_poverty_food_i, tet_poverty_food_i, hib_poverty_food_i)
rm(meas_poverty_food_i, tet_poverty_food_i, hib_poverty_food_i)

intxtests_poverty_food_i$dv = factor(intxtests_poverty_food_i$dv,
                                     levels = c("notprot_meas", "notprot_tet", "notprot_hib"),
                                     labels = c("Measles", "Tetanus", "HIB"))


## Getting stratified estimates ##
meas_foodpov = getStratEstimatesME(df=dr_meas_5yr, outcome="notprot_meas", exposures=exp_list, 
                                   adjustment = adjset_meas, l1="Above poverty level", l2="At or below poverty level", stratby="poverty_food_i")

tet_foodpov = getStratEstimatesME(df=dr_tet_5yr, outcome="notprot_tet", exposures=exp_list, 
                                  adjustment = adjset_tet, l1="Above poverty level", l2="At or below poverty level", stratby="poverty_food_i")

hib_foodpov = getStratEstimatesME(df=dr_tet_5yr, outcome="notprot_hib", exposures=exp_list, 
                                  adjustment = adjset_tet, l1="Above poverty level", l2="At or below poverty level", stratby="poverty_food_i")

mod_foodpov = rbind(meas_foodpov, tet_foodpov, hib_foodpov)

mod_foodpov$dv = factor(mod_foodpov$dv,
                        levels = c("notprot_meas", "notprot_tet", "notprot_hib"),
                        labels = c("Measles", "Tetanus", "HIB"))



## Saving Stratified Plot ##
plotStratEstimates(df=mod_foodpov,sigdf=intxtests_poverty_food_i, savename="foodpov_5yr")

rm(meas_foodpov, tet_foodpov, hib_foodpov, mod_foodpov, intxtests_poverty_food_i)




## >TETANUS DOSE-RESPONSE @ EACH TIME #################################################################


## >>In 3.5 Year Data ----------------------------------------------------------------


## Getting complete case dataset with datadist ##
f1_cc = dr_tet_3yr[,colnames(dr_tet_3yr) %in% c("log_cisdbca_sg", adjset_tet_lintime, "notprot_tet")]
f1_cc = f1_cc[complete.cases(f1_cc),]
dd = rms::datadist(f1_cc)
options(datadist="dd")


## Poisson regression model ##
f1_m = rms::Glm(data=f1_cc,
                as.formula(paste("notprot_tet", "~", "rms::rcs(log_cisdbca_sg,3)", "+", paste(adjset_tet, collapse="+"))),
                family=poisson)

## Getting predicted values ##
f1_pred = rms::Predict(f1_m, log_cisdbca_sg, conf.type = "simultaneous")

q25 = unname(quantile(f1_cc$log_cisdbca_sg, probs=c(.25)))

f1_contrast = rms::contrast(f1_m, 
                            list(log_cisdbca_sg = f1_pred$log_cisdbca_sg),
                            list(log_cisdbca_sg = DescTools::Closest(f1_pred$log_cisdbca_sg, q25)),
                            conf.type="individual",
                            usebootcoef=T, boot.type="basic") #see contrast.rms in https://cran.r-project.org/web/packages/rms/rms.pdf

f1_pred$or = exp(f1_contrast$Contrast)
f1_pred$lo_or = exp(f1_contrast$Lower)
f1_pred$hi_or = exp(f1_contrast$Upper)


## Plotting Curve ##
tet_cisdbca_plot_3yr = 
  ggplot() +
  geom_ribbon(data=f1_pred, aes(x=log_cisdbca_sg, ymin=lo_or, ymax=hi_or), alpha=0.2, show.legend=F) +
  scale_fill_manual(values=c(
    ghibli_palettes$MarnieMedium2[1:3], ghibli_palettes$PonyoMedium[4], ghibli_palettes$MarnieMedium2[5:7])) + #Change this line for different numbers of cutpoints
  geom_line(data=f1_pred, aes(x=log_cisdbca_sg, y=or), linetype=1, show.legend = F) +
  scale_color_manual(values=c(
    ghibli_palettes$MarnieMedium2[1:3], ghibli_palettes$PonyoMedium[4], ghibli_palettes$MarnieMedium2[5:7])) + #Change this line for different numbers of cutpoints
  scale_x_continuous(expand = c(0, 0)) + #these lil bbs get rid of the padding within the plot!
  scale_y_continuous(trans = 'log', breaks=c(0.1, 1,2,4)) +
  geom_hline(yintercept = 1, linetype="dashed") +
  
  labs(
    y="Relative Risk per 10-fold increase",
    x=expression(paste("log(", italic("cis-"),"DBCA)"))
  ) +
  theme_light() +
  theme(
    text = element_text(family = "Century Gothic"),
    legend.title = element_blank(),
    legend.text = element_text(size=20),
    legend.key.height = unit(3, "cm"),
    legend.key.width = unit(3, "cm"),
    axis.text = element_text(size=20),
    axis.text.x = element_text(size=20),
    axis.title = element_text(size=25),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) 
ggsave(tet_cisdbca_plot_3yr, file="tet_cisdbca_plot_3yr.png", width=7, height=7)


rm(f1_cc, f1_m, f1_pred)


## >>In 5 Year Data ----------------------------------------------------------------


## Getting complete case dataset with datadist ##
f1_cc = dr_tet_5yr[,colnames(dr_tet_5yr) %in% c("log_cisdbca_sg", adjset_tet_lintime, "notprot_tet")]
f1_cc = f1_cc[complete.cases(f1_cc),]
dd = rms::datadist(f1_cc)
options(datadist="dd")


## Poisson regression model ##
f1_m = rms::Glm(data=f1_cc,
                as.formula(paste("notprot_tet", "~", "rms::rcs(log_cisdbca_sg,3)", "+", paste(adjset_tet, collapse="+"))),
                family=poisson)

## Getting predicted values ##
f1_pred = rms::Predict(f1_m, log_cisdbca_sg, conf.type = "simultaneous")

q25 = unname(quantile(f1_cc$log_cisdbca_sg, probs=c(.25)))

f1_contrast = rms::contrast(f1_m, 
                            list(log_cisdbca_sg = f1_pred$log_cisdbca_sg),
                            list(log_cisdbca_sg = DescTools::Closest(f1_pred$log_cisdbca_sg, q25)),
                            conf.type="individual",
                            usebootcoef=T, boot.type="basic") #see contrast.rms in https://cran.r-project.org/web/packages/rms/rms.pdf

f1_pred$or = exp(f1_contrast$Contrast)
f1_pred$lo_or = exp(f1_contrast$Lower)
f1_pred$hi_or = exp(f1_contrast$Upper)


## Plotting Curve ##
tet_cisdbca_plot_5yr = 
  ggplot() +
  geom_ribbon(data=f1_pred, aes(x=log_cisdbca_sg, ymin=lo_or, ymax=hi_or), alpha=0.2, show.legend=F) +
  scale_fill_manual(values=c(
    ghibli_palettes$MarnieMedium2[1:3], ghibli_palettes$PonyoMedium[4], ghibli_palettes$MarnieMedium2[5:7])) + #Change this line for different numbers of cutpoints
  geom_line(data=f1_pred, aes(x=log_cisdbca_sg, y=or), linetype=1, show.legend = F) +
  scale_color_manual(values=c(
    ghibli_palettes$MarnieMedium2[1:3], ghibli_palettes$PonyoMedium[4], ghibli_palettes$MarnieMedium2[5:7])) + #Change this line for different numbers of cutpoints
  scale_x_continuous(expand = c(0, 0)) + #these lil bbs get rid of the padding within the plot!
  scale_y_continuous(trans = 'log', breaks=c(0.1, 1,2,4)) +
  geom_hline(yintercept = 1, linetype="dashed") +
  
  labs(
    y="Relative Risk per 10-fold increase",
    x=expression(paste("log(", italic("cis-"),"DBCA)"))
  ) +
  theme_light() +
  theme(
    text = element_text(family = "Century Gothic"),
    legend.title = element_blank(),
    legend.text = element_text(size=20),
    legend.key.height = unit(3, "cm"),
    legend.key.width = unit(3, "cm"),
    axis.text = element_text(size=20),
    axis.text.x = element_text(size=20),
    axis.title = element_text(size=25),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) 
ggsave(tet_cisdbca_plot_5yr, file="tet_cisdbca_plot_5yr.png", width=7, height=7)





## >GEE MAIN EFFECTS #########################################################################################

## >>Table of Main Effects -----------------------------------------------------------------------

## Measles ##
measles_gee_lin = runGEE(outcome="logmev", exposures = exp_list, adjustment = adjset_meas, df = dr_meas)
measles_gee_bin = runGEEbinary(outcome="notprot_meas", exposures = exp_list, adjustment = adjset_meas, df = dr_meas)


## Tetanus ##
tetanus_gee_lin = runGEE(outcome="logtet", exposures = exp_list, adjustment = adjset_meas, df = dr_tet)
tetanus_gee_bin = runGEEbinary(outcome="notprot_tet", exposures = exp_list, adjustment = adjset_meas, df = dr_tet)


## HIB ##
hib_gee_lin = runGEE(outcome="loghib", exposures = exp_list, adjustment = adjset_meas, df = dr_tet)
hib_gee_bin = runGEEbinary(outcome="notprot_hib", exposures = exp_list, adjustment = adjset_meas, df = dr_tet)


## Combining All and Outputting ##
t3_gee = rbind(cbind(measles_gee_bin, measles_gee_lin),
               cbind(tetanus_gee_bin, tetanus_gee_lin),
               cbind(hib_gee_bin, hib_gee_lin))

write.xlsx(t3_gee, file="analyses_2019-10-30.xlsx",  sheetName="MEsGEE", append=T)




## >>Plotting CONTINUOUS Main Effects --------------------------------------------------------------------------------

meas_me = getGEEOverallEstsCont(outcome="logmev", exposures=exp_list, adjustment=adjset_meas,  df=dr_meas)
tet_me = getGEEOverallEstsCont(outcome="logtet", exposures=exp_list, adjustment=adjset_tet,  df=dr_tet)
hib_me = getGEEOverallEstsCont(outcome="loghib", exposures=exp_list, adjustment=adjset_tet,  df=dr_tet)

mes = rbind(meas_me, tet_me, hib_me)

mes$dv = factor(mes$dv,
                        levels = c("logmev", "logtet", "loghib"),
                        labels = c("Measles", "Tetanus", "Hib"))


produceStratPlot3(d1=mes[mes$dv=="Measles",], d3=mes[mes$dv %in% c("Measles", "Tetanus"),], 
                  d4=mes, name="overall_gee_continuous", cuts=c(-35,80))


## >>Exploring HIB association w PBA ----------------------------------------------------------------------------------

f_hib_pba = as.formula(paste("loghib", "~", "rms::rcs(log_pba_sg,3)", "+", "time", "+", paste(adjset_tet_lintime, collapse="+"))) #Step 2: Create the GEE regression formula 

zmodel = Zelig::zelig(formula=f_hib_pba, model = "normal.gee", id="hsn", data=na.omit(dr_tet), corstr="exchangeable") #using zelig to run GEE so can get the predicted more easily

s1 = createPredictDf(expvect=dr_tet$log_pba_sg, expname="log_pba_sg")

zpred = c(do.call("cbind", Zelig::predict(zmodel, s1)))
s1$estimated = zpred


ggplot() +
  geom_point(data=dr_tet, aes(x=log_pba_sg, y=loghib, colour=as.factor(time))) +
  geom_line(data=s1, aes(x=log_pba_sg, y=estimated))

summary(dr_tet$loghib[dr_tet$time==0])
summary(dr_tet$loghib[dr_tet$time==1])

summary(df$hib)
summary(df$loghib)
summary(yr5$hib_5yr)



## >>Exploratory Splines for MEs -----------------------------------------------------------------------------------


## Tetanus outcome ##
preddf = createPredictDf(expvect=dr_tet$ln_ppddt_i, expname="ln_ppddt_i")
f = logtet ~rms::rcs(ln_ppddt_i,3) + time + mage + income_pers + parity + hiv_combined_i + csex + time_dtp4
zmodel = Zelig::zelig(f, model = "normal.gee", id="hsn", data=na.omit(dr_tet), corstr="exchangeable") #using zelig to run GEE so can get the predicted more easily
zpred = data.frame(pred = do.call("cbind",Zelig::predict(zmodel, preddf)), exp=preddf$ln_ppddt_i) #note for above line, should eventually not use na.omit so broadly - just doing for exploration
ggsave(ggplot() + geom_line(data=zpred, aes(exp, pred)), file="lnppddt-logtet.png", height=4, width=4)

preddf = createPredictDf(expvect=dr_tet$ln_ppdde_i, expname="ln_ppdde_i")
f = logtet ~rms::rcs(ln_ppdde_i,3) + time + mage + income_pers + parity + hiv_combined_i + csex + time_dtp4
zmodel = Zelig::zelig(f, model = "normal.gee", id="hsn", data=na.omit(dr_tet), corstr="exchangeable") #using zelig to run GEE so can get the predicted more easily
zpred = data.frame(pred = do.call("cbind",Zelig::predict(zmodel, preddf)), exp=preddf$ln_ppdde_i) #note for above line, should eventually not use na.omit so broadly - just doing for exploration
ggsave(ggplot() + geom_line(data=zpred, aes(exp, pred)), file="lnppdde-logtet.png", height=4, width=4)

preddf = createPredictDf(expvect=dr_tet$log_cisdbca_sg, expname="log_cisdbca_sg")
f = logtet ~rms::rcs(log_cisdbca_sg,3) + time + mage + income_pers + parity + hiv_combined_i + csex + time_dtp4
zmodel = Zelig::zelig(f, model = "normal.gee", id="hsn", data=na.omit(dr_tet), corstr="exchangeable") #using zelig to run GEE so can get the predicted more easily
zpred = data.frame(pred = do.call("cbind",Zelig::predict(zmodel, preddf)), exp=preddf$log_cisdbca_sg) #note for above line, should eventually not use na.omit so broadly - just doing for exploration
ggsave(ggplot() + geom_line(data=zpred, aes(exp, pred)), file="logcisdbca-logtet.png", height=4, width=4)

preddf = createPredictDf(expvect=dr_tet$log_cisdcca_sg, expname="log_cisdcca_sg")
f = logtet ~rms::rcs(log_cisdcca_sg,3) + time + mage + income_pers + parity + hiv_combined_i + csex + time_dtp4
zmodel = Zelig::zelig(f, model = "normal.gee", id="hsn", data=na.omit(dr_tet), corstr="exchangeable") #using zelig to run GEE so can get the predicted more easily
zpred = data.frame(pred = do.call("cbind",Zelig::predict(zmodel, preddf)), exp=preddf$log_cisdcca_sg) #note for above line, should eventually not use na.omit so broadly - just doing for exploration
ggsave(ggplot() + geom_line(data=zpred, aes(exp, pred)), file="logcisdcca-logtet.png", height=4, width=4)

preddf = createPredictDf(expvect=dr_tet$log_transdcca_sg, expname="log_transdcca_sg")
f = logtet ~rms::rcs(log_transdcca_sg,3) + time + mage + income_pers + parity + hiv_combined_i + csex + time_dtp4
zmodel = Zelig::zelig(f, model = "normal.gee", id="hsn", data=na.omit(dr_tet), corstr="exchangeable") #using zelig to run GEE so can get the predicted more easily
zpred = data.frame(pred = do.call("cbind",Zelig::predict(zmodel, preddf)), exp=preddf$log_transdcca_sg) #note for above line, should eventually not use na.omit so broadly - just doing for exploration
ggsave(ggplot() + geom_line(data=zpred, aes(exp, pred)), file="logtransdcca-logtet.png", height=4, width=4)

preddf = createPredictDf(expvect=dr_tet$log_pba_sg, expname="log_pba_sg")
f = logtet ~rms::rcs(log_pba_sg,3) + time + mage + income_pers + parity + hiv_combined_i + csex + time_dtp4
zmodel = Zelig::zelig(f, model = "normal.gee", id="hsn", data=na.omit(dr_tet), corstr="exchangeable") #using zelig to run GEE so can get the predicted more easily
zpred = data.frame(pred = do.call("cbind",Zelig::predict(zmodel, preddf)), exp=preddf$log_pba_sg) #note for above line, should eventually not use na.omit so broadly - just doing for exploration
ggsave(ggplot() + geom_line(data=zpred, aes(exp, pred)), file="logpba-logtet.png", height=4, width=4)


## Hib outcome ##
preddf = createPredictDf(expvect=dr_tet$ln_ppddt_i, expname="ln_ppddt_i")
f = loghib ~rms::rcs(ln_ppddt_i,3) + time + mage + income_pers + parity + hiv_combined_i + csex + time_dtp4
zmodel = Zelig::zelig(f, model = "normal.gee", id="hsn", data=na.omit(dr_tet), corstr="exchangeable") #using zelig to run GEE so can get the predicted more easily
zpred = data.frame(pred = do.call("cbind",Zelig::predict(zmodel, preddf)), exp=preddf$ln_ppddt_i) #note for above line, should eventually not use na.omit so broadly - just doing for exploration
ggsave(ggplot() + geom_line(data=zpred, aes(exp, pred)), file="lnppddt-loghib.png", height=4, width=4)

preddf = createPredictDf(expvect=dr_tet$ln_ppdde_i, expname="ln_ppdde_i")
f = loghib ~rms::rcs(ln_ppdde_i,3) + time + mage + income_pers + parity + hiv_combined_i + csex + time_dtp4
zmodel = Zelig::zelig(f, model = "normal.gee", id="hsn", data=na.omit(dr_tet), corstr="exchangeable") #using zelig to run GEE so can get the predicted more easily
zpred = data.frame(pred = do.call("cbind",Zelig::predict(zmodel, preddf)), exp=preddf$ln_ppdde_i) #note for above line, should eventually not use na.omit so broadly - just doing for exploration
ggsave(ggplot() + geom_line(data=zpred, aes(exp, pred)), file="lnppdde-loghib.png", height=4, width=4)

preddf = createPredictDf(expvect=dr_tet$log_cisdbca_sg, expname="log_cisdbca_sg")
f = loghib ~rms::rcs(log_cisdbca_sg,3) + time + mage + income_pers + parity + hiv_combined_i + csex + time_dtp4
zmodel = Zelig::zelig(f, model = "normal.gee", id="hsn", data=na.omit(dr_tet), corstr="exchangeable") #using zelig to run GEE so can get the predicted more easily
zpred = data.frame(pred = do.call("cbind",Zelig::predict(zmodel, preddf)), exp=preddf$log_cisdbca_sg) #note for above line, should eventually not use na.omit so broadly - just doing for exploration
ggsave(ggplot() + geom_line(data=zpred, aes(exp, pred)), file="logcisdbca-loghib.png", height=4, width=4)

preddf = createPredictDf(expvect=dr_tet$log_cisdcca_sg, expname="log_cisdcca_sg")
f = loghib ~rms::rcs(log_cisdcca_sg,3) + time + mage + income_pers + parity + hiv_combined_i + csex + time_dtp4
zmodel = Zelig::zelig(f, model = "normal.gee", id="hsn", data=na.omit(dr_tet), corstr="exchangeable") #using zelig to run GEE so can get the predicted more easily
zpred = data.frame(pred = do.call("cbind",Zelig::predict(zmodel, preddf)), exp=preddf$log_cisdcca_sg) #note for above line, should eventually not use na.omit so broadly - just doing for exploration
ggsave(ggplot() + geom_line(data=zpred, aes(exp, pred)), file="logcisdcca-loghib.png", height=4, width=4)

preddf = createPredictDf(expvect=dr_tet$log_transdcca_sg, expname="log_transdcca_sg")
f = loghib ~rms::rcs(log_transdcca_sg,3) + time + mage + income_pers + parity + hiv_combined_i + csex + time_dtp4
zmodel = Zelig::zelig(f, model = "normal.gee", id="hsn", data=na.omit(dr_tet), corstr="exchangeable") #using zelig to run GEE so can get the predicted more easily
zpred = data.frame(pred = do.call("cbind",Zelig::predict(zmodel, preddf)), exp=preddf$log_transdcca_sg) #note for above line, should eventually not use na.omit so broadly - just doing for exploration
ggsave(ggplot() + geom_line(data=zpred, aes(exp, pred)), file="logtransdcca-loghib.png", height=4, width=4)

preddf = createPredictDf(expvect=dr_tet$log_pba_sg, expname="log_pba_sg")
f = loghib ~rms::rcs(log_pba_sg,3) + time + mage + income_pers + parity + hiv_combined_i + csex + time_dtp4
zmodel = Zelig::zelig(f, model = "normal.gee", id="hsn", data=na.omit(dr_tet), corstr="exchangeable") #using zelig to run GEE so can get the predicted more easily
zpred = data.frame(pred = do.call("cbind",Zelig::predict(zmodel, preddf)), exp=preddf$log_pba_sg) #note for above line, should eventually not use na.omit so broadly - just doing for exploration
ggsave(ggplot() + geom_line(data=zpred, aes(exp, pred)), file="logpba-loghib.png", height=4, width=4)



## >PLOTTING [SPLINE] MAIN EFFECTS ####################################################################################################


## >>Cis DBCA --> Tetanus, Linear -----------------------------------------------------------------

dr_tet_model = dr_tet[,colnames(dr_tet) %in% 
                            c("notprot_tet", "log_cisdbca_sg", "time", "mage",
                              "income_pers", "parity", "hiv_combined_i", "csex", "time_dtp4")] #creating dataframe with only model-relevant variables so na.omit will not kick everything out 

s1 = createPredictDf(expvect=dr_tet$log_cisdbca_sg, expname="log_cisdbca_sg") #Step 1: Create a prediction dataframe (with values to predict outcome from model)

s2 =as.formula(paste("logtet", "~", "log_cisdbca_sg", "+", "time", "+", paste(adjset_tet_lintime, collapse="+"))) #Step 2: Create the GEE regression formula 

s3 = boot::boot(data=na.omit(dr_tet_model), statistic=runGEEAndPredict, R=10000, f=s2, preddf=s1) #Step 3: Run GEE and boot the predicted values (simulate R sets of predicted values)

s4 = bootCIsGEE(results=s3, range_exp=s1$log_cisdbca_sg) #Boot confidence intervals



## Plotting the Linear Effect ##
tet_cisdbca_lin_plot = 
  ggplot() +
  geom_ribbon(data=s4, aes(x=exp, ymin=lcl, ymax=ucl), alpha=0.2, show.legend=F, fill=c(ghibli_palettes$SpiritedMedium[3])) +
  #scale_fill_manual(values=c(ghibli_palettes$SpiritedMedium[3])) + #Change this line for different numbers of cutpoints
  geom_line(data=s4, aes(x=exp, y=estq), linetype=1, show.legend = F, colour = c(ghibli_palettes$SpiritedMedium[3])) +
  #scale_color_manual(values=c(ghibli_palettes$SpiritedMedium[3])) + #Change this line for different numbers of cutpoints
  scale_x_continuous(expand = c(0, 0)) + #these lil bbs get rid of the padding within the plot!
  #scale_y_continuous(trans = 'log', breaks=c(0.1, 1,2,4)) +
  #geom_hline(yintercept = 1, linetype="dashed") +
  
  labs(
    y=expression("Change in log(tetanus antibodies)"),
    x=expression(paste("log(", italic("cis-"),"DBCA)"))
  ) +
  theme_light() +
  theme(
    text = element_text(family = "Century Gothic"),
    legend.title = element_blank(),
    legend.text = element_text(size=20),
    legend.key.height = unit(3, "cm"),
    legend.key.width = unit(3, "cm"),
    axis.text = element_text(size=20),
    axis.text.x = element_text(size=20),
    axis.title = element_text(size=25),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) 
ggsave(tet_cisdbca_lin_plot, file="tet_cisdbca_lin_plot_gee.png", width=7, height=7)

rm (s1, s2, s3, s4, dr_tet_model)


## >>Cis DBCA --> Tetanus, Spline on Exposure -----------------------------------------------------------------

dr_tet_model = dr_tet[,colnames(dr_tet) %in% 
                        c("hsn" ,"logtet", "log_cisdbca_sg", "time", "mage",
                          "income_pers", "parity", "hiv_combined_i", "csex", "time_dtp4")] #creating dataframe with only model-relevant variables so na.omit will not kick everything out 


s1 = createPredictDf(expvect=dr_tet$log_cisdbca_sg, expname="log_cisdbca_sg") #Step 1: Create a prediction dataframe (with values to predict outcome from model)

s2 = as.formula(paste("logtet", "~", "rms::rcs(log_cisdbca_sg,3)", "+", "time", "+", paste(adjset_tet_lintime, collapse="+"))) #Step 2: Create the GEE regression formula 

s3 = boot::boot(data=na.omit(dr_tet_model), statistic=runGEEAndPredict, R=5000, f=s2, preddf=s1) #Step 3: Run GEE and boot the predicted values (simulate R sets of predicted values)

s4 = bootCIsGEE(results=s3, range_exp=s1$log_cisdbca_sg) #Boot confidence intervals


s4$exp = 10^(s4$exp) #back-transforming to original values
s4$lcl = 10^(s4$lcl)
s4$ucl = 10^(s4$ucl)
s4$estq = 10^(s4$estq)



## Plotting the Linear Effect ##
tet_cisdbca_spline_plot = 
  ggplot() +
  geom_ribbon(data=s4, aes(x=exp, ymin=lcl, ymax=ucl), alpha=0.2, show.legend=F) +
  #scale_fill_manual(values=c(ghibli_palettes$SpiritedMedium[3])) + #Change this line for different numbers of cutpoints
  geom_line(data=s4, aes(x=exp, y=estq), linetype=1, show.legend = F) +
  geom_rug(data=dr_tet, aes(10^log_cisdbca_sg), alpha=0.5) +
  #scale_color_manual(values=c(ghibli_palettes$SpiritedMedium[3])) + #Change this line for different numbers of cutpoints
  
  scale_x_continuous(trans= 'log10',expand = c(0, 0)) + #these lil bbs get rid of the padding within the plot!
  scale_y_continuous(trans = 'log10', expand=c(0,0)) +
  
  #geom_hline(yintercept = 1, linetype="dashed") +
  
  labs(
    y=expression("Tetanus Antibodies (IU/mL)"),
    x=expression(paste(italic("cis-"),"DBCA"))
  ) +
  theme_light() +
  theme(
    text = element_text(family = "Century Gothic"),
    legend.title = element_blank(),
    legend.text = element_text(size=20),
    legend.key.height = unit(3, "cm"),
    legend.key.width = unit(3, "cm"),
    axis.text = element_text(size=20),
    axis.text.x = element_text(size=20),
    axis.title = element_text(size=25),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) 
ggsave(tet_cisdbca_spline_plot, file="tet_cisdbca_spline_plot_gee_transformed.png", width=7, height=7)

rm (s1, s2, s3, s4)




## >>Cis DCCA --> Tetanus, Spline on Exposure -----------------------------------------------------------------


s1 = createPredictDf(expvect=dr_tet$log_cisdcca_sg, expname="log_cisdcca_sg") #Step 1: Create a prediction dataframe (with values to predict outcome from model)

s2 = as.formula(paste("logtet", "~", "rms::rcs(log_cisdcca_sg,3)", "+", "time", "+", paste(adjset_tet_lintime, collapse="+"))) #Step 2: Create the GEE regression formula 

s3 = boot::boot(data=na.omit(dr_tet), statistic=runGEEAndPredict, R=5000, f=s2, preddf=s1) #Step 3: Run GEE and boot the predicted values (simulate R sets of predicted values)

s4 = bootCIsGEE(results=s3, range_exp=s1$log_cisdcca_sg) #Boot confidence intervals


s4$exp = 10^(s4$exp) #back-transforming to original values
s4$lcl = 10^(s4$lcl)
s4$ucl = 10^(s4$ucl)
s4$estq = 10^(s4$estq)



## Plotting the Linear Effect ##
tet_cisdcca_spline_plot = 
  ggplot() +
  geom_ribbon(data=s4, aes(x=exp, ymin=lcl, ymax=ucl), alpha=0.2, show.legend=F, fill=c(ghibli_palettes$SpiritedMedium[3])) +
  #scale_fill_manual(values=c(ghibli_palettes$SpiritedMedium[3])) + #Change this line for different numbers of cutpoints
  geom_line(data=s4, aes(x=exp, y=estq), linetype=1, show.legend = F, colour = c(ghibli_palettes$SpiritedMedium[3])) +
  geom_rug(data=dr_tet, aes(10^log_cisdcca_sg), alpha=0.5) +
  #scale_color_manual(values=c(ghibli_palettes$SpiritedMedium[3])) + #Change this line for different numbers of cutpoints
  
  scale_x_continuous(trans= 'log10',expand = c(0, 0),
                     limits = c(10^(unname(quantile(dr_tet$log_cisdcca_sg, probs = c(0.01, 0.99), na.rm=T))))) + #these lil bbs get rid of the padding within the plot!
  scale_y_continuous(trans = 'log10', expand=c(0,0)) +
  
  #geom_hline(yintercept = 1, linetype="dashed") +
  
  labs(
    y=expression("Tetanus Antibodies (IU/mL)"),
    x=expression(paste(italic("cis-"),"DCCA"))
  ) +
  theme_light() +
  theme(
    text = element_text(family = "Century Gothic"),
    legend.title = element_blank(),
    legend.text = element_text(size=20),
    legend.key.height = unit(3, "cm"),
    legend.key.width = unit(3, "cm"),
    axis.text = element_text(size=20),
    axis.text.x = element_text(size=20),
    axis.title = element_text(size=25),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) 
ggsave(tet_cisdcca_spline_plot, file="tet_cisdcca_spline_plot_gee_transformed.png", width=7, height=7)

rm (s1, s2, s3, s4)



## >>Trans DCCA --> Tetanus, Spline on Exposure -----------------------------------------------------------------

ggplot()+
  geom_freqpoly(data=d, aes(tet))

s1 = createPredictDf(expvect=dr_tet$log_transdcca_sg, expname="log_transdcca_sg") #Step 1: Create a prediction dataframe (with values to predict outcome from model)

s2 = as.formula(paste("logtet", "~", "rms::rcs(log_transdcca_sg,3)", "+", "time", "+", paste(adjset_tet_lintime, collapse="+"))) #Step 2: Create the GEE regression formula 

s3 = boot::boot(data=na.omit(dr_tet), statistic=runGEEAndPredict, R=5000, f=s2, preddf=s1) #Step 3: Run GEE and boot the predicted values (simulate R sets of predicted values)

s4 = bootCIsGEE(results=s3, range_exp=s1$log_transdcca_sg) #Boot confidence intervals


s4$exp = 10^(s4$exp) #back-transforming to original values
s4$lcl = 10^(s4$lcl)
s4$ucl = 10^(s4$ucl)
s4$estq = 10^(s4$estq)



## Plotting the Linear Effect ##
tet_transdcca_spline_plot = 
  ggplot() +
  geom_ribbon(data=s4, aes(x=exp, ymin=lcl, ymax=ucl), alpha=0.2, show.legend=F, fill=c(ghibli_palettes$SpiritedMedium[3])) +
  #scale_fill_manual(values=c(ghibli_palettes$SpiritedMedium[3])) + #Change this line for different numbers of cutpoints
  geom_line(data=s4, aes(x=exp, y=estq), linetype=1, show.legend = F, colour = c(ghibli_palettes$SpiritedMedium[3])) +
  geom_rug(data=dr_tet, aes(10^log_transdcca_sg), alpha=0.5) +
  #scale_color_manual(values=c(ghibli_palettes$SpiritedMedium[3])) + #Change this line for different numbers of cutpoints
  
  scale_x_continuous(trans= 'log10',expand = c(0, 0),
                     limits = c(10^(unname(quantile(dr_tet$log_transdcca_sg, probs = c(0.01, 0.99), na.rm=T))))) + #these lil bbs get rid of the padding within the plot!
  scale_y_continuous(trans = 'log10', expand=c(0,0)) +
  
  #geom_hline(yintercept = 1, linetype="dashed") +
  
  labs(
    y=expression("Tetanus Antibodies (IU/mL)"),
    x=expression(paste(italic("trans-"),"DCCA"))
  ) +
  theme_light() +
  theme(
    text = element_text(family = "Century Gothic"),
    legend.title = element_blank(),
    legend.text = element_text(size=20),
    legend.key.height = unit(3, "cm"),
    legend.key.width = unit(3, "cm"),
    axis.text = element_text(size=20),
    axis.text.x = element_text(size=20),
    axis.title = element_text(size=25),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) 
ggsave(tet_transdcca_spline_plot, file="tet_transdcca_spline_plot_gee_transformed.png", width=7, height=7)

rm (s1, s2, s3, s4)




## >>3-PBA --> Tetanus, Spline on Exposure -----------------------------------------------------------------


s1 = createPredictDf(expvect=dr_tet$log_pba_sg, expname="log_pba_sg") #Step 1: Create a prediction dataframe (with values to predict outcome from model)

s2 = as.formula(paste("logtet", "~", "rms::rcs(log_pba_sg,3)", "+", "time", "+", paste(adjset_tet_lintime, collapse="+"))) #Step 2: Create the GEE regression formula 

s3 = boot::boot(data=na.omit(dr_tet), statistic=runGEEAndPredict, R=5000, f=s2, preddf=s1) #Step 3: Run GEE and boot the predicted values (simulate R sets of predicted values)

s4 = bootCIsGEE(results=s3, range_exp=s1$log_pba_sg) #Boot confidence intervals


s4$exp = 10^(s4$exp) #back-transforming to original values
s4$lcl = 10^(s4$lcl)
s4$ucl = 10^(s4$ucl)
s4$estq = 10^(s4$estq)



## Plotting the Linear Effect ##
tet_3pba_spline_plot = 
  ggplot() +
  geom_ribbon(data=s4, aes(x=exp, ymin=lcl, ymax=ucl), alpha=0.2, show.legend=F, fill=c(ghibli_palettes$SpiritedMedium[3])) +
  #scale_fill_manual(values=c(ghibli_palettes$SpiritedMedium[3])) + #Change this line for different numbers of cutpoints
  geom_line(data=s4, aes(x=exp, y=estq), linetype=1, show.legend = F, colour = c(ghibli_palettes$SpiritedMedium[3])) +
  geom_rug(data=dr_tet, aes(10^log_pba_sg), alpha=0.5) +
  #scale_color_manual(values=c(ghibli_palettes$SpiritedMedium[3])) + #Change this line for different numbers of cutpoints
  
  scale_x_continuous(trans= 'log10',expand = c(0, 0),
                     limits = c(10^(unname(quantile(dr_tet$log_pba_sg, probs = c(0.01, 0.99), na.rm=T))))) + #these lil bbs get rid of the padding within the plot!
  scale_y_continuous(trans = 'log10', expand=c(0,0)) +
  
  #geom_hline(yintercept = 1, linetype="dashed") +
  
  labs(
    y=expression("Tetanus Antibodies (IU/mL)"),
    x="3-PBA"
  ) +
  theme_light() +
  theme(
    text = element_text(family = "Avenir Next Cyr W04 Regular"),
    legend.title = element_blank(),
    legend.text = element_text(size=20),
    legend.key.height = unit(3, "cm"),
    legend.key.width = unit(3, "cm"),
    axis.text = element_text(size=20),
    axis.text.x = element_text(size=20),
    axis.title = element_text(size=25),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) 
ggsave(tet_3pba_spline_plot, file="tet_3pba_spline_plot_gee_transformed.png", width=7, height=7)

rm (s1, s2, s3, s4)






## >>Cis DBCA --> Not Protected Tetanus, Spline on Exposure, 25th ref-----------------------------------------------------------------


dr_tet_model = dr_tet[,colnames(dr_tet) %in% 
                        c("hsn" ,"notprot_tet", "log_cisdbca_sg", "time", "mage",
                          "income_pers", "parity", "hiv_combined_i", "csex", "time_dtp4")] #creating dataframe with only model-relevant variables so na.omit will not kick everything out 

test=na.omit(dr_tet_model)

s1 = createPredictDf(expvect=dr_tet$log_cisdbca_sg, expname="log_cisdbca_sg") #Step 1: Create a prediction dataframe (with values to predict outcome from model)

s2 = as.formula(paste("notprot_tet", "~", "rms::rcs(log_cisdbca_sg,3)", "+", "time", "+", paste(adjset_tet_lintime, collapse="+"))) #Step 2: Create the GEE regression formula 

s3 = boot::boot(data=na.omit(dr_tet_model), statistic=runGEEAndPredictBinaryLogit, R=10000, f=s2, preddf=s1, q25=unname(quantile(dr_tet_model$log_cisdbca_sg, probs=c(.25), na.rm=T))) #Step 3: Run GEE and boot the predicted values (simulate R sets of predicted values)


#s3 = boot::boot(data=na.omit(dr_tet), statistic=runGEEAndPredictBinaryLogit, R=10000, f=s2, preddf=s1, q25=unname(quantile(dr_tet$log_cisdbca_sg, probs=c(.25), na.rm=T))) #Step 3: Run GEE and boot the predicted values (simulate R sets of predicted values)

s4 = bootCIsGEEBinary(results=s3, range_exp=s1$log_cisdbca_sg, q25=unname(quantile(dr_tet$log_cisdbca_sg, probs=c(.25), na.rm=T))) #Boot confidence intervals


s4$exp = 10^(s4$exp) #back-transforming to original values



summary(dr_tet$log_cisdbca_sg)
unname(quantile(dr_tet$log_cisdbca_sg, probs = c(0.025, 0.975), na.rm=T))

## Plotting the Linear Effect ##
notprottet_cisdbca_spline_plot = 
  ggplot() +
  geom_ribbon(data=s4, aes(x=exp, ymin=lcl, ymax=ucl), alpha=0.2, show.legend=F) +
  #scale_fill_manual(values=c(ghibli_palettes$SpiritedMedium[3])) + #Change this line for different numbers of cutpoints
  geom_line(data=s4, aes(x=exp, y=estq), linetype=1, show.legend = F) +
  #scale_color_manual(values=c(ghibli_palettes$SpiritedMedium[3])) + #Change this line for different numbers of cutpoints
  geom_rug(data=dr_tet, aes(10^log_cisdbca_sg), alpha=0.5) +
  scale_x_continuous(expand = c(0, 0),trans='log10',
                     limits = c(10^(unname(quantile(dr_tet$log_cisdbca_sg, probs = c(0.025, 0.975), na.rm=T))))) + #these lil bbs get rid of the padding within the plot!
  scale_y_continuous(trans = 'log', breaks=c(1,2), limits = c(0.5,2.5)) +
  geom_hline(yintercept = 1, linetype="dashed") +
  
  labs(
    y=expression("Relative Risk"),
    x=expression(paste(italic("cis-"),"DBCA"))
  ) +
  theme_light() +
  theme(
    text = element_text(family = "Century Gothic"),
    legend.title = element_blank(),
    legend.text = element_text(size=20),
    legend.key.height = unit(3, "cm"),
    legend.key.width = unit(3, "cm"),
    axis.text = element_text(size=20),
    axis.text.x = element_text(size=20),
    axis.title = element_text(size=25),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) 
ggsave(notprottet_cisdbca_spline_plot, file="notprottet_cisdbca_spline_plot_gee-transformed.png", width=7, height=7)



## >CATEGORICAL INTERACTION EFFECTS ##############################################################################################



## >>Modifier: Food Poverty ------------------------------------------------------------------------------------


## Running Interaction Tests ##
meas_pov = runGEEIntxBinary(outcome="notprot_meas", exposures=exp_list, adjustment=adjset_meas, modifier="poverty_food_i", df=dr_meas)
tet_pov = runGEEIntxBinary(outcome="notprot_tet", exposures=exp_list, adjustment=adjset_tet, modifier="poverty_food_i", df=dr_tet)
hib_pov = runGEEIntxBinary(outcome="notprot_hib", exposures=exp_list, adjustment=adjset_tet, modifier="poverty_food_i", df=dr_tet)

intxtests_poverty_food_i = rbind(meas_pov, tet_pov, hib_pov)
rm(meas_pov, tet_pov, hib_pov)

intxtests_poverty_food_i$dv = factor(intxtests_poverty_food_i$dv,
                                     levels = c("notprot_meas", "notprot_tet", "notprot_hib"),
                                     labels = c("Measles", "Tetanus", "HIB"))



## Getting Stratified Estimates ##
meas_foodpov = getGEEStratEsts(outcome="notprot_meas", exposures=exp_list, adjustment=adjset_meas, 
                               stratby="poverty_food_i", l1="Above poverty level", l2="At or below poverty level", df=dr_meas)

tet_foodpov = getGEEStratEsts(outcome="notprot_tet", exposures=exp_list, adjustment=adjset_tet, 
                               stratby="poverty_food_i", l1="Above poverty level", l2="At or below poverty level", df=dr_tet)

hib_foodpov = getGEEStratEsts(outcome="notprot_hib", exposures=exp_list, adjustment=adjset_tet, 
                              stratby="poverty_food_i", l1="Above poverty level", l2="At or below poverty level", df=dr_tet)

mod_foodpov = rbind(meas_foodpov, tet_foodpov, hib_foodpov)

mod_foodpov$dv = factor(mod_foodpov$dv,
                        levels = c("notprot_meas", "notprot_tet", "notprot_hib"),
                        labels = c("Measles", "Tetanus", "HIB"))


mod_foodpov$level = factor(mod_foodpov$level, #this is a very stupid fix to keep the graph/legend aligned...want to do a better one when have more time
                           levels = c("Above poverty level", "At or below poverty level"),
                           labels = c("Above poverty level", "At or below\npoverty level"))


## Producing strata plots ##
produceStratPlot(d1=mod_foodpov[mod_foodpov$level=="Above poverty level",], 
                 d2=intxtests_poverty_food_i, d3=mod_foodpov, name="strat_foodpov_gee")



## >>Modifier: Child Sex ------------------------------------------------------------------------------------


## Running Interaction Tests ##
meas_csex = runGEEIntxBinary(outcome="notprot_meas", exposures=exp_list, adjustment=adjset_meas[!(adjset_meas=="csex")], modifier="csex", df=dr_meas)
tet_csex = runGEEIntxBinary(outcome="notprot_tet", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="csex")], modifier="csex", df=dr_tet)
hib_csex = runGEEIntxBinary(outcome="notprot_hib", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="csex")], modifier="csex", df=dr_tet)

intxtests_csex = rbind(meas_csex, tet_csex, hib_csex)
rm(meas_csex, tet_csex, hib_csex)

intxtests_csex$dv = factor(intxtests_csex$dv,
                           levels = c("notprot_meas", "notprot_tet", "notprot_hib"),
                           labels = c("Measles", "Tetanus", "HIB"))



## Getting Stratified Estimates ##
meas_csex = getGEEStratEsts(outcome="notprot_meas", exposures=exp_list, adjustment=adjset_meas[!(adjset_meas=="csex")], 
                               stratby="csex", l1="Boy", l2="Girl", df=dr_meas)

tet_csex = getGEEStratEsts(outcome="notprot_tet", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="csex")], 
                              stratby="csex", l1="Boy", l2="Girl",  df=dr_tet)

hib_csex = getGEEStratEsts(outcome="notprot_hib", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="csex")], 
                              stratby="csex", l1="Boy", l2="Girl",  df=dr_tet)

mod_csex = rbind(meas_csex, tet_csex, hib_csex)

mod_csex$dv = factor(mod_csex$dv,
                     levels = c("notprot_meas", "notprot_tet", "notprot_hib"),
                     labels = c("Measles", "Tetanus", "HIB"))

mod_csex$level = factor(mod_csex$level, 
                        levels = c("Boy", "Girl"),
                        labels = c("Boy ", "Girl"))


## Producing strata plots ##
produceStratPlot(d1=mod_csex[mod_csex$level  == "Boy ",], 
                 d2=intxtests_csex, d3=mod_csex, name="strat_csex_gee")



## >>Modifier: Maternal HIV ------------------------------------------------------------------------------------


## Running Interaction Tests ##
meas_hiv = runGEEIntxBinary(outcome="notprot_meas", exposures=exp_list, adjustment=adjset_meas[!(adjset_meas=="hiv_combined_i")], modifier="hiv_combined_i", df=dr_meas)
tet_hiv = runGEEIntxBinary(outcome="notprot_tet", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="hiv_combined_i")], modifier="hiv_combined_i", df=dr_tet)
hib_hiv = runGEEIntxBinary(outcome="notprot_hib", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="hiv_combined_i")], modifier="hiv_combined_i", df=dr_tet)

intxtests_hiv = rbind(meas_hiv, tet_hiv, hib_hiv)
rm(meas_hiv, tet_hiv, hib_hiv)

intxtests_hiv$dv = factor(intxtests_hiv$dv,
                           levels = c("notprot_meas", "notprot_tet", "notprot_hib"),
                           labels = c("Measles", "Tetanus", "HIB"))



## Getting Stratified Estimates ##
meas_hiv = getGEEStratEsts(outcome="notprot_meas", exposures=exp_list, adjustment=adjset_meas[!(adjset_meas=="hiv_combined_i")], 
                            stratby="hiv_combined_i", l1="HIV Negative", l2="HIV Positive", df=dr_meas)

tet_hiv = getGEEStratEsts(outcome="notprot_tet", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="hiv_combined_i")], 
                           stratby="hiv_combined_i", l1="HIV Negative", l2="HIV Positive",  df=dr_tet)

hib_hiv = getGEEStratEsts(outcome="notprot_hib", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="hiv_combined_i")], 
                           stratby="hiv_combined_i", l1="HIV Negative", l2="HIV Positive",  df=dr_tet)

mod_hiv = rbind(meas_hiv, tet_hiv, hib_hiv)

mod_hiv$dv = factor(mod_hiv$dv,
                     levels = c("notprot_meas", "notprot_tet", "notprot_hib"),
                     labels = c("Measles", "Tetanus", "HIB"))



## Plotting ##
mod_hiv_plot = 
  ggplot() + 
  geom_pointrange(data=mod_hiv,
                  aes(x=predictor, y=rr, ymin=lo_rr, ymax=hi_rr, fill=level), pch=21,
                  position=position_dodge(width=0.6), size=1.5) +
  geom_text(data=intxtests_hiv, aes(x=predictor, y=2, label=siglabel), size=20) + 
  facet_grid(dv~.) + 
  scale_y_continuous(trans="log", breaks=c(0.25, 0.5, 1, 2, 4)) +
  #coord_cartesian(ylim=c(0.25,6)) +
  scale_fill_manual(values=c(ghibli_palettes$MononokeMedium[c(3,5)])) + 
  scale_x_discrete(labels = c(expression(paste(italic("o,p'"), "-DDE")),
                              expression(paste(italic("p,p'"), "-DDE")),
                              expression(paste(italic("p,p'"), "-DDT")),
                              expression(paste(italic("cis"), "-DBCA")),
                              expression(paste(italic("cis"), "-DCCA")),
                              expression(paste(italic("trans"), "-DCCA")),
                              "PBA")) +
  geom_hline(yintercept=1, linetype='dashed') +
  labs(
    y="Relative Risk per 10-fold increase\n",
    x="\nExposure") +
  theme_light() +
  theme(
    text = element_text(family = "Avenir Next Cyr W04 Regular"),
    legend.title = element_blank(),
    legend.text = element_text(size=18),
    legend.key.height = unit(3, "cm"),
    legend.key.width = unit(3, "cm"),
    legend.position = "right",
    axis.text = element_text(size=20),
    axis.text.x = element_text(angle = 45, size=20, hjust=1),
    axis.title = element_text(size=25),
    plot.caption = element_text(size=20, hjust=0),
    strip.background = element_blank(),
    strip.text = element_text(size=20, colour="black")
  ) 
ggsave(mod_hiv_plot, file="strat_by_hiv_gee.png", height=8, width =12)




## >>Modifier: Food Security ------------------------------------------------------------------------------------


## Running Interaction Tests ##
meas_insecure_low = runGEEIntxBinary(outcome="notprot_meas", exposures=exp_list, adjustment=adjset_meas, modifier="insecure_low", df=dr_meas)
tet_insecure_low = runGEEIntxBinary(outcome="notprot_tet", exposures=exp_list, adjustment=adjset_tet, modifier="insecure_low", df=dr_tet)
hib_insecure_low = runGEEIntxBinary(outcome="notprot_hib", exposures=exp_list, adjustment=adjset_tet, modifier="insecure_low", df=dr_tet)

intxtests_insecure_low = rbind(meas_insecure_low, tet_insecure_low, hib_insecure_low)
rm(meas_insecure_low, tet_insecure_low, hib_insecure_low)

intxtests_insecure_low$dv = factor(intxtests_insecure_low$dv,
                                   levels = c("notprot_meas", "notprot_tet", "notprot_hib"),
                                   labels = c("Measles", "Tetanus", "HIB"))



## Getting Stratified Estimates ##
meas_foodsec = getGEEStratEsts(outcome="notprot_meas", exposures=exp_list, adjustment=adjset_meas, 
                           stratby="insecure_low", l1="Food Secure", l2="Low or Very Low Food Security", df=dr_meas)

tet_foodsec = getGEEStratEsts(outcome="notprot_tet", exposures=exp_list, adjustment=adjset_tet, 
                          stratby="insecure_low", l1="Food Secure", l2="Low or Very Low Food Security",  df=dr_tet)

hib_foodsec = getGEEStratEsts(outcome="notprot_hib", exposures=exp_list, adjustment=adjset_tet, 
                          stratby="insecure_low", l1="Food Secure", l2="Low or Very Low Food Security",  df=dr_tet)


mod_foodsec = rbind(meas_foodsec, tet_foodsec, hib_foodsec)

mod_foodsec$dv = factor(mod_foodsec$dv,
                        levels = c("notprot_meas", "notprot_tet", "notprot_hib"),
                        labels = c("Measles", "Tetanus", "HIB"))



## Plotting ##
mod_foodsec_plot = 
  ggplot() + 
  geom_pointrange(data=mod_foodsec,
                  aes(x=predictor, y=rr, ymin=lo_rr, ymax=hi_rr, fill=level), pch=21,
                  position=position_dodge(width=0.6), size=1.5) +
  geom_text(data=intxtests_insecure_low, aes(x=predictor, y=2, label=siglabel), size=20) + 
  facet_grid(dv~.) + 
  scale_y_continuous(trans="log", breaks=c(0.25, 0.5, 1, 2, 4)) +
  #coord_cartesian(ylim=c(0.25,6)) +
  scale_fill_manual(values=c(ghibli_palettes$MononokeMedium[c(3,5)])) + 
  scale_x_discrete(labels = c(expression(paste(italic("o,p'"), "-DDE")),
                              expression(paste(italic("p,p'"), "-DDE")),
                              expression(paste(italic("p,p'"), "-DDT")),
                              expression(paste(italic("cis"), "-DBCA")),
                              expression(paste(italic("cis"), "-DCCA")),
                              expression(paste(italic("trans"), "-DCCA")),
                              "PBA")) +
  geom_hline(yintercept=1, linetype='dashed') +
  labs(
    y="Relative Risk per 10-fold increase\n",
    x="\nExposure") +
  theme_light() +
  theme(
    text = element_text(family = "Avenir Next Cyr W04 Regular"),
    legend.title = element_blank(),
    legend.text = element_text(size=18),
    legend.key.height = unit(3, "cm"),
    legend.key.width = unit(3, "cm"),
    legend.position = "right",
    axis.text = element_text(size=20),
    axis.text.x = element_text(angle = 45, size=20, hjust=1),
    axis.title = element_text(size=25),
    plot.caption = element_text(size=20, hjust=0),
    strip.background = element_blank(),
    strip.text = element_text(size=20, colour="black")
  ) 
ggsave(mod_foodsec_plot, file="strat_by_foodsec_gee.png", height=8, width =12)



## >CONTINUOUS INTERACTION EFFECTS #################################################################################


## >>Modifier: Food Poverty ------------------------------------------------------------------------------------

## Running Interaction Tests ##
meas_pov = runGEEIntxCont(outcome="logmev", exposures=exp_list, adjustment=adjset_meas, modifier="poverty_food_i", df=dr_meas)
tet_pov = runGEEIntxCont(outcome="logtet", exposures=exp_list, adjustment=adjset_tet, modifier="poverty_food_i", df=dr_tet)
hib_pov = runGEEIntxCont(outcome="loghib", exposures=exp_list, adjustment=adjset_tet, modifier="poverty_food_i", df=dr_tet)

intxtests_poverty_food_i_cont = rbind(meas_pov, tet_pov, hib_pov)
rm(meas_pov, tet_pov, hib_pov)

intxtests_poverty_food_i_cont$dv = factor(intxtests_poverty_food_i_cont$dv,
                                     levels = c("logmev", "logtet", "loghib"),
                                     labels = c("Measles", "Tetanus", "HIB"))



## Getting Stratified Estimates ##
meas_foodpov = getGEEStratEstsCont(outcome="logmev", exposures=exp_list, adjustment=adjset_meas, 
                               stratby="poverty_food_i", l1="Above poverty level", l2="At or below poverty level", df=dr_meas)

tet_foodpov = getGEEStratEstsCont(outcome="logtet", exposures=exp_list, adjustment=adjset_tet, 
                              stratby="poverty_food_i", l1="Above poverty level", l2="At or below poverty level", df=dr_tet)

hib_foodpov = getGEEStratEstsCont(outcome="loghib", exposures=exp_list, adjustment=adjset_tet, 
                              stratby="poverty_food_i", l1="Above poverty level", l2="At or below poverty level", df=dr_tet)

mod_foodpov_cont = rbind(meas_foodpov, tet_foodpov, hib_foodpov)

mod_foodpov_cont$dv = factor(mod_foodpov_cont$dv,
                        levels = c("logmev", "logtet", "loghib"),
                        labels = c("Measles", "Tetanus", "HIB"))


mod_foodpov_cont$level = factor(mod_foodpov_cont$level, 
                           levels = c("Above poverty level", "At or below poverty level"),
                           labels = c("Above poverty level", "At or below poverty level"))


## Plotting ## 
produceStratPlot(d1=mod_foodpov_cont[mod_foodpov_cont$level=="Above poverty level",], 
                 d2=intxtests_poverty_food_i_cont, d3=mod_foodpov_cont, name="strat_foodpov_cont_gee")







## >>Modifier: Child Sex ------------------------------------------------------------------------------------


## Running Interaction Tests ##
meas_csex = runGEEIntxCont(outcome="logmev", exposures=exp_list, adjustment=adjset_meas[!(adjset_meas=="csex")], modifier="csex", df=dr_meas)
tet_csex = runGEEIntxCont(outcome="logtet", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="csex")], modifier="csex", df=dr_tet)
hib_csex = runGEEIntxCont(outcome="loghib", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="csex")], modifier="csex", df=dr_tet)

intxtests_csex_cont = rbind(meas_csex, tet_csex, hib_csex)
rm(meas_csex, tet_csex, hib_csex)

intxtests_csex_cont$dv = factor(intxtests_csex_cont$dv,
                           levels = c("logmev", "logtet", "loghib"),
                           labels = c("Measles", "Tetanus", "HIB"))



## Getting Stratified Estimates ##
meas_csex = getGEEStratEstsCont(outcome="logmev", exposures=exp_list, adjustment=adjset_meas[!(adjset_meas=="csex")], 
                            stratby="csex", l1="Boy", l2="Girl", df=dr_meas)

tet_csex = getGEEStratEstsCont(outcome="logtet", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="csex")], 
                           stratby="csex", l1="Boy", l2="Girl",  df=dr_tet)

hib_csex = getGEEStratEstsCont(outcome="loghib", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="csex")], 
                           stratby="csex", l1="Boy", l2="Girl",  df=dr_tet)

mod_csex_cont = rbind(meas_csex, tet_csex, hib_csex)

mod_csex_cont$dv = factor(mod_csex_cont$dv,
                     levels = c("logmev", "logtet", "loghib"),
                     labels = c("Measles", "Tetanus", "HIB"))

mod_csex_cont$level = factor(mod_csex_cont$level, 
                        levels = c("Boy", "Girl"),
                        labels = c("Boy ", "Girl")) #again, stupid way to make sure spacing works out, but do not have time to do it better right now

## Producing strata plots ##
produceStratPlot(d1=mod_csex_cont[mod_csex_cont$level  == "Boy ",], 
                 d2=intxtests_csex_cont, d3=mod_csex_cont, name="strat_csex_cont_gee", lims=c(-30, 35), ast=20, cuts=c(-20,0,20))



## >>Modifier: Maternal HIV ------------------------------------------------------------------------------------


## Running Interaction Tests ##
meas_hiv = runGEEIntxCont(outcome="logmev", exposures=exp_list, adjustment=adjset_meas[!(adjset_meas=="hiv_combined_i")], modifier="hiv_combined_i", df=dr_meas)
tet_hiv = runGEEIntxCont(outcome="logtet", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="hiv_combined_i")], modifier="hiv_combined_i", df=dr_tet)
hib_hiv = runGEEIntxCont(outcome="loghib", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="hiv_combined_i")], modifier="hiv_combined_i", df=dr_tet)

intxtests_hiv_cont = rbind(meas_hiv, tet_hiv, hib_hiv)
rm(meas_hiv, tet_hiv, hib_hiv)

intxtests_hiv_cont$dv = factor(intxtests_hiv_cont$dv,
                               levels = c("logmev", "logtet", "loghib"),
                               labels = c("Measles", "Tetanus", "HIB"))



## Getting Stratified Estimates ##
meas_hiv = getGEEStratEstsCont(outcome="logmev", exposures=exp_list, adjustment=adjset_meas[!(adjset_meas=="hiv_combined_i")], 
                               stratby="hiv_combined_i", l1="HIV Negative", l2="HIV Positive", df=dr_meas)

tet_hiv = getGEEStratEstsCont(outcome="logtet", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="hiv_combined_i")], 
                              stratby="hiv_combined_i", l1="HIV Negative", l2="HIV Positive",  df=dr_tet)

hib_hiv = getGEEStratEstsCont(outcome="loghib", exposures=exp_list, adjustment=adjset_tet[!(adjset_tet=="hiv_combined_i")], 
                              stratby="hiv_combined_i", l1="HIV Negative", l2="HIV Positive", df=dr_tet)

mod_hiv_cont = rbind(meas_hiv, tet_hiv, hib_hiv)

mod_hiv_cont$dv = factor(mod_hiv_cont$dv,
                         levels = c("logmev", "logtet", "loghib"),
                         labels = c("Measles", "Tetanus", "HIB"))


## Producing strata plots ##
produceStratPlot(d1=mod_hiv_cont[mod_hiv_cont$level  == "HIV Negative",], 
                 d2=intxtests_hiv_cont, d3=mod_hiv_cont, name="strat_hiv_cont_gee", lims=c(-0.4,0.5), cuts=c(-0.4,-0.2, 0, 0.2,0.4))





