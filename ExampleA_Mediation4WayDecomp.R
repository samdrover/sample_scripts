## >FILE INFO ############################################################################

# AUTHOR:               Samantha Drover
# PURPOSE:              Models for neurodevelopment paper, including mediation  
# LAST MODIFIED:        2019-10-04


# DATA DICTIONARY:      
    # [Redacted path to data dictionary]


# MODIFICATIONS:
    # 2019-08-01: First created 
    # 2019-08-26: Adding in plots of the mediation effect. Also restricting to mediation effects where temporality is clearly established
    # 2019-10-04: Focusing mediation on specific timepoints. Also integrating the 3.5 year CBCL data and visualizing that. 

# HELPFUL LINKS:
    # https://cran.r-project.org/web/packages/ggalluvial/vignettes/ggalluvial.html
    # https://www.rdocumentation.org/packages/mediation/versions/4.4.7/topics/mediate
    # https://stats.stackexchange.com/questions/104692/comprehending-output-from-mediation-analysis-in-r
    # https://github.com/Unalmut/4way-decomposition 
          # Above link includes all the functions I am using for the mediation


# PAPER/ONLINE RESOURCE USED (TO CITE IN PAPER): 

    # https://github.com/Unalmut/4way-decomposition 
    
    #Citation
    #If you use this R code, please refer to the paper (including supplement) of VanderWeele (2014): "A unification of mediation and interaction: a four-way decomposition."
    
    #License
    #This project is licensed under GNU GPL v3.
    
    #Authors
    #Unal Mutlu (Departments of Epidemiology and Ophthalmology, Erasmus MC, Rotterdam, Netherlands)
    #Gennady V. Roshchupkin (Department of Epidemiology, Radiology and Medical Informatics, Erasmus MC, Rotterdam, Netherlands)
    #Gautam Sajeev (Department of Epidemiology, Harvard T.H. Chan School of Public Health, Boston, Massachusetts)
    #M. Kamran Ikram (Departments of Epidemiology and Neurology, Erasmus MC, Rotterdam, Netherlands)
    #M. Arfan Ikram (Departments of Epidemiology, Neurology and Radiology, Erasmus MC, Rotterdam, Netherlands
     #               Contacts
     #                If you have any questions/suggestions/comments or problems do not hesitate to contact us: u.mutlu@erasmusmc.nl, g.roshchupkin@erasmusmc.nl


## >PATHS & LIBRARIES #############################################################################

rm(list=ls()) #clears workspace

library(cowplot)
library(devtools)
library(dplyr)
library(extrafont)
#library(flipPlots) #not available for newest version of R
library(ggalluvial)
library(ggplot2)
library(ghibli)
library(gridExtra)
library(psych)
library(mediation)
library(rms)
library(dplyr)
library(boot)

#extrafont::font_import() #for importing fonts - takes a while
extrafont::loadfonts(device="win")

neuro = read.csv("[REDACTED PATH TO NEURODEVELOPMENT DATA - e.g., behavioral assessments")
neuro_labels = read.csv("[REDACTED PATH TO NEURODEVELOPMENT DATA LABELS]")
nih = read.csv("[REDACTED PATH TO NIH ASSESSMENT DATA]")
cbcl35 = read.csv("[REDACTED PATH TO CBCL AT 3.5 YEARS DATA]")
cbcl35 = cbcl35[, colnames(cbcl35) != "X"] #dropping meaningless variable/irrelevant variable

setwd("[REDACTED PATH TO OUTPUT DIRECTORY")


## >MERGING DATA #####################################################################################


colnames(nih)[colnames(nih)=="pin"] = "hsn" #first renaming  pin to be hsn

neuro = neuro[order(neuro$hsn),]
nih = nih[order(nih$hsn),]
df = merge(neuro, nih, by="hsn", all=T)

nih$hsn[!(nih$hsn %in% neuro$hsn)] #note HSN that are in the NIH toolbox data but are not in the main dataset: 1064  3132 11765 11953

length(df$hsn[!is.na(df$attn_age.corrected.standard.score)])
length(df$hsn[!is.na(df$cbcl_adhd_2yr)])


cbcl35 = cbcl35[order(cbcl35$hsn),]
df = df[order(df$hsn), ]
d = merge(df, cbcl35, all=T)

df$hsn[!(df$hsn %in% cbcl35$hsn)] #list of HSNs in the main dataset that do not have the 3.5 year CBCL data (quite a lot)


## >RECODING/CREATING VARIABLES ############################################################################

## >>Food security --------------------------------------------------------------------------

## Coding 3-level as factors ##
d$food_cat_dqx =   factor(d$food_cat_dqx, levels=c(1,2,3), labels=c("Secure", "Low", "Very Low"))
d$food_cat_i_1yr = factor(d$food_cat_i_1yr, levels=c(1,2,3), labels=c("Secure", "Low", "Very Low"))
d$food_cat_i_2yr = factor(d$food_cat_i_2yr, levels=c(1,2,3), labels=c("Secure", "Low", "Very Low"))


## >>Poverty Recoding -------------------------------------------------------------------

## Delivery timepoint ##
table(d$poverty_food_i_dqx)

d$poverty_food_i_dqx = factor(d$poverty_food_i_dqx,
                              levels = c(0,1),
                              labels = c("Above", "At or Below")) 

## Year 1 timepoint ##
table(d$poverty_food_1yr)

d$poverty_food_1yr = factor(d$poverty_food_1yr,
                              levels = c(0,1),
                              labels = c("Above", "At or Below"))


## Year 2 timepoint ##
table(d$poverty_food_2yr)

d$poverty_food_2yr = factor(d$poverty_food_2yr,
                            levels = c(0,1),
                            labels = c("Above", "At or Below"))



## >>Wealth Index Recoding ---------------------------------------------------------------------

table(d$wealth_low) #categorical version
summary(d$wealth_i) #continuous 

psych::describeBy(d$wealth_i, d$wealth_low) #lower wealth scores are in wealth_low=1

table(d$poverty_food_i_dqx, d$wealth_low) #looks like 1 is probably low wealth 

d$wealth_low = factor(d$wealth_low,
                      levels = c(0,1),
                      labels = c("High Wealth", "Low Wealth"))



## >>Child Diet Diversity -------------------------------------------------------------------

## Function to Relevel & Reformat Diet Values to Sum More easily ##
relevelFoodItem = function(input) {
  
  input[input=="Refused"] = NA #setting refused to missing 
  input = droplevels(input) #dropping the refused level
  
  output = factor(input,
                  levels = c("No", "Yes"),
                  labels = c(0,1))
  
  output = as.numeric(as.character(output)) #changing to numeric so can sum the 1s 
  
  return(output)
}


## Year 1 Diversity ##
dietdiversity_yr1 = rowSums(cbind(relevelFoodItem(d$bqcmz_1yr),
                                  relevelFoodItem(d$bqcfrt_1yr),
                                  relevelFoodItem(d$bqcvg_1yr),
                                  relevelFoodItem(d$bqcbn_1yr),
                                  relevelFoodItem(d$bqcptsw_1yr),
                                  relevelFoodItem(d$bqcsp_1yr),
                                  relevelFoodItem(d$bqcpnut_1yr),
                                  relevelFoodItem(d$bqxmash_1yr),
                                  relevelFoodItem(d$bqcpmpkn_1yr),
                                  relevelFoodItem(d$bqcbhl_1yr),
                                  relevelFoodItem(d$bqcmkhh_1yr),
                                  relevelFoodItem(d$bqcredmt_1yr),
                                  relevelFoodItem(d$bqcchklv_1yr),
                                  relevelFoodItem(d$bqceggs_1yr),
                                  relevelFoodItem(d$bqcfish_1yr),
                                  relevelFoodItem(d$bqcmlktea_1yr)
                                  ))
summary(dietdiversity_yr1)


## Year 2 Diversity ##
dietdiversity_yr2 = rowSums(cbind(relevelFoodItem(d$bqcmz_2yr),
                                  relevelFoodItem(d$bqcfrt_2yr),
                                  relevelFoodItem(d$bqcvg_2yr),
                                  relevelFoodItem(d$bqcbn_2yr),
                                  relevelFoodItem(d$bqcptsw_2yr),
                                  relevelFoodItem(d$bqcsp_2yr),
                                  relevelFoodItem(d$bqcpnut_2yr),
                                  relevelFoodItem(d$bqxmash_2yr),
                                  relevelFoodItem(d$bqcpmpkn_2yr),
                                  relevelFoodItem(d$bqcbhl_2yr),
                                  relevelFoodItem(d$bqcmkhh_2yr),
                                  relevelFoodItem(d$bqcredmt_2yr),
                                  relevelFoodItem(d$bqcchklv_2yr),
                                  relevelFoodItem(d$bqceggs_2yr),
                                  relevelFoodItem(d$bqcfish_2yr),
                                  relevelFoodItem(d$bqcmlktea_2yr),
                                  relevelFoodItem(d$bqccan_2yr),
                                  relevelFoodItem(d$bqccand_2yr)    
                                  ))
summary(dietdiversity_yr2) #lines up with previous paper [SPECIFICS REDACTED]

d$dietdiversity_yr1 = dietdiversity_yr1
d$dietdiversity_yr2 = dietdiversity_yr2

cbp2 = c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #colour-blind friendly colour palette


## >VARIABLE SETS ##################################################################################

bayleyvars = c(as.character(neuro_labels$V1[grep("bayley|Bayley", neuro_labels$new_labels)]))

bayleyvars_1yr = bayleyvars[c(grep("_1yr", bayleyvars))]
bayleyvars_2yr = bayleyvars[c(grep("_2yr", bayleyvars))]

cbclvars = c(colnames(d)[grep("cbcl", colnames(d))]) #list of all the CBCL variables 

cbclvars_2yr = cbclvars[c(grep("_2yr", cbclvars))]
cbclvars_3yr = cbclvars[c(grep("_3.5yr", cbclvars))]
cbclvars_3yr =c("cbcl_agg_3.5yr", "cbcl_anx_dep_3.5yr", "cbcl_att_3.5yr", "cbcl_with_3.5yr",
                "cbcl_adhd_3.5yr", "cbcl_aff_3.5yr", "cbcl_anx_3.5yr", "cbcl_odd_3.5yr", "cbcl_pdd_3.5yr",
                "cbcl_er_3.5yr","cbcl_som_3.5yr","cbcl_sleep_3.5yr"  ) #ordering so in same order as the 2 year              
  
  
nihvars  = c("peg_fully.corrected.t.scores.dominant", "speed_fully.corrected.t.score", 
             "attn_fully.corrected.t.score", "card_fully.corrected.t.score")

foodsecvars = c("insecure_low_dqx", "insecure_low_i_1yr", "insecure_low_i_2yr")

foodsecvars_alt = c("food_cat_dqx", "food_cat_i_1yr", "food_cat_i_2yr")

foodpovvars = c("poverty_food_i_dqx", "poverty_food_1yr", "poverty_food_2yr")

wealthvars = c("wealth_i", "wealth_low")

matnutrvars = c("femg", "alcoholg", "haemironmg", "nonhaemironmg", "mgmg", "znmg", "semcg",
                "imcg", "vitaminb6mg", "vitaminb12mcg", "folatemcg", "tryptophang")

t1_vars=c("mage_cat", "edu_cat", "poverty_food_i", "insecure_low_dqx", "parity2",
          "hiv_combined_i", "smoke_ever", "alcohol", "mal_iom_high_act",
          "csex", "bwt_cat", "ga_cat", "deliv_method") #add other variables -- for now use these 
            #look at Eskenazi paper for which variables to add 
            #especialyl need to add the breastfeeding variable


adj1 = c("edu_cat", "married", "rms::rcs(mage,3)", "parity2")
adj2 = c("edu_cat", "married", "mage", "parity2")
adj3 = c("mage")


preddf_adj2 = data.frame(
  edu_cat = "<12th Grade",
  married = "No",
  mage = 26.43,
  parity2 = "0"
) #prediction dataframe for the second adjustment set 

## >RESTRICTIN DATA FRAME ######################################################################

# For now , restricting dataframe to any child with at least one neurodevelopment outcome

neurodata = d[,colnames(d) %in% c(bayleyvars, cbclvars, nihvars)] #getting dataset of only neuro data
neurodata$rowsum =rowSums(neurodata, na.rm=T) #creating a sum across rows -- should be NA when there are no available neuro scores

dr = d[neurodata$rowsum>0,] #restricting to subset where have at least one entry for a given neuro variable 



## >FUNCTIONS ###########################################################################


## For Getting Descriptives of Continuous Outcomes ##
simpleDescripT = function(varlist, data = dr) {
  
  iloop = list() 
  for (i in 1:length(varlist)) {
    
    mean_sd = paste(
      round(mean(data[,varlist[i]], na.rm=T),2),
      " (",
      round(sd(data[,varlist[i]], na.rm=T),2),
      ")",
      sep=""
    )
    
    min_max = paste(
      round(min(data[,varlist[i]], na.rm=T),2),
      "-",
      round(max(data[,varlist[i]], na.rm=T),2),
      sep=""
    )
    
    nobs = sum(!is.na(data[,varlist[i]]))
    
    variable = varlist[i]
    
    iloop[i] = list(
      as.data.frame(cbind(variable, nobs, mean_sd, min_max))
    )
    
    
  }
  ibound = do.call("rbind", iloop)
  return(ibound)
  
}


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


## For making table 1 including tests of association ##
catTable1 = function(df, varlist, stratby) { 
  
  t1 = list()
  for (i in seq_along(varlist)) {
    
    strat = levels(stratby)
    iv = df[,colnames(df)==varlist[i]]
    
    ## for stratum 1 ##
    levels = as.data.frame(table(df[stratby==strat[1],colnames(df)==varlist[i]])) #Getting frequencies
    sum=sum(levels$Freq)
    levels$Pct = round(levels$Freq/sum*100,1) #Getting percentages
    
    levels$Var1 = as.character(levels$Var1) #Converting all to characters to make merge easier
    levels$Freq = as.character(levels$Freq)
    levels$Pct = as.character(levels$Pct)
    
    name = as.data.frame(cbind("Var1" = varlist[i], "Freq" = NA, "Pct" = NA)) #getting name of the variable to put on top
    name$Var1 = as.character(name$Var1)
    name$Freq = as.character(name$Freq)
    name$Pct = as.character(name$Pct)
    
    ## for stratum 2 ##
    levels2 = as.data.frame(table(df[stratby==strat[2],colnames(df)==varlist[i]])) #Getting frequencies
    sum=sum(levels2$Freq)
    levels2$Pct = round(levels2$Freq/sum*100,1) #Getting percentages
    
    levels2$Var1 = as.character(levels2$Var1) #Converting all to characters to make merge easier
    levels2$Freq = as.character(levels2$Freq)
    levels2$Pct = as.character(levels2$Pct)
    
    ## chi square test ##
    c = chisq.test(iv,stratby)
    chi_p = paste(round(c$statistic,1), " (", round(c$p.value,2), ")", sep="")
  
    ## binding all together ##
    t1[i] = list(
      cbind(rbind(name,levels), 
            rbind(name,levels2),
            chi_p)
    )
  }
  
  out = as.data.frame(do.call("rbind", t1))
  
  return(out)
  
}


## For producing main effects (categorical exposures) ##
runCatExpReg = function(exposure, outcome, adjset=adj2, preddf=preddf_adj2, df=dr){ 
  
  iloop = list()
  for(i in 1:length(outcome)) {
    
    f = as.formula(paste(outcome[i], "~", exposure, "+", paste(adjset, collapse="+") )) #regression formula
    
    m = lm(data=df, formula=f) #running regression
    
    av = as.data.frame(anova(m))
    p = av$`Pr(>F)`[rownames(av)==exposure] #getting p value for the overall effect of the exposure
    siglabel = ifelse(p<.05, "*", "")
    
    exp_levels = levels(df[,exposure]) #getting categorical values of factor
    
    forpred = as.data.frame(cbind(exp_levels, preddf))
    forpred$exp_levels = factor(forpred$exp_levels, levels=exp_levels, labels=exp_levels)
    colnames(forpred)[1] = exposure #renaming the exposure levels as the exposure name 
    
    outpred = as.data.frame(predict.lm(m, forpred, interval="confidence")) #getting predicted values at each level
    outpred$outcome = outcome[i]
    outpred$exposure = exposure
    outpred$siglabel = siglabel
    outpred$exp_level = exp_levels
    outpred$exp_level = factor(outpred$exp_level, levels=exp_levels, labels=exp_levels)
    
    outpred = outpred[, c("exposure", "outcome", "exp_level", "fit", "lwr", "upr", "siglabel")] #ordering output 
    
    iloop[i] = list(outpred)
    
  }
  ibound = do.call("rbind", iloop)
  return(ibound)
  
}


## For producing main effects (continuous exposures) and outputting results neatly ##
runContExpReg = function(exposure, outcome, adjset=adj2, df=dr){ 
  
  iloop = list()
  for(i in 1:length(outcome)) {
    
    f = as.formula(paste(outcome[i], "~", exposure, "+", paste(adjset, collapse="+") )) #regression formula
    
    m = lm(data=df, formula=f) #running regression
    
    s = as.data.frame(summary(m)$coefficients) #getting coefficients
    s2 = s[rownames(s)==exposure,]
    s2$est_95 = paste(round(s2$Estimate,2), " [",
                      round((s2$Estimate - (1.96*s2$`Std. Error`) ),2 ),
                      ", ",
                      round((s2$Estimate + (1.96*s2$`Std. Error`) ),2 ),
                      "]",sep="") #getting the estimate and 95% CI
      
    s2$outcome = outcome[i]
    s2$exposure = exposure
    
    s3 = s2[, c("exposure", "outcome", "est_95")] 
    
    iloop[i] = list(s3)
    
  }
  ibound = do.call("rbind", iloop)
  return(ibound)
  
}


## Function to Save MEs after Running ME ##
saveMEPlot = function(gd, plotname, ylabel, height=10){
  p = 
    ggplot() + 
    geom_pointrange(data=gd, aes(x=exp_level, y=fit, ymin=lwr, ymax=upr), 
                    position=position_dodge(0.6)) +
    geom_text(data=gd, aes(x=Inf, label = siglabel, y = Inf), size=8, vjust=1, hjust=1, colour="red") +
    facet_grid(outcome~exposure, scales="free") +
    labs(
      x="",
      y=ylabel
    ) +
    theme(
      strip.text.y = element_text(angle=0)
    )
  
  ggsave(p, file=paste(plotname, ".png", sep=""), height=height, width=6)
}


## >>Four-Way Decompisition Mediation Functions --------------------------------------------------

#The below functions call each other, so need all three. Call the third when using (the third calls the previous two)

## Restrict the Dataframe to be Complete Case for Exposure, Mediator, and Outcome ##
restrictData = function(df, exp, med, out) {
  
  expvect = df[, c(exp)] #creating vectors of exposure, mediator, and outcome
  medvect = df[, c(med)] 
  outvect = df[, c(out)]
  
  newdf = df[!is.na(medvect) & !is.na(expvect) & !is.na(outvect),] #using dataframe where not missing mediator or exposure or outcome
  
  return(newdf)
  
}


## Mediation 4-Way Decomposition (Based on VanderWeele 2014 supplement formulae) ##
mediation4WayDecomp = function(data, indices, exp, med, out, adjset, index_a, mstar) {
  
  newdf = data[indices,]
  
  ## Formulae for Regression Models ##
  f_mediator = as.formula(paste(med, "~", paste(c(exp,adjset), collapse="+")))
  f_outcome = as.formula(paste(out, "~", exp, "*", med, "+ ",paste(c(med, exp, adjset), collapse="+")))
  
  
  ## Regression Models for Mediator and Outcome ##
  m_mediator = glm(data = newdf, formula = f_mediator, family=binomial(link="logit"))
  m_outcome = lm(data = newdf, formula = f_outcome)
  
  
  ## Getting Coefficients from Models ##
  med_coefs = as.data.frame(coef(m_mediator)) #mediator model coefficients
  med_coefs$param = rownames(med_coefs)
  
  out_coefs = as.data.frame(coef(m_outcome))
  out_coefs$param = rownames(out_coefs)
  # Note: Coefficients from mediator model will be refered to as b (for beta)
  #Coefficients from outcome model will be refered t0 as t (for theta)
  
  
  ## Cycling Through Covariate Set to Get Mode/Mean Values ##
  
  b2ck = list() #list to add the products of the betas and covariate values for each covariate 
  for (k in seq_along(adjset)) {
    
    
    if (is.numeric(newdf[,c(adjset[k])])) { 
      
      b2k = med_coefs$`coef(m_mediator)`[grep(adjset[k], med_coefs$param)] #pulling out coefficient for continuous covariate
      
      b2ck[k] = list(b2k * mean(newdf[,c(adjset[k])], na.rm=T)) #beta * covariate mean 
      
    }
    
    else if (is.factor(newdf[,c(adjset[k])])) {
      
      coefs = med_coefs[grep(adjset[k], med_coefs$param),] #get only the coefficients pertaining to the factor covariate
      
      mode = Mode(newdf[,c(adjset[k])]) #getting mode of the categorical variable
      
      mode_match_param = grep(mode, coefs$param) #getting the index where the mode matches one of the dummy variables
      
      if (length(mode_match_param)==0) {
        values = rep(0,nrow(coefs))
        b2ck[k] = sum(values*coefs$`coef(m_mediator)`)
      }
      
      else {
        values = rep(0,nrow(coefs))
        values[mode_match_param] = 1
        b2ck[k] = sum(values*coefs$`coef(m_mediator)`)
      }
      
      
    }
    
  }
  b2cksum = sum(do.call("cbind", b2ck)) #summing the product of the covariate coefficients times the covariate mean or mode values
  
  
  ## Pulling Coefficients for Exposure and Mediator from Models ##
  t1 = out_coefs$`coef(m_outcome)`[out_coefs$param %in% out_coefs$param[grep(exp, out_coefs$param)] & 
                                     !(out_coefs$param %in% out_coefs$param[grep(med,out_coefs$param)])] #coefficient for exposure in exposure model
  
  t2 = out_coefs$`coef(m_outcome)`[out_coefs$param %in% out_coefs$param[grep(med, out_coefs$param)] & 
                                     !(out_coefs$param %in% out_coefs$param[grep(exp,out_coefs$param)])] #coefficient for mediator in exposure model
  
  t3 = out_coefs$`coef(m_outcome)`[out_coefs$param %in% out_coefs$param[grep(exp, out_coefs$param)] & 
                                     (out_coefs$param %in% out_coefs$param[grep(med,out_coefs$param)])] #coefficient for exposure X mediatior interaction in exposure model
  
  b0 = med_coefs$`coef(m_mediator)`[grep("Intercept", med_coefs$param)] #intercept coefficient for mediator model
  
  b1 = med_coefs$`coef(m_mediator)`[grep(exp, med_coefs$param)]
  
  
  ## Getting a and a* values ##
  if (length(grep(index_a, out_coefs$param))>0) {
    a = 1
    astar = 0 #basically, if the parameter is coded so that the dummy variable is giving the coefficient is giving the comparison of the specified index to the referent, this is what we want
  } else {
    a = 0
    astar = 1
  } 
  
  
  ## Computing 4-Way Decomposition Effects ##
  
  cde = (t1 + (t3*mstar))*(a - astar) #Controlled Direct Effect at mediator = mstar
  
  pie1 = t2 + (t3*astar) #getting each part of the PIE equation separately 
  pie2 = exp(b0 + (b1*a) + b2cksum)
  pie3 = exp(b0 + (b1*astar) + b2cksum)
  pie = pie1*((pie2/(1+pie2)) - (pie3/(1+pie3))) #Pure Indirect Effect - roughly lines up with ACMEcontrol from mediate()
  
  intref1 = t3*(a - astar)
  intref2 = exp(b0 + (b1*astar) + b2cksum)
  intref = intref1*((intref2/(1+intref2)) - mstar) #eference Interaction
  
  intmed1 = t3*(a - astar)
  intmed2 = exp(b0 + (b1*a) + b2cksum)
  intmed3 = exp(b0 + (b1*astar) + b2cksum) #Mediated Interaction
  intmed = intmed1*((intmed2/(1+intmed2))-(intmed3/(1+intmed3))) #lines up with NIE-PIE (ACMEtx - ACMEcontrol from mediate() ) 
  
  te = cde + pie + intref + intmed #total effect slightly different than what got from mediate()
  
  ## Computing the Proportion Mediated/Eliminated ##
  prop_cde = cde/te
  prop_intref = intref/te
  prop_intmed = intmed/te
  prop_pie = pie/te
  prop_med = (pie + intmed)/te
  prop_int = (intmed + intref)/te
  prop_elm = (pie + intmed + intref)/te
  
  ## Putting All Results Together ##
  output = c(te, cde, intref, intmed, pie, prop_cde, prop_intref, prop_intmed, prop_pie, prop_med, prop_int, prop_elm)
  
  return(output)
}


## Mediation 4-Way Decomposition (Based on VanderWeele 2014 supplement formulae) ##
mediation4WayDecompCrude = function(data, indices, exp, med, out, index_a, mstar) {
  
  newdf = data[indices,]
  
  ## Formulae for Regression Models ##
  f_mediator = as.formula(paste(med, "~", exp))
  f_outcome = as.formula(paste(out, "~", exp, "*", med, "+ ",paste(c(med, exp), collapse="+")))
  
  
  ## Regression Models for Mediator and Outcome ##
  m_mediator = glm(data = newdf, formula = f_mediator, family=binomial(link="logit"))
  m_outcome = lm(data = newdf, formula = f_outcome)
  
  
  ## Getting Coefficients from Models ##
  med_coefs = as.data.frame(coef(m_mediator)) #mediator model coefficients
  med_coefs$param = rownames(med_coefs)
  
  out_coefs = as.data.frame(coef(m_outcome))
  out_coefs$param = rownames(out_coefs)
  # Note: Coefficients from mediator model will be refered to as b (for beta)
  #Coefficients from outcome model will be refered t0 as t (for theta)
  
  
  ## Cycling Through Covariate Set to Get Mode/Mean Values -- no covariates, therefore sum is 0? ##
  # 
  # b2ck = list() #list to add the products of the betas and covariate values for each covariate 
  # for (k in seq_along(adjset)) {
  #   
  #   
  #   if (is.numeric(newdf[,c(adjset[k])])) { 
  #     
  #     b2k = med_coefs$`coef(m_mediator)`[grep(adjset[k], med_coefs$param)] #pulling out coefficient for continuous covariate
  #     
  #     b2ck[k] = list(b2k * mean(newdf[,c(adjset[k])], na.rm=T)) #beta * covariate mean 
  #     
  #   }
  #   
  #   else if (is.factor(newdf[,c(adjset[k])])) {
  #     
  #     coefs = med_coefs[grep(adjset[k], med_coefs$param),] #get only the coefficients pertaining to the factor covariate
  #     
  #     mode = Mode(newdf[,c(adjset[k])]) #getting mode of the categorical variable
  #     
  #     mode_match_param = grep(mode, coefs$param) #getting the index where the mode matches one of the dummy variables
  #     
  #     if (length(mode_match_param)==0) {
  #       values = rep(0,nrow(coefs))
  #       b2ck[k] = sum(values*coefs$`coef(m_mediator)`)
  #     }
  #     
  #     else {
  #       values = rep(0,nrow(coefs))
  #       values[mode_match_param] = 1
  #       b2ck[k] = sum(values*coefs$`coef(m_mediator)`)
  #     }
  #     
  #     
  #   }
  #   
  # }
  b2cksum = 0 #summing the product of the covariate coefficients times the covariate mean or mode values
  
  
  ## Pulling Coefficients for Exposure and Mediator from Models ##
  t1 = out_coefs$`coef(m_outcome)`[out_coefs$param %in% out_coefs$param[grep(exp, out_coefs$param)] & 
                                     !(out_coefs$param %in% out_coefs$param[grep(med,out_coefs$param)])] #coefficient for exposure in exposure model
  
  t2 = out_coefs$`coef(m_outcome)`[out_coefs$param %in% out_coefs$param[grep(med, out_coefs$param)] & 
                                     !(out_coefs$param %in% out_coefs$param[grep(exp,out_coefs$param)])] #coefficient for mediator in exposure model
  
  t3 = out_coefs$`coef(m_outcome)`[out_coefs$param %in% out_coefs$param[grep(exp, out_coefs$param)] & 
                                     (out_coefs$param %in% out_coefs$param[grep(med,out_coefs$param)])] #coefficient for exposure X mediatior interaction in exposure model
  
  b0 = med_coefs$`coef(m_mediator)`[grep("Intercept", med_coefs$param)] #intercept coefficient for mediator model
  
  b1 = med_coefs$`coef(m_mediator)`[grep(exp, med_coefs$param)]
  
  
  ## Getting a and a* values ##
  if (length(grep(index_a, out_coefs$param))>0) {
    a = 1
    astar = 0 #basically, if the parameter is coded so that the dummy variable is giving the coefficient is giving the comparison of the specified index to the referent, this is what we want
  } else {
    a = 0
    astar = 1
  } 
  
  
  ## Computing 4-Way Decomposition Effects ##
  
  cde = (t1 + (t3*mstar))*(a - astar) #Controlled Direct Effect at mediator = mstar
  
  pie1 = t2 + (t3*astar) #getting each part of the PIE equation separately 
  pie2 = exp(b0 + (b1*a) + b2cksum)
  pie3 = exp(b0 + (b1*astar) + b2cksum)
  pie = pie1*((pie2/(1+pie2)) - (pie3/(1+pie3))) #Pure Indirect Effect - roughly lines up with ACMEcontrol from mediate()
  
  intref1 = t3*(a - astar)
  intref2 = exp(b0 + (b1*astar) + b2cksum)
  intref = intref1*((intref2/(1+intref2)) - mstar) #eference Interaction
  
  intmed1 = t3*(a - astar)
  intmed2 = exp(b0 + (b1*a) + b2cksum)
  intmed3 = exp(b0 + (b1*astar) + b2cksum) #Mediated Interaction
  intmed = intmed1*((intmed2/(1+intmed2))-(intmed3/(1+intmed3))) #lines up with NIE-PIE (ACMEtx - ACMEcontrol from mediate() ) 
  
  te = cde + pie + intref + intmed #total effect slightly different than what got from mediate()
  
  ## Computing the Proportion Mediated/Eliminated ##
  prop_cde = cde/te
  prop_intref = intref/te
  prop_intmed = intmed/te
  prop_pie = pie/te
  prop_med = (pie + intmed)/te
  prop_int = (intmed + intref)/te
  prop_elm = (pie + intmed + intref)/te
  
  ## Putting All Results Together ##
  output = c(te, cde, intref, intmed, pie, prop_cde, prop_intref, prop_intmed, prop_pie, prop_med, prop_int, prop_elm)
  
  return(output)
}


## Mediation 4-Way Decomposition for Continuous Mediator (Based on VanderWeele 2014 supplement formulae) ##
mediation4WayDecompContM = function(data, indices, exp, med, out, adjset, index_a, mstar) {
  
  newdf = data[indices,]
  
  ## Formulae for Regression Models ##
  f_mediator = as.formula(paste(med, "~", paste(c(exp,adjset), collapse="+")))
  f_outcome = as.formula(paste(out, "~", exp, "*", med, "+ ",paste(c(med, exp, adjset), collapse="+")))
  
  
  ## Regression Models for Mediator and Outcome ##
  m_mediator = lm(data = newdf, formula = f_mediator) #linear model for mediator bc it is continuous
  m_outcome = lm(data = newdf, formula = f_outcome)
  
  
  ## Getting Coefficients from Models ##
  med_coefs = as.data.frame(coef(m_mediator)) #mediator model coefficients
  med_coefs$param = rownames(med_coefs)
  
  out_coefs = as.data.frame(coef(m_outcome))
  out_coefs$param = rownames(out_coefs)
  # Note: Coefficients from mediator model will be refered to as b (for beta)
  #Coefficients from outcome model will be refered t0 as t (for theta)
  
  
  ## Cycling Through Covariate Set to Get Mode/Mean Values ##
  
  b2ck = list() #list to add the products of the betas and covariate values for each covariate 
  for (k in seq_along(adjset)) {
    
    
    if (is.numeric(newdf[,c(adjset[k])])) { 
      
      b2k = med_coefs$`coef(m_mediator)`[grep(adjset[k], med_coefs$param)] #pulling out coefficient for continuous covariate
      
      b2ck[k] = list(b2k * mean(newdf[,c(adjset[k])], na.rm=T)) #beta * covariate mean 
      
    }
    
    else if (is.factor(newdf[,c(adjset[k])])) {
      
      coefs = med_coefs[grep(adjset[k], med_coefs$param),] #get only the coefficients pertaining to the factor covariate
      
      mode = Mode(newdf[,c(adjset[k])]) #getting mode of the categorical variable
      
      mode_match_param = grep(mode, coefs$param) #getting the index where the mode matches one of the dummy variables
      
      if (length(mode_match_param)==0) {
        values = rep(0,nrow(coefs))
        b2ck[k] = sum(values*coefs$`coef(m_mediator)`)
      }
      
      else {
        values = rep(0,nrow(coefs))
        values[mode_match_param] = 1
        b2ck[k] = sum(values*coefs$`coef(m_mediator)`)
      }
      
      
    }
    
  }
  b2cksum = sum(do.call("cbind", b2ck)) #summing the product of the covariate coefficients times the covariate mean or mode values
  
  
  ## Pulling Coefficients for Exposure and Mediator from Models ##
  t1 = out_coefs$`coef(m_outcome)`[out_coefs$param %in% out_coefs$param[grep(exp, out_coefs$param)] & 
                                     !(out_coefs$param %in% out_coefs$param[grep(med,out_coefs$param)])] #coefficient for exposure in exposure model
  
  t2 = out_coefs$`coef(m_outcome)`[out_coefs$param %in% out_coefs$param[grep(med, out_coefs$param)] & 
                                     !(out_coefs$param %in% out_coefs$param[grep(exp,out_coefs$param)])] #coefficient for mediator in exposure model
  
  t3 = out_coefs$`coef(m_outcome)`[out_coefs$param %in% out_coefs$param[grep(exp, out_coefs$param)] & 
                                     (out_coefs$param %in% out_coefs$param[grep(med,out_coefs$param)])] #coefficient for exposure X mediatior interaction in exposure model
  
  b0 = med_coefs$`coef(m_mediator)`[grep("Intercept", med_coefs$param)] #intercept coefficient for mediator model
  
  b1 = med_coefs$`coef(m_mediator)`[grep(exp, med_coefs$param)]
  
  
  ## Getting a and a* values ##
  if (length(grep(index_a, out_coefs$param))>0) {
    a = 1
    astar = 0 #basically, if the parameter is coded so that the dummy variable is giving the coefficient is giving the comparison of the specified index to the referent, this is what we want
  } else {
    a = 0
    astar = 1
  } 
  
  
  ## Computing 4-Way Decomposition Effects ##
  
  cde = (t1 + (t3*mstar))*(a - astar) #Controlled Direct Effect at mediator = mstar
  
  pie1 = (t2*b1) + (t3*b1*astar) #getting each part of the PIE equation separately
  pie2 = (a - astar)
  pie = pie1*pie2 #Pure Indirect Effect 
  
  intref1 = t3*b1
  intref2 = a - astar
  intref = intref1*intref2*intref2 #reference Interaction
  
  intmed1 = t3
  intmed2 = b0 + (b1*astar) + b2cksum - mstar
  intmed3 = a - astar
  intmed = intmed1*intmed2*intmed3 #lines up with NIE-PIE (ACMEtx - ACMEcontrol from mediate() ) 
  
  te = cde + pie + intref + intmed #total effect slightly different than what got from mediate()
  
  ## Computing the Proportion Mediated/Eliminated ##
  prop_cde = cde/te
  prop_intref = intref/te
  prop_intmed = intmed/te
  prop_pie = pie/te
  prop_med = (pie + intmed)/te
  prop_int = (intmed + intref)/te
  prop_elm = (pie + intmed + intref)/te
  
  ## Putting All Results Together ##
  output = c(te, cde, intref, intmed, pie, prop_cde, prop_intref, prop_intmed, prop_pie, prop_med, prop_int, prop_elm)
  
  return(output)
}


## Function to Bootstrap 4-Way Decomposition Statistics (get CI)##
runMediationBoot = function(inputdata, exposure, mediator, outcome, adjustment, index_exp, m_forcde, nsim=5000) {
  
  ## Restrict Data Frame to Complete Case for Mediator, Outcome, Exposure ##
  restricted = restrictData(df=inputdata, exp=exposure, med=mediator, out=outcome) #Covariates are same in both models, so this will make sure models use same data (mediator and outcome models)
  
  ## Boot the Mediation 4-Way Decomposition Statistics ##
  results = boot::boot(data=restricted, statistic= mediation4WayDecomp,
                       R=nsim, ##Increase this to ~10 000 once done testing and ready to run properly
                       exp=exposure, med=mediator, out=outcome,
                       adjset=adjustment, index_a=index_exp, mstar=m_forcde)
  
  ## Preparing to Bootstrap CIs ##
  df_list<-list() #create list to contain the output
  
  estimand_names<-c(
    "Total Effect", "CDE", "INTref", "INTmed",
    "PIE","Proportion CDE","Proportion INTref",
    "Proportion INTmed", "Proportion PIE", "Overall Proportion Mediated",
    "Overall Proportion Attributable to Interaction","Overall Proportion Eliminated"
  ) #list of the statistics to be bootstrapped
  index = length(estimand_names) #index to cycle through the estimand
  
  
  ## Loop to Boostrap CIs ##
  for (i  in 1:index){
    r = boot.ci(results, type = "basic", index=i) 
    estq<-r[[2]]
    cis<- r[[4]]
    colnames(cis)<-NULL
    lcl<-cis[1,4]
    ucl<-cis[1,5]
    estimand<-estimand_names[i]
    rdf<-data.frame(estimand, estq, lcl, ucl)
    names(rdf)<-c("Estimand", "Estimate", "LCL", "UCL")
    df_list[[i]]<-rdf
    #rm(r, estq, cis, colnames, lcl, ucl, estimand, rdf)
  }
  
  ## Combining all Bootstrapped Effects ##
  allresults = do.call("rbind", df_list)
  allresults$exposure = exposure
  allresults$mediator = mediator
  allresults$outcome = outcome
  
  return(allresults) #retyrn these results -- tidy up in a separate program
  
}


## Function to Bootstrap 4-Way Decomposition Statistics (get CI) FOR CRUDE MODEL##
runMediationBootCrude = function(inputdata, exposure, mediator, outcome, index_exp, m_forcde, nsim=5000) {
  
  ## Restrict Data Frame to Complete Case for Mediator, Outcome, Exposure ##
  restricted = restrictData(df=inputdata, exp=exposure, med=mediator, out=outcome) #Covariates are same in both models, so this will make sure models use same data (mediator and outcome models)
  
  ## Boot the Mediation 4-Way Decomposition Statistics ##
  results = boot::boot(data=restricted, statistic= mediation4WayDecompCrude,
                       R=nsim, ##Increase this to ~10 000 once done testing and ready to run properly
                       exp=exposure, med=mediator, out=outcome,
                       index_a=index_exp, mstar=m_forcde)
  
  ## Preparing to Bootstrap CIs ##
  df_list<-list() #create list to contain the output
  
  estimand_names<-c(
    "Total Effect", "CDE", "INTref", "INTmed",
    "PIE","Proportion CDE","Proportion INTref",
    "Proportion INTmed", "Proportion PIE", "Overall Proportion Mediated",
    "Overall Proportion Attributable to Interaction","Overall Proportion Eliminated"
  ) #list of the statistics to be bootstrapped
  index = length(estimand_names) #index to cycle through the estimand
  
  
  ## Loop to Boostrap CIs ##
  for (i  in 1:index){
    r = boot.ci(results, type = "basic", index=i) 
    estq<-r[[2]]
    cis<- r[[4]]
    colnames(cis)<-NULL
    lcl<-cis[1,4]
    ucl<-cis[1,5]
    estimand<-estimand_names[i]
    rdf<-data.frame(estimand, estq, lcl, ucl)
    names(rdf)<-c("Estimand", "Estimate", "LCL", "UCL")
    df_list[[i]]<-rdf
    #rm(r, estq, cis, colnames, lcl, ucl, estimand, rdf)
  }
  
  ## Combining all Bootstrapped Effects ##
  allresults = do.call("rbind", df_list)
  allresults$exposure = exposure
  allresults$mediator = mediator
  allresults$outcome = outcome
  
  return(allresults) #retyrn these results -- tidy up in a separate program
  
}



## Function to Bootstrap 4-Way Decomposition Statistics (get CI) for Continuous Mediator##
runMediationBootContM = function(inputdata, exposure, mediator, outcome, adjustment, index_exp, m_forcde, nsim=5000) {
  
  ## Restrict Data Frame to Complete Case for Mediator, Outcome, Exposure ##
  restricted = restrictData(df=inputdata, exp=exposure, med=mediator, out=outcome) #Covariates are same in both models, so this will make sure models use same data (mediator and outcome models)
  
  ## Boot the Mediation 4-Way Decomposition Statistics ##
  results = boot::boot(data=restricted, statistic= mediation4WayDecompContM,
                       R=nsim, ##Increase this to ~10 000 once done testing and ready to run properly
                       exp=exposure, med=mediator, out=outcome,
                       adjset=adjustment, index_a=index_exp, mstar=m_forcde)
  
  ## Preparing to Bootstrap CIs ##
  df_list<-list() #create list to contain the output
  
  estimand_names<-c(
    "Total Effect", "CDE", "INTref", "INTmed",
    "PIE","Proportion CDE","Proportion INTref",
    "Proportion INTmed", "Proportion PIE", "Overall Proportion Mediated",
    "Overall Proportion Attributable to Interaction","Overall Proportion Eliminated"
  ) #list of the statistics to be bootstrapped
  index = length(estimand_names) #index to cycle through the estimand
  
  ## Loop to Boostrap CIs ##
  for (i  in 1:index){
    r = boot.ci(results, type = "basic", index=i) 
    estq<-r[[2]]
    cis<- r[[4]]
    colnames(cis)<-NULL
    lcl<-cis[1,4]
    ucl<-cis[1,5]
    estimand<-estimand_names[i]
    rdf<-data.frame(estimand, estq, lcl, ucl)
    names(rdf)<-c("Estimand", "Estimate", "LCL", "UCL")
    df_list[[i]]<-rdf
  }
  
  ## Combining all Bootstrapped Effects ##
  allresults = do.call("rbind", df_list)
  allresults$exposure = exposure
  allresults$mediator = mediator
  allresults$outcome = outcome
  
  return(allresults) #retyrn these results -- tidy up in a separate program
  
}


## >>Function to Tidy Four-Way Decomposition Results--------------------------------------------------


## Cleaning Estimates for Table ##
cleanMediationResults = function(results) {
  
  ## Getting the Estimates ##
  ests = results[results$Estimand %in% c("Total Effect", "CDE", "INTref", "INTmed", "PIE"),]
  ests$est = round(ests$Estimate,3)
  ests$est_ci = paste("[", round(ests$LCL,3), ", ", round(ests$UCL,3), "]", sep="")
  
  ## Getting Proportions of Mediated ##
  props = results[results$Estimand %in% c("Proportion CDE", "Proportion INTref", "Proportion INTmed", "Proportion PIE"),]
  props$pct = round(100*props$Estimate,1)
  pctmed = c(NA, props$pct) #adding NA at top to merge with the estimates (which has total effect at the top)
  
  ## Putting Results Together
  out = as.data.frame(cbind(ests$exposure, ests$mediator, ests$outcome,
                            as.character(ests$Estimand), ests$est, ests$est_ci, pctmed))
  colnames(out) = c("Exposure", "Mediator", "Outcome","Estimand", "Estimate", "[95% CI]", "% of Effect")
  
  return(out)
  
}


## Cleaning Estimates for Plotting ##
cleanMediationResultsPlot = function(results) {
  
  ## Getting the Estimates ##
  ests = results[results$Estimand %in% c("Total Effect", "CDE", "INTref", "INTmed", "PIE"),]
  
  ## Getting Proportions of Mediated ##
  props = results[results$Estimand %in% c("Proportion CDE", "Proportion INTref", "Proportion INTmed", "Proportion PIE"),]
  props$pct = round(100*props$Estimate,1)
  pctmed = c(100, props$pct) #adding 100 at top to merge with the estimates (which has total effect at the top -- want that to be the referent size)
  
  ## Putting Results Together
  out = as.data.frame(cbind(ests$exposure, ests$mediator, ests$outcome,
                            as.character(ests$Estimand), ests$Estimate, ests$LCL, ests$UCL, pctmed ))
  colnames(out) = c("Exposure", "Mediator", "Outcome","Estimand", "Estimate", "LCL", "UCL", "%ofEffect")
  out$Estimate = as.numeric(as.character(out$Estimate)) #formtting numeric for easier plotting
  out$LCL = as.numeric(as.character(out$LCL))
  out$UCL = as.numeric(as.character(out$UCL))
  out$`%ofEffect` = as.numeric(as.character(out$`%ofEffect`))
  out$Estimand = factor(out$Estimand,
                        levels = c("Total Effect", "CDE", "INTref", "INTmed", "PIE"),
                        labels = c("Total Effect", "CDE", "INTref", "INTmed", "PIE"))
  
  return(out)
  
}


## For Plotting Mediation Effects ## 
makeMedPlot = function(df, name) {
  
  df$effect = as.factor(ifelse(df$Estimand=="Total Effect", 1, 0))
  df$siglabel = ifelse(df$UCL >0 & df$LCL <0 , " ", "*")
  
  
  ## Creating a Total Value for where to plot asterisks ##
  totals = df[df$Estimand != "Total Effect" & df$`%ofEffect` >= 0,] %>% #restricting to the positive non-total %effects
    group_by(Outcome, model) %>%
    dplyr::summarize(total = sum(`%ofEffect`))
  
  
  ## Sorting Data & Merging ##
  df=df[order(df$model),]
  df=df[order(df$Outcome),]
  
  totals = totals[order(totals$model),]
  totals = totals[order(totals$Outcome),]
  
  df2 = merge(df, totals, by=c("model", "Outcome"))
  
  
  ## Creating a Label for Total Effect ##
  df3 = df2[df2$Estimand=="Total Effect",]
  df3$te = paste(
    round(df3$Estimate,2),
    "\n[",
    round(df3$LCL,2),
    ", \n",
    round(df3$UCL,2),
    "]",
    df3$siglabel,
    sep=""
  )
  
  df3 = df3[,colnames(df3) %in% c("model", "Outcome", "te")]
  
  
  ## Sorting and Merging Again ## 
  df2 = df2[order(df2$model),]
  df2 = df2[order(df2$Outcome),]
  
  df3 = df3[order(df3$model),]
  df3 = df3[order(df3$Outcome),]
  
  df4 = merge(df2, df3, by=c("model","Outcome"))
  
  
  ## Restricting to non-total effects ##
  df5 = df4[df4$Estimand != "Total Effect",]
  
  
  ## Plotting ##
  p1 =
    ggplot(df5, aes(x=Outcome, y=`%ofEffect`)) +
    geom_col(aes(colour=Estimand, fill=Estimand), position=position_stack()) +
    #geom_text(aes(label=df5$siglabel, y=df5$total + 10, colour=Estimand), 
    # size=8, position=position_dodge(0.6)) + #this makes asterisks that float above the bars - really nice :)
    #geom_text(aes(label=df5$te, y=df5$total+80), size=3) + #added the significant labels above the bars
    facet_grid(~model, scales="free") +
    #coord_flip() +
    labs(y = "% of Total Effect") +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_fill_manual(values = c(cbp2[4], cbp2[6], cbp2[7], cbp2[8])) +
    scale_colour_manual(values = c(cbp2[4], cbp2[6], cbp2[7], cbp2[8])) +
    theme_light() +
    theme(
      strip.background = element_rect(fill=NULL),
      strip.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=12),
      axis.text.x = element_text(size=11)
    )
  
  
  p2 = 
    ggplot() +
    geom_pointrange(data=df4, aes(x=Outcome, y=Estimate, ymin=LCL, ymax=UCL, colour=Estimand), position=position_dodge(0.8), size=0.3) +
    geom_text(data =df4, aes(label=siglabel, x=Outcome, y=UCL+0.1, colour=Estimand), size=5, position=position_dodge(0.8), show.legend = F) +
    facet_grid(~model) +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_colour_manual(values = c(cbp2[1], cbp2[4], cbp2[6], cbp2[7], cbp2[8])) +
    labs(
      y = "Effect Estimate [95% CI]"
    ) +
    theme_light() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      strip.text = element_text(size=12),
      strip.background = element_rect(fill="grey24"),
      axis.title.y = element_text(size=12),
      legend.title = element_blank(),
      legend.text = element_text(size=12)
    )
  
  legend = cowplot::get_legend(p2 + theme(legend.direction = "horizontal"))
  
  
  ## Arranging Plots Vertically ##
  plots = cowplot::align_plots(p2 + theme(legend.position = "none"), 
                               p1 + theme(legend.position = "none"), align="v")
  
  povsec_cowplot = cowplot::plot_grid(legend, plots[[1]], plots[[2]], ncol=1, rel_heights = c(0.5, 3,3) )
  
  cowplot::save_plot(povsec_cowplot, file=paste(name,"_cowplot.png", sep=""), ncol=1, base_height=8, base_width = 12)
  
}


## >T1: DESCRIPTIVES #######################################################################

t1 = catTable1(df=dr, varlist=t1_vars, stratby=as.factor(dr$insecure_low_dqx)) #stratified by baseline food security (maternal) 
write.csv(t1, "t1_2019-10-06.csv") 


## >EXPOSURE DESCRIPTION ###################################################################


## >>Poverty X Food Security ---------------------------------------------------------------

library(gmodels)

gmodels::CrossTable(neuro$poverty_food_i_dqx, neuro$insecure_low_dqx, chisq=T,
                    prop.r = F, prop.c = F, prop.t = F, prop.chisq = F)


gmodels::CrossTable(neuro$poverty_food_1yr, neuro$insecure_low_i_1yr, chisq=T,
                    prop.r = F, prop.c = F, prop.t = F, prop.chisq = F)

gmodels::CrossTable(neuro$poverty_food_2yr, neuro$insecure_low_i_2yr, chisq=T,
                    prop.r = F, prop.c = F, prop.t = F, prop.chisq = F)



## >>Food security tables across time --------------------------------------------------------

table(dr$insecure_low_dqx, useNA="always")
table(dr$insecure_low_i_1yr, useNA="always")
table(dr$insecure_low_i_2yr, useNA="always")

addmargins(table(dr$insecure_low_dqx, dr$insecure_low_i_1yr, useNA="always"))
addmargins(table(dr$insecure_low_i_1yr, dr$insecure_low_i_2yr, useNA="always"))


## >>Sankey Plot of Exposure Trajectory ---------------------------------------------------------


## Preparing data ## 
delivery =  c("Secure", "Insecure", "Insecure", "Secure", "Secure", "Insecure", "Insecure", "Secure")
year1 =     c("Secure", "Insecure", "Secure", "Insecure", "Secure", "Insecure", "Secure", "Insecure")
year2 =     c("Secure", "Insecure", "Secure", "Secure", "Insecure", "Secure", "Insecure", "Insecure")
forsankey = as.data.frame(cbind(delivery,year1,year2)) #basically making a dataframe of all possible trajectories


## getting the frequncies of each pattern 3#
sss = as.numeric(nrow(dr[dr$insecure_low_dqx==0 &
                           dr$insecure_low_i_1yr==0 &
                           dr$insecure_low_i_2yr==0 &
                           !is.na(dr$insecure_low_dqx) & !is.na(dr$insecure_low_i_1yr) & !is.na(dr$insecure_low_i_2yr),]))

iii = as.numeric(nrow(dr[dr$insecure_low_dqx==1 &
                           dr$insecure_low_i_1yr==1 &
                           dr$insecure_low_i_2yr==1 &
                           !is.na(dr$insecure_low_dqx) & !is.na(dr$insecure_low_i_1yr) & !is.na(dr$insecure_low_i_2yr),]))

iss = as.numeric(nrow(dr[dr$insecure_low_dqx==1 &
                           dr$insecure_low_i_1yr==0 &
                           dr$insecure_low_i_2yr==0 &
                           !is.na(dr$insecure_low_dqx) & !is.na(dr$insecure_low_i_1yr) & !is.na(dr$insecure_low_i_2yr),]))

sis = as.numeric(nrow(dr[dr$insecure_low_dqx==0 &
                           dr$insecure_low_i_1yr==1 &
                           dr$insecure_low_i_2yr==0 &
                           !is.na(dr$insecure_low_dqx) & !is.na(dr$insecure_low_i_1yr) & !is.na(dr$insecure_low_i_2yr),]))

ssi = as.numeric(nrow(dr[dr$insecure_low_dqx==0 &
                           dr$insecure_low_i_1yr==0 &
                           dr$insecure_low_i_2yr==1 &
                           !is.na(dr$insecure_low_dqx) & !is.na(dr$insecure_low_i_1yr) & !is.na(dr$insecure_low_i_2yr),]))

iis = as.numeric(nrow(dr[dr$insecure_low_dqx==1 &
                           dr$insecure_low_i_1yr==1 &
                           dr$insecure_low_i_2yr==0 &
                           !is.na(dr$insecure_low_dqx) & !is.na(dr$insecure_low_i_1yr) & !is.na(dr$insecure_low_i_2yr),]))

isi = as.numeric(nrow(dr[dr$insecure_low_dqx==1 &
                           dr$insecure_low_i_1yr==0 &
                           dr$insecure_low_i_2yr==1 &
                           !is.na(dr$insecure_low_dqx) & !is.na(dr$insecure_low_i_1yr) & !is.na(dr$insecure_low_i_2yr),]))

sii = as.numeric(nrow(dr[dr$insecure_low_dqx==0 &
                           dr$insecure_low_i_1yr==1 &
                           dr$insecure_low_i_2yr==1 &
                           !is.na(dr$insecure_low_dqx) & !is.na(dr$insecure_low_i_1yr) & !is.na(dr$insecure_low_i_2yr),]))

forsankey$freq = c(sss, iii, iss, sis, ssi, iis, isi, sii)
forsankey$pattern = as.character(c(1:8))
sum(forsankey$freq)



## Using flipPlot SankeyDiagram method ##
flipPlots::SankeyDiagram(data=forsankey,
                         weights=forsankey$freq) #basically come back to this to refine


## Using ggalluvial ##
ggalluvial::is_alluvia_form(forsankey) #oh good

exp_trajectory_plot = 
  ggplot(forsankey, aes(y=freq, axis1=delivery, axis2=year1, axis3=year2)) +
  geom_alluvium(aes(fill=pattern), width=1/12, show.legend = F) +
  geom_stratum(width=1/10, color=c("black"), size=2) +
  geom_text(stat="stratum", label.strata=T, angle=90, size=9) +
  scale_x_discrete(limits = c("Delivery", "Year 1", "Year 2"), labels=c("\nDelivery", "\nYear 1", "\nYear 2"), 
                   expand=c(0.01, 0.01)) +
  scale_fill_brewer(type="qual", palette="Set1") + 
  scale_y_continuous(expand=c(0,0)) +
  labs(
    y="Count\n",
    x="\nTime Point") +
  theme_light() +
  theme(
    text = element_text(family = "Avenir Next Cyr W04 Regular"),
    legend.title = element_blank(),
    legend.text = element_text(size=20),
    legend.key.height = unit(3, "cm"),
    legend.key.width = unit(3, "cm"),
    legend.position = "top",
    axis.text = element_text(size=20),
    axis.text.x = element_text(size=20),
    axis.title = element_text(size=25)
  ) 
ggsave(exp_trajectory_plot, file="exp_trajectory_plot.png", height=8, width=10)


## >>Associations of Food Security with Other SES Indices -----------------------------------------------------------------


## Delivery timepoint ##
exp_assocs = catTable1(df=dr, varlist = c(foodpovvars,"wealth_low"), stratby=as.factor(dr$insecure_low_dqx))
write.csv(exp_assocs, "pov_vars_by_food_sec_deliv.csv")
rm(exp_assocs)


## Year 1 timepoint ##
exp_assocs = catTable1(df=dr, varlist = c(foodpovvars,"wealth_low"), stratby=as.factor(dr$insecure_low_i_1yr))
write.csv(exp_assocs, "pov_vars_by_food_sec_yr1.csv")
rm(exp_assocs)


## Year 2 timepoint ##
exp_assocs = catTable1(df=dr, varlist = c(foodpovvars,"wealth_low"), stratby=as.factor(dr$insecure_low_i_2yr))
write.csv(exp_assocs, "pov_vars_by_food_sec_yr2.csv")
rm(exp_assocs)



## >OUTCOME DESCRIPTION #######################################################################################


## Mean, SD, Range, N nonmissing ##
write.csv(
  rbind(
    simpleDescripT(varlist=cbclvars_2yr, data=neuro), #restricting to df for 2 year analyses -- somehow get greater N in larger dataframe? 
    simpleDescripT(varlist=cbclvars_3yr, data=cbcl35)),
  "cbcl_descriptives.csv")


## Correlations ##
quickCor = function(sc, vect1, vect2) { #quick function to nicely put correlatiosn together 
  r = cor.test(vect1, vect2)
  cor_ci = paste(round(r$estimate,2),
                 " [", round(r$conf.int[[1]],2),
                 ", ", round(r$conf.int[[2]],2),
                 "]", sep="")
  scale = sc
  out = as.data.frame(cbind(scale, cor_ci))
  
}

cbcl_corrs = rbind(
  quickCor("agg", dr$cbcl_agg_2yr, dr$cbcl_agg_3.5yr),
  quickCor("anx_dep", dr$cbcl_anx_dep_2yr, dr$cbcl_anx_dep_3.5yr),
  quickCor("att", dr$cbcl_att_2yr, dr$cbcl_att_3.5yr),
  quickCor("with", dr$cbcl_with_2yr, dr$cbcl_with_3.5yr),
  quickCor("adhd", dr$cbcl_adhd_2yr, dr$cbcl_adhd_3.5yr),
  quickCor("aff", dr$cbcl_aff_2yr, dr$cbcl_aff_3.5yr),
  quickCor("anx", dr$cbcl_anx_2yr, dr$cbcl_anx_3.5yr),
  quickCor("odd", dr$cbcl_odd_2yr, dr$cbcl_odd_3.5yr),
  quickCor("pdd", dr$cbcl_pdd_2yr, dr$cbcl_pdd_3.5yr)
)

write.csv(cbcl_corrs, "cbcl_corrs.csv")



## >MAIN EFFECTS ##########################################################################################################


## >>Main effect of food security (3 level) ------------------------------------------------------------------


## On CBCL at 2 Years ##
foodsec_cbcl = rbind(
  runCatExpReg(exposure="food_cat_dqx", outcome=cbclvars_2yr),
  runCatExpReg(exposure="food_cat_i_1yr", outcome=cbclvars_2yr),
  runCatExpReg(exposure="food_cat_i_2yr", outcome=cbclvars_2yr)
) ##Running ME models

foodsec_cbcl$exposure = factor(foodsec_cbcl$exposure, 
                               levels = c("food_cat_dqx","food_cat_i_1yr", "food_cat_i_2yr" ),
                               labels = c("Delivery","Year 1", "Year 2" ))

foodsec_cbcl$outcome = factor(foodsec_cbcl$outcome,
                              levels = c("cbcl_agg_2yr", "cbcl_anx_dep_2yr", "cbcl_att_2yr", "cbcl_with_2yr", "cbcl_ext_2yr",
                                         "cbcl_adhd_2yr", "cbcl_aff_2yr", "cbcl_anx_2yr", "cbcl_odd_2yr", "cbcl_pdd_2yr"),
                              labels = c("Aggressive", "Anxious/Depressed", "Attention", "Withdrawn", "Externalizing",
                                         "ADHD Problems", "Affective", "Anxiety", "ODD Problems", "Pervisive Developmental\nProblems"))

saveMEPlot(gd=foodsec_cbcl, ylabel="Outcome Score\nAt 2 Years", plotname="foodsec_cbcl2")


## On CBCL at 3.5 Years ##
foodsec_cbcl3 = rbind(
  runCatExpReg(exposure="food_cat_dqx", outcome=cbclvars_3yr),
  runCatExpReg(exposure="food_cat_i_1yr", outcome=cbclvars_3yr),
  runCatExpReg(exposure="food_cat_i_2yr", outcome=cbclvars_3yr)
) ##Running ME models

foodsec_cbcl3$exposure = factor(foodsec_cbcl3$exposure, 
                               levels = c("food_cat_dqx","food_cat_i_1yr", "food_cat_i_2yr" ),
                               labels = c("Delivery","Year 1", "Year 2" ))

foodsec_cbcl3$outcome = factor(foodsec_cbcl3$outcome,
                              levels = c("cbcl_agg_3.5yr", "cbcl_anx_dep_3.5yr", "cbcl_att_3.5yr", "cbcl_with_3.5yr", "cbcl_ext_3.5yr",
                                         "cbcl_adhd_3.5yr", "cbcl_aff_3.5yr", "cbcl_anx_3.5yr", "cbcl_odd_3.5yr", "cbcl_pdd_3.5yr",
                                         "cbcl_er_3.5yr", "cbcl_som_3.5yr", "cbcl_sleep_3.5yr"),
                              labels = c("Aggressive", "Anxious/Depressed", "Attention", "Withdrawn", "Externalizing",
                                         "ADHD Problems", "Affective", "Anxiety", "ODD Problems", "Pervisive Developmental\nProblems",
                                         "Emotional Regulation", "Somatic Complaints", "Sleep Problems"))



saveMEPlot(gd=foodsec_cbcl3, ylabel="Outcome Score\nAt 3.5 Years", plotname="foodsec_cbcl3", height=12)




## On NIH Toolbox ##
foodsec_nih = rbind(
  runCatExpReg(exposure="food_cat_dqx", outcome=nihvars),
  runCatExpReg(exposure="food_cat_i_1yr", outcome=nihvars),
  runCatExpReg(exposure="food_cat_i_2yr", outcome=nihvars)
)

foodsec_nih$exposure = factor(foodsec_nih$exposure, 
                               levels = c("food_cat_dqx","food_cat_i_1yr", "food_cat_i_2yr" ),
                               labels = c("Delivery","Year 1", "Year 2" ))

foodsec_nih$outcome = factor(foodsec_nih$outcome,
                             levels=c("peg_fully.corrected.t.scores.dominant",
                                      "speed_fully.corrected.t.score",
                                      "attn_fully.corrected.t.score",
                                      "card_fully.corrected.t.score"),
                             labels = c("Dexterity",
                                        "Processing\nSpeed",
                                        "Inhibition\nand Attention",
                                        "Cognitive\nFlexibility"))

saveMEPlot(gd=foodsec_nih, ylabel="Outcome Score\nAt 5 Years", plotname="foodsec_nih5", height=5)



## On BSID at 1 year ##
foodsec_bsid1 = rbind(
  runCatExpReg(exposure="food_cat_dqx", outcome=bayleyvars_1yr),
  runCatExpReg(exposure="food_cat_i_1yr", outcome=bayleyvars_1yr),
  runCatExpReg(exposure="food_cat_i_2yr", outcome=bayleyvars_1yr)
)

foodsec_bsid1$exposure = factor(foodsec_bsid1$exposure, 
                              levels = c("food_cat_dqx","food_cat_i_1yr", "food_cat_i_2yr" ),
                              labels = c("Delivery","Year 1", "Year 2" ))

foodsec_bsid1$outcome = factor(foodsec_bsid1$outcome,
                               levels = c("cbcog_sc_1yr",
                                          "cbrc_sc_1yr",
                                          "cbec_sc_1yr",
                                          "cbfm_sc_1yr",
                                          "cbgm_sc_1yr",
                                          "cblang_comp_1yr",
                                          "cbmotor_comp_1yr",   
                                          "bayleyse_scaled_1yr"),
                               labels = c("Cognitive",
                                          "Receptive\ncommunication",
                                          "Expressive\ncommunication",
                                          "Fine\nmotor",
                                          "Gross\nmotor",
                                          "Language\ncomposite",
                                          "Motor\ncomposite",
                                          "Combined\nSocio-emotional" 
                               ))

saveMEPlot(gd=foodsec_bsid1, plotname="foodsec_bsid1", ylabel="Outcome Score\nAt 1 Year")



## On BSID at 2 Years ##
foodsec_bsid2 = rbind(
  runCatExpReg(exposure="food_cat_dqx", outcome=bayleyvars_2yr),
  runCatExpReg(exposure="food_cat_i_1yr", outcome=bayleyvars_2yr),
  runCatExpReg(exposure="food_cat_i_2yr", outcome=bayleyvars_2yr)
)

foodsec_bsid2$exposure = factor(foodsec_bsid2$exposure, 
                                levels = c("food_cat_dqx","food_cat_i_1yr", "food_cat_i_2yr" ),
                                labels = c("Delivery","Year 1", "Year 2" ))

foodsec_bsid2$outcome = factor(foodsec_bsid2$outcome,
                               levels = c("cbcog_sc_2yr",
                                          "cbrc_sc_2yr",
                                          "cbec_sc_2yr",
                                          "cbfm_sc_2yr",
                                          "cbgm_sc_2yr",
                                          "cblang_comp_2yr",
                                          "cbmotor_comp_2yr"),
                               labels = c("Cognitive",
                                          "Receptive\ncommunication",
                                          "Expressive\ncommunication",
                                          "Fine\nmotor",
                                          "Gross\nmotor",
                                          "Language\ncomposite",
                                          "Motor\ncomposite" 
                               ))

saveMEPlot(gd=foodsec_bsid2, plotname="foodsec_bsid2", ylabel="Outcome Score\nAt 2 Years")



## >>Main Effect of Poverty ------------------------------------------------------------------------------------


## On CBCL at 2 Years ##
pov_cbcl = rbind(
  runCatExpReg(exposure="poverty_food_i_dqx", outcome=cbclvars_2yr),
  runCatExpReg(exposure="poverty_food_1yr", outcome=cbclvars_2yr),
  runCatExpReg(exposure="poverty_food_2yr", outcome=cbclvars_2yr)
) ##Running ME models

pov_cbcl$exposure = factor(pov_cbcl$exposure, 
                               levels = c("poverty_food_i_dqx","poverty_food_1yr", "poverty_food_2yr" ),
                               labels = c("Delivery","Year 1", "Year 2" ))

pov_cbcl$outcome = factor(pov_cbcl$outcome,
                              levels = c("cbcl_agg_2yr", "cbcl_anx_dep_2yr", "cbcl_att_2yr", "cbcl_with_2yr", "cbcl_ext_2yr",
                                         "cbcl_adhd_2yr", "cbcl_aff_2yr", "cbcl_anx_2yr", "cbcl_odd_2yr", "cbcl_pdd_2yr"),
                              labels = c("Aggressive", "Anxious/Depressed", "Attention", "Withdrawn", "Externalizing",
                                         "ADHD Problems", "Affective", "Anxiety", "ODD Problems", "Pervisive Developmental\nProblems"))


saveMEPlot(gd=pov_cbcl, ylabel="Outcome Score\nAt 2 Years", plotname="pov_cbcl2")


## On CBCL at 3.5 Years ##
pov_cbcl3 = rbind(
  runCatExpReg(exposure="poverty_food_i_dqx", outcome=cbclvars_3yr),
  runCatExpReg(exposure="poverty_food_1yr", outcome=cbclvars_3yr),
  runCatExpReg(exposure="poverty_food_2yr", outcome=cbclvars_3yr)
) ##Running ME models

pov_cbcl3$exposure = factor(pov_cbcl3$exposure, 
                           levels = c("poverty_food_i_dqx","poverty_food_1yr", "poverty_food_2yr" ),
                           labels = c("Delivery","Year 1", "Year 2" ))

pov_cbcl3$outcome = factor(pov_cbcl3$outcome,
                           levels = c("cbcl_agg_3.5yr", "cbcl_anx_dep_3.5yr", "cbcl_att_3.5yr", "cbcl_with_3.5yr", "cbcl_ext_3.5yr",
                                      "cbcl_adhd_3.5yr", "cbcl_aff_3.5yr", "cbcl_anx_3.5yr", "cbcl_odd_3.5yr", "cbcl_pdd_3.5yr",
                                      "cbcl_er_3.5yr", "cbcl_som_3.5yr", "cbcl_sleep_3.5yr"),
                           labels = c("Aggressive", "Anxious/Depressed", "Attention", "Withdrawn", "Externalizing",
                                      "ADHD Problems", "Affective", "Anxiety", "ODD Problems", "Pervisive Developmental\nProblems",
                                      "Emotional Regulation", "Somatic Complaints", "Sleep Problems"))


saveMEPlot(gd=pov_cbcl3, ylabel="Outcome Score\nAt 3.5 Years", plotname="pov_cbcl3", height=11)



## On NIH Toolbox at 5 Years ##
pov_nih = rbind(
  runCatExpReg(exposure="poverty_food_i_dqx", outcome=nihvars),
  runCatExpReg(exposure="poverty_food_1yr", outcome=nihvars),
  runCatExpReg(exposure="poverty_food_2yr", outcome=nihvars)
) ##Running ME models

pov_nih$exposure = factor(pov_nih$exposure, 
                           levels = c("poverty_food_i_dqx","poverty_food_1yr", "poverty_food_2yr" ),
                           labels = c("Delivery","Year 1", "Year 2" ))


pov_nih$outcome = factor(pov_nih$outcome,
                             levels=c("peg_fully.corrected.t.scores.dominant",
                                      "speed_fully.corrected.t.score",
                                      "attn_fully.corrected.t.score",
                                      "card_fully.corrected.t.score"),
                             labels = c("Dexterity",
                                        "Processing\nSpeed",
                                        "Inhibition\nand Attention",
                                        "Cognitive\nFlexibility"))

saveMEPlot(gd=pov_nih, ylabel="Outcome Score\nAt 5 Years", plotname="pov_nih", height=5)


## On BSID at 1 year ##
pov_bsid1 = rbind(
  runCatExpReg(exposure="poverty_food_i_dqx", outcome=bayleyvars_1yr),
  runCatExpReg(exposure="poverty_food_1yr", outcome=bayleyvars_1yr),
  runCatExpReg(exposure="poverty_food_2yr", outcome=bayleyvars_1yr)
)

pov_bsid1$exposure = factor(pov_bsid1$exposure, 
                            levels = c("poverty_food_i_dqx","poverty_food_1yr", "poverty_food_2yr" ),
                            labels = c("Delivery","Year 1", "Year 2" ))

pov_bsid1$outcome = factor(pov_bsid1$outcome,
                               levels = c("cbcog_sc_1yr",
                                          "cbrc_sc_1yr",
                                          "cbec_sc_1yr",
                                          "cbfm_sc_1yr",
                                          "cbgm_sc_1yr",
                                          "cblang_comp_1yr",
                                          "cbmotor_comp_1yr",   
                                          "bayleyse_scaled_1yr"),
                               labels = c("Cognitive",
                                          "Receptive\ncommunication",
                                          "Expressive\ncommunication",
                                          "Fine\nmotor",
                                          "Gross\nmotor",
                                          "Language\ncomposite",
                                          "Motor\ncomposite",
                                          "Combined\nSocio-emotional" 
                               ))

saveMEPlot(gd=pov_bsid1, plotname="pov_bsid1", ylabel="Outcome Score\nAt 1 Year")


## On BSID at 2 year ##
pov_bsid2 = rbind(
  runCatExpReg(exposure="poverty_food_i_dqx", outcome=bayleyvars_2yr),
  runCatExpReg(exposure="poverty_food_1yr", outcome=bayleyvars_2yr),
  runCatExpReg(exposure="poverty_food_2yr", outcome=bayleyvars_2yr)
)

pov_bsid2$exposure = factor(pov_bsid2$exposure, 
                            levels = c("poverty_food_i_dqx","poverty_food_1yr", "poverty_food_2yr" ),
                            labels = c("Delivery","Year 1", "Year 2" ))

pov_bsid2$outcome = factor(pov_bsid2$outcome,
                           levels = c("cbcog_sc_2yr",
                                      "cbrc_sc_2yr",
                                      "cbec_sc_2yr",
                                      "cbfm_sc_2yr",
                                      "cbgm_sc_2yr",
                                      "cblang_comp_2yr",
                                      "cbmotor_comp_2yr"),
                           labels = c("Cognitive",
                                      "Receptive\ncommunication",
                                      "Expressive\ncommunication",
                                      "Fine\nmotor",
                                      "Gross\nmotor",
                                      "Language\ncomposite",
                                      "Motor\ncomposite" 
                           ))

saveMEPlot(gd=pov_bsid2, plotname="pov_bsid2", ylabel="Outcome Score\nAt 1 Year")



## >>Main Effect of Wealth Index ------------------------------------------------------------------------------------------

wealth_alloutcomes = rbind(
  runContExpReg(exposure="wealth_i", outcome=cbclvars),
  runContExpReg(exposure="wealth_i", outcome=bayleyvars),
  runContExpReg(exposure="wealth_i", outcome=nihvars)
)

write.csv(wealth_alloutcomes, "wealth_mes.csv")


## >>NEED TO UPDATE AND RUN WITH CBCL 2/3.5 YEAR -- Main Effect of 2-Level Food Security ----------------------------------------------------------------------------------

dr$insecure_low_dqx = as.factor(dr$insecure_low_dqx) #need to convert to factors to run models using functions I created
dr$insecure_low_i_1yr = as.factor(dr$insecure_low_i_1yr)
dr$insecure_low_i_2yr = as.factor(dr$insecure_low_i_2yr)

## On CBCL at 2 Years ##
foodsec2_cbcl = rbind(
  runCatExpReg(exposure="insecure_low_dqx", outcome=cbclvars),
  runCatExpReg(exposure="insecure_low_i_1yr", outcome=cbclvars),
  runCatExpReg(exposure="insecure_low_i_2yr", outcome=cbclvars)
) ##Running ME models

foodsec2_cbcl$exposure = factor(foodsec2_cbcl$exposure, 
                               levels = c("insecure_low_dqx","insecure_low_i_1yr", "insecure_low_i_2yr" ),
                               labels = c("Delivery","Year 1", "Year 2" ))

foodsec2_cbcl$outcome = factor(foodsec2_cbcl$outcome,
                              levels = c("cbcl_agg_2yr", "cbcl_anx_dep_2yr", "cbcl_att_2yr", "cbcl_with_2yr", "cbcl_ext_2yr",
                                         "cbcl_adhd_2yr", "cbcl_aff_2yr", "cbcl_anx_2yr", "cbcl_odd_2yr", "cbcl_pdd_2yr"),
                              labels = c("Aggressive", "Anxious/Depressed", "Attention", "Withdrawn", "Externalizing",
                                         "ADHD Problems", "Affective", "Anxiety", "ODD Problems", "Pervisive Developmental\nProblems"))


saveMEPlot(gd=foodsec2_cbcl, ylabel="Outcome Score\nAt 2 Years", plotname="foodsec2_cbcl2")


## On NIH Toolbox ##
foodsec2_nih = rbind(
  runCatExpReg(exposure="insecure_low_dqx", outcome=nihvars),
  runCatExpReg(exposure="insecure_low_i_1yr", outcome=nihvars),
  runCatExpReg(exposure="insecure_low_i_2yr", outcome=nihvars)
)

foodsec2_nih$exposure = factor(foodsec2_nih$exposure, 
                              levels = c("insecure_low_dqx","insecure_low_i_1yr", "insecure_low_i_2yr" ),
                              labels = c("Delivery","Year 1", "Year 2" ))

foodsec2_nih$outcome = factor(foodsec2_nih$outcome,
                             levels=c("peg_fully.corrected.t.scores.dominant",
                                      "speed_fully.corrected.t.score",
                                      "attn_fully.corrected.t.score",
                                      "card_fully.corrected.t.score"),
                             labels = c("Dexterity",
                                        "Processing\nSpeed",
                                        "Inhibition\nand Attention",
                                        "Cognitive\nFlexibility"))

saveMEPlot(gd=foodsec2_nih, ylabel="Outcome Score\nAt 5 Years", plotname="foodsec2_nih5", height=5)



## On BSID at 1 year ##
foodsec2_bsid1 = rbind(
  runCatExpReg(exposure="insecure_low_dqx", outcome=bayleyvars_1yr),
  runCatExpReg(exposure="insecure_low_i_1yr", outcome=bayleyvars_1yr),
  runCatExpReg(exposure="insecure_low_i_2yr", outcome=bayleyvars_1yr)
)

foodsec_bsid1$exposure = factor(foodsec2_bsid1$exposure, 
                                levels = c("insecure_low_dqx","insecure_low_i_1yr", "insecure_low_i_2yr" ),
                                labels = c("Delivery","Year 1", "Year 2" ))

foodsec2_bsid1$outcome = factor(foodsec2_bsid1$outcome,
                               levels = c("cbcog_sc_1yr",
                                          "cbrc_sc_1yr",
                                          "cbec_sc_1yr",
                                          "cbfm_sc_1yr",
                                          "cbgm_sc_1yr",
                                          "cblang_comp_1yr",
                                          "cbmotor_comp_1yr",   
                                          "bayleyse_scaled_1yr"),
                               labels = c("Cognitive",
                                          "Receptive\ncommunication",
                                          "Expressive\ncommunication",
                                          "Fine\nmotor",
                                          "Gross\nmotor",
                                          "Language\ncomposite",
                                          "Motor\ncomposite",
                                          "Combined\nSocio-emotional" 
                               ))

saveMEPlot(gd=foodsec2_bsid1, plotname="foodsec2_bsid1", ylabel="Outcome Score\nAt 1 Year")



## On BSID at 2 Years ##
foodsec2_bsid2 = rbind(
  runCatExpReg(exposure="insecure_low_dqx", outcome=bayleyvars_2yr),
  runCatExpReg(exposure="insecure_low_i_1yr", outcome=bayleyvars_2yr),
  runCatExpReg(exposure="insecure_low_i_2yr", outcome=bayleyvars_2yr)
)

foodsec2_bsid2$exposure = factor(foodsec2_bsid2$exposure, 
                                levels = c("insecure_low_dqx","insecure_low_i_1yr", "insecure_low_i_2yr" ),
                                labels = c("Delivery","Year 1", "Year 2" ))

foodsec2_bsid2$outcome = factor(foodsec2_bsid2$outcome,
                               levels = c("cbcog_sc_2yr",
                                          "cbrc_sc_2yr",
                                          "cbec_sc_2yr",
                                          "cbfm_sc_2yr",
                                          "cbgm_sc_2yr",
                                          "cblang_comp_2yr",
                                          "cbmotor_comp_2yr"),
                               labels = c("Cognitive",
                                          "Receptive\ncommunication",
                                          "Expressive\ncommunication",
                                          "Fine\nmotor",
                                          "Gross\nmotor",
                                          "Language\ncomposite",
                                          "Motor\ncomposite" 
                               ))

saveMEPlot(gd=foodsec2_bsid2, plotname="foodsec2_bsid2", ylabel="Outcome Score\nAt 2 Years")




## >4-WAY DECOMPOSITION MEDIATION USING HARD CODE, INCLUDING PLOTS - 2019-08-26 Version ###############################################



## >>Poverty-->Security--> Outcomes (with fancy plot..) -------------------------------------------------------------------------

## Note that food security 2-level was used because of the way mediation is computed...
      #BUT, main effects using the 2-level food security are very similar 

# Now only loo=oking at mediation at the same timepoint 

## Poverty @ 0, Security @ 0 ##
iloop = list()
iloop2 = list()
for (i in seq_along(cbclvars)) {
  
  medresults = runMediationBoot(inputdata=dr, exposure="poverty_food_i_dqx", mediator="insecure_low_dqx", outcome=cbclvars[i],
                                adjustment=adj2, index_exp = "At or Below", m_forcde=0, nsim=1000) #INCREASE TO AT LEAST 5000 ONCE CODE SET
  
  iloop[i] = list(cleanMediationResults(medresults)) #to create formatted tables for writign to csv
  
  iloop2[i] = list(cleanMediationResultsPlot(medresults)) #to create a dataframe for plotting 
  
}
write.csv(do.call("rbind", iloop), "pov0_sec0_4way.csv")
pov0sec0_plot = do.call("rbind", iloop2)
rm(iloop, iloop2)


## Poverty @ 1, Security @ 1 ##
iloop = list()
iloop2 = list()
for (i in seq_along(cbclvars)) {
  
  medresults = runMediationBoot(inputdata=dr, exposure="poverty_food_1yr", mediator="insecure_low_i_1yr", outcome=cbclvars[i],
                                adjustment=adj2, index_exp = "At or Below", m_forcde=0, nsim=1000) #INCREASE TO AT LEAST 5000 ONCE CODE SET
  
  iloop[i] = list(cleanMediationResults(medresults)) #to create formatted tables for writign to csv
  
  iloop2[i] = list(cleanMediationResultsPlot(medresults)) #to create a dataframe for plotting 
  
}
write.csv(do.call("rbind", iloop), "pov1_sec1_4way.csv")
pov1sec1_plot = do.call("rbind", iloop2)
rm(iloop, iloop2)



## Poverty @ 2, Security @ 2 ##
iloop = list()
iloop2 = list()
for (i in seq_along(cbclvars)) {
  
  medresults = runMediationBoot(inputdata=dr, exposure="poverty_food_2yr", mediator="insecure_low_i_2yr", outcome=cbclvars[i],
                                adjustment=adj2, index_exp = "At or Below", m_forcde=0, nsim=1000) #INCREASE TO AT LEAST 5000 ONCE CODE SET
  
  iloop[i] = list(cleanMediationResults(medresults)) #to create formatted tables for writign to csv
  
  iloop2[i] = list(cleanMediationResultsPlot(medresults)) #to create a dataframe for plotting 
  
}
write.csv(do.call("rbind", iloop), "pov2_sec2_4way.csv")
pov2sec2_plot = do.call("rbind", iloop2)
rm(iloop, iloop2)


## Poverty @ 2, Security @ 2  - CRUDE ##
iloop = list()
iloop2 = list()
for (i in seq_along(cbclvars)) {
  
  medresults = runMediationBootCrude(inputdata=dr, exposure="poverty_food_2yr", mediator="insecure_low_i_2yr", outcome=cbclvars[i],
                                index_exp = "At or Below", m_forcde=0, nsim=1000) #INCREASE TO AT LEAST 5000 ONCE CODE SET
  
  iloop[i] = list(cleanMediationResults(medresults)) #to create formatted tables for writign to csv
  
  iloop2[i] = list(cleanMediationResultsPlot(medresults)) #to create a dataframe for plotting 
  
}
write.csv(do.call("rbind", iloop), "pov2_sec2_4way_Crude.csv")
pov2sec2_plot = do.call("rbind", iloop2)
rm(iloop, iloop2)




## Plotting Concurrent timepoints ##
pov0sec0_plot$model = "Exposure: Poverty at Delivery\nMediator: Food Security at Delivery" #creating a factor for each model
pov1sec1_plot$model = "Exposure: Poverty at 1 Year\nMediator: Food Security at 1 Year"
pov2sec2_plot$model = "Exposure: Poverty at 2 Years\nMediator: Food Security at 2 Years"

povsec_plotc = rbind(pov0sec0_plot, pov1sec1_plot, pov2sec2_plot)
povsec_plotc$model = as.factor(povsec_plotc$model)


povsec_plotc$Outcome = factor(povsec_plotc$Outcome,
                             levels = c("cbcl_agg_2yr", "cbcl_anx_dep_2yr", "cbcl_att_2yr", "cbcl_with_2yr", "cbcl_ext_2yr",
                                        "cbcl_adhd_2yr", "cbcl_aff_2yr", "cbcl_anx_2yr", "cbcl_odd_2yr", "cbcl_pdd_2yr"),
                             labels = c("Aggressive", "Anxious/\nDepressed", "Attention", "Withdrawn", "Externalizing",
                                        "ADHD\nProblems", "Affective", "Anxiety", "ODD\nProblems", "Pervisive\nDevelopmental\nProblems"))

## Make Plot ##
makeMedPlot(df=povsec_plotc, name="povsec_cbcl2_concurrent")




## >>Security-->Nutrition--> Outcomes (with fancy plot..) -------------------------------------------------------------------------

## Note that food security 2-level was used because of the way mediation is computed...
#BUT, main effects using the 2-level food security are very similar 

# Note: only looking at associations with 1 and 2 year food security since that was the only place there were effects



## Security @ 1 year, Diet Diversity @ 1 years ## 
iloop = list()
iloop2 = list()
for (i in seq_along(cbclvars)) {
  
  medresults = runMediationBootContM(inputdata=dr, exposure="insecure_low_i_1yr", mediator="dietdiversity_yr1", outcome=cbclvars[i],
                                     adjustment=adj2, index_exp = 1, m_forcde=mean(dr$dietdiversity_yr2, na.rm=T), nsim=1000) #INCREASE TO AT LEAST 5000 ONCE CODE SET
  
  iloop[i] = list(cleanMediationResults(medresults)) #to create formatted tables for writign to csv
  
  iloop2[i] = list(cleanMediationResultsPlot(medresults)) #to create a dataframe for plotting 
  
}
write.csv(do.call("rbind", iloop), "sec1_dd1_4way.csv")
sec1dd1_plot = do.call("rbind", iloop2)
rm(iloop, iloop2)


## Security @ 1 year, Diet Diversity @ 2 years ## 
iloop = list()
iloop2 = list()
for (i in seq_along(cbclvars)) {
  
  medresults = runMediationBootContM(inputdata=dr, exposure="insecure_low_i_1yr", mediator="dietdiversity_yr2", outcome=cbclvars[i],
                                adjustment=adj2, index_exp = 1, m_forcde=mean(dr$dietdiversity_yr2, na.rm=T), nsim=1000) #INCREASE TO AT LEAST 5000 ONCE CODE SET
  
  iloop[i] = list(cleanMediationResults(medresults)) #to create formatted tables for writign to csv
  
  iloop2[i] = list(cleanMediationResultsPlot(medresults)) #to create a dataframe for plotting 
  
}
write.csv(do.call("rbind", iloop), "sec1_dd2_4way.csv")
sec1dd2_plot = do.call("rbind", iloop2)
rm(iloop, iloop2)


## Security @ 2 year, Diet Diversity @ 2 years ## 
iloop = list()
iloop2 = list()
for (i in seq_along(cbclvars)) {
  
  medresults = runMediationBootContM(inputdata=dr, exposure="insecure_low_i_2yr", mediator="dietdiversity_yr2", outcome=cbclvars[i],
                                     adjustment=adj2, index_exp = 1, m_forcde=mean(dr$dietdiversity_yr2, na.rm=T), nsim=1000) #INCREASE TO AT LEAST 5000 ONCE CODE SET
  
  iloop[i] = list(cleanMediationResults(medresults)) #to create formatted tables for writign to csv
  
  iloop2[i] = list(cleanMediationResultsPlot(medresults)) #to create a dataframe for plotting 
  
}
write.csv(do.call("rbind", iloop), "sec2_dd2_4way.csv")
sec2dd2_plot = do.call("rbind", iloop2)
rm(iloop, iloop2)



## Plotting ##
sec1dd1_plot$model = "Exposure: Food Security at 1 Year\nMediator: Diet Diversity at 1 Years"
sec1dd2_plot$model = "Exposure: Food Security at 1 Year\nMediator: Diet Diversity at 2 Years"
sec2dd2_plot$model = "Exposure: Food Security at 2 Year\nMediator: Diet Diversity at 2 Years"

secdd_plot = rbind(sec1dd1_plot, sec1dd2_plot, sec2dd2_plot)

secdd_plot$model = as.factor(secdd_plot$model)

secdd_plot$Outcome = factor(secdd_plot$Outcome,
                             levels = c("cbcl_agg_2yr", "cbcl_anx_dep_2yr", "cbcl_att_2yr", "cbcl_with_2yr", "cbcl_ext_2yr",
                                        "cbcl_adhd_2yr", "cbcl_aff_2yr", "cbcl_anx_2yr", "cbcl_odd_2yr", "cbcl_pdd_2yr"),
                             labels = c("Aggressive", "Anxious/\nDepressed", "Attention", "Withdrawn", "Externalizing",
                                        "ADHD\nProblems", "Affective", "Anxiety", "ODD\nProblems", "Pervisive\nDevelopmental\nProblems"))

## Make Plot ##
makeMedPlot(df=secdd_plot, name="secdd_cbcl")


## >>Seurity-->Depression--> Outcomes (with fancy plot..) -------------------------------------------------------------------------

## Note that food security 2-level was used because of the way mediation is computed...
#BUT, main effects using the 2-level food security are very similar 

dr$srq6_i_1yr = as.factor(dr$srq6_i_1yr) #making depression a factor for code
dr$srq6_2yr = as.factor(dr$srq6_2yr) 
levels(dr$srq6_i_1yr)
levels(dr$srq6_2yr)


## Security @ 1, Depression @ 1 ##
iloop = list()
iloop2 = list()
for (i in seq_along(cbclvars)) {
  
  medresults = runMediationBoot(inputdata=dr, "insecure_low_i_1yr", mediator="srq6_i_1yr", outcome=cbclvars[i],
                                adjustment=adj2, index_exp = 1, m_forcde=0, nsim=1000) #INCREASE TO AT LEAST 5000 ONCE CODE SET
  
  iloop[i] = list(cleanMediationResults(medresults)) #to create formatted tables for writign to csv
  
  iloop2[i] = list(cleanMediationResultsPlot(medresults)) #to create a dataframe for plotting 
  
}
write.csv(do.call("rbind", iloop), "sec1_dep1_4way.csv")
sec1dep1_plot = do.call("rbind", iloop2)
rm(iloop, iloop2)



## Security @ 1, Depression @ 2 ##
iloop = list()
iloop2 = list()
for (i in seq_along(cbclvars)) {
  
  medresults = runMediationBoot(inputdata=dr, "insecure_low_i_1yr", mediator="srq6_2yr", outcome=cbclvars[i],
                                adjustment=adj2, index_exp = 1, m_forcde=0, nsim=1000) #INCREASE TO AT LEAST 5000 ONCE CODE SET
  
  iloop[i] = list(cleanMediationResults(medresults)) #to create formatted tables for writign to csv
  
  iloop2[i] = list(cleanMediationResultsPlot(medresults)) #to create a dataframe for plotting 
  
}
write.csv(do.call("rbind", iloop), "sec1_dep2_4way.csv")
sec1dep2_plot = do.call("rbind", iloop2)
rm(iloop, iloop2)


## Security @ 2, Depression @ 2 ##
iloop = list()
iloop2 = list()
for (i in seq_along(cbclvars)) {
  
  medresults = runMediationBoot(inputdata=dr, "insecure_low_i_2yr", mediator="srq6_2yr", outcome=cbclvars[i],
                                adjustment=adj2, index_exp = 1, m_forcde=0, nsim=1000) #INCREASE TO AT LEAST 5000 ONCE CODE SET
  
  iloop[i] = list(cleanMediationResults(medresults)) #to create formatted tables for writign to csv
  
  iloop2[i] = list(cleanMediationResultsPlot(medresults)) #to create a dataframe for plotting 
  
}
write.csv(do.call("rbind", iloop), "sec2_dep2_4way.csv")
sec2dep2_plot = do.call("rbind", iloop2)
rm(iloop, iloop2)





## Plotting ##
sec1dep1_plot$model = "Exposure: Food Security at 1 Year\nMediator: Maternal Depression at 1 Years"
sec1dep2_plot$model = "Exposure: Food Security at 1 Year\nMediator: Maternal Depression at 2 Years"
sec2dep2_plot$model = "Exposure: Food Security at 2 Year\nMediator: Maternal Depression at 2 Years"

secdep_plot = rbind(sec1dep1_plot, sec1dep2_plot, sec2dep2_plot)

secdep_plot$model = as.factor(secdep_plot$model)
levels(secdep_plot$model)

secdep_plot$Outcome = factor(secdep_plot$Outcome,
                              levels = c("cbcl_agg_2yr", "cbcl_anx_dep_2yr", "cbcl_att_2yr", "cbcl_with_2yr", "cbcl_ext_2yr",
                                         "cbcl_adhd_2yr", "cbcl_aff_2yr", "cbcl_anx_2yr", "cbcl_odd_2yr", "cbcl_pdd_2yr"),
                              labels = c("Aggressive", "Anxious/\nDepressed", "Attention", "Withdrawn", "Externalizing",
                                         "ADHD\nProblems", "Affective", "Anxiety", "ODD\nProblems", "Pervisive\nDevelopmental\nProblems"))


## Make Plot ##
makeMedPlot(df=secdep_plot, name="secdep_cbcl")









## >SOME REGRESSION MODELS FOR COMPARING WITH STATA ########################################################### 


## Running Distributions to Check Against STATA ##
table(dr$edu_cat, useNA="always")
table(dr$married, useNA="always")
summary(dr$mage)
table(dr$parity2, useNA="always")
summary(dr$cbcl_anx_dep_2yr, useNA="always")
table(dr$food_cat_i_1yr, useNA="always")
table(dr$poverty_food_1yr, useNA="always")
table(dr$insecure_low_i_1yr, useNA="always")
summary(dr$cbcl_agg_2yr, useNA="always")


## Main Effect Regressions ##

crude1 = lm(data=dr, 
   formula =
     cbcl_anx_dep_2yr ~ food_cat_i_1yr)
summary(crude1)


me1 = lm(data=dr,
   formula = 
     cbcl_anx_dep_2yr ~ food_cat_i_1yr + edu_cat + married + mage + parity2) 
summary(me1)


me2 = lm(data=dr,
         formula = 
           cbcl_anx_dep_2yr ~ food_cat_i_1yr + poverty_food_1yr + edu_cat + married + mage + parity2 )
summary(me2)


me3 = lm(data=dr,
         formula = 
           cbcl_anx_dep_2yr ~ insecure_low_i_1yr + poverty_food_1yr + edu_cat + married + mage + parity2 )
summary(me3)



## Mediation ##

cleanMediationResults(runMediationBoot(inputdata=dr, exposure="poverty_food_1yr", mediator="insecure_low_i_1yr", outcome="cbcl_anx_dep_2yr",
                                       adjustment=adj2, index_exp="At or Below", m_forcde=0, nsim=1000)
)


cleanMediationResults(runMediationBoot(inputdata=dr, exposure="poverty_food_1yr", mediator="insecure_low_i_1yr", outcome="cbcl_agg_2yr",
                                       adjustment=adj2, index_exp="At or Below", m_forcde=0, nsim=1000))








