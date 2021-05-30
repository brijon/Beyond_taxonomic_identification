#The below code shows how to model a single OTU's pH response with HOF
#and then use this model to assign OTU pH classification based on HOF model optima


#script assumes you have the following variables in working environment
#1) otu_table data frame/ matrix with otus as colnames and samples as rownames 

#2) a variable called sample_ph with pH of samples, with the same sample order as in otu_table.

#3) otu variable specifying OTU of interest
#===========================================================================================================================================================================================
#load in eHOF library
library(eHOF)
#===========================================================================================================================================================================================

#run HOF with poisson family , assessing model fit with AIC and bootstrapping using model shapes I:V
mod=HOF(otu_table[,otu],family=poisson, sample_pH,test="AIC",model=c("I","II","III","IV","V"))

#get parameters of model including model optima
param=Para(mod)

#simple function for assigning pH class from single model optima.
ph_class=function(optima){
  if(optima<=5.2){
    classification="Acid"
  }
  if(optima<7 && optima>5.2){
    classification="Mid"
  }
  
  if(optima>=7){
    classification="Neutral"
  }
  return(classification)
}


#function to assign pH class taking into account whether there are multiple optima
ph_class_multi_optima=function(param){
#if only one optima
  if(length(param$opt)==1){
#and this optima is not NA (in the case of model II,IV and V)
    if(!is.na(param$opt[[1]])){
#get pH class for optima      
      return(ph_class(param$opt[[1]]))
#if no optima  (in the case of model I) return none as class cant be assigned     
    }else(return("none"))}
#if two optima (in the case of model III)
  else{
#if ph class for both optima are the same
    if(ph_class(param$opt[[1]])==ph_class(param$opt[[2]])){
#return ph class for first optima      
      return(ph_class(param$opt[[1]]))
    }
#if the two optima are not assigned the same pH class return classes for both optima
    else{ return(paste(ph_class(param$opt[[1]]),"to",ph_class(param$opt[[2]]),sep=" "))}
  }
}


#run function on param object
ph_class_multi_optima(param)