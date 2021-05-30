#This code shows example of how model objects were plotted in Fig.3 and ID-TaxER tool

#script assumes you have the following variables in working environment
#1) mod containing HOF model for single OTU

#Note you can also simply use plot(mod), the below code was used to ensure that majority of OTU abundance data is clearly visible.
#by setting max ylim value (OTU read number) to probabilty 0.999 in quantile function. 
#this prevents y axis from being much larger than necessary due to outliers.
#===========================================================================================================================================================================================
#load in eHOF library
library(eHOF)
#===========================================================================================================================================================================================

#get model choice
mod_choice<-Para(mod)$model
#get fit for best model
fitted=mod$models[mod_choice][[1]][9]
mod_stats=as.data.frame(cbind(mod$y,mod$x,fitted$fitted))
#order by ph
mod_stats=mod_stats[order(mod_stats$V2),]

#plot model
plot(mod_stats$V2,mod_stats$V1,ylim=c(0,quantile(mod_stats$V1,0.999)),cex=0.5,pch=19,xlab="pH",ylab="Number of Reads",col="#2B2D2F",cex.lab=1.5,cex.axis=1.4)
lines(mod_stats$V2,mod_stats$V3,ylim=c(0,quantile(mod_stats$V1,0.999)),cex=0.07,lwd=1.2,type="o")

