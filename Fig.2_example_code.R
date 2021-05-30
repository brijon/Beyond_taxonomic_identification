#Example script for Fig.2
#code to plot query OTU rank abundance class vs proportion of query OTU's with match to reference database
#script aims to gain insights in the dominance of query OTUS not found in the reference database
#===========================================================================================================================================================================================
#load library
library(dplyr)
#===========================================================================================================================================================================================
#script assumes you have the following variables in working environment
#1)  query_hits dataframe, with each row representing individual Query OTU
#rel_abund collumn should contain query OTU relative abundance within query dataset
#hit collumn  should  denote whether OTU has hit in the reference database (based on 97% blast allignment), represented by binary yes/no values.

#add rank abundance quantile based on query OTU relative abundance 
#1000 representing the most abundant rank abundance class, 1 being the least
query_hits$quantile <- ntile(query_hits$rel_abund, 1000)  

#get table with number of yes/no  values in hit collumn per abundance quantile
query_hits_table<-table(as.factor(query_hits$quantile),query_hits$hit)
query_hits_table<-data.frame(query_hits_table)
#subset table to just get proportion of OTU's with hit within abundance quantile
pos_query_hits_table<-query_hits_table[which(query_hits_table$Var2=="yes"),]

#calculate proportion of OTU's with hit in reference database within each abundance quantile
pos_query_hits_table$prophits=pos_query_hits_table$Freq/(nrow(query_hits)/1000)

#plot proportion of matches to CS database vs rank abunsance class of query OTU
par(mar=c(5,3,2,2), pty = "s")
par(mgp=c(2.5,0.5,0))
plot(rev(pos_query_hits_table$prophits), xlab="Rank Abundance class of query OTU", pch=16,ylab="Proportion of matches to CS database",cex=0.7)
lines(lowess(rev(pos_query_hits_table$prophits)), col="red")
