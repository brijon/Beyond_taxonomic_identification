#Example code for Fig.1 
#Script to generate species accumulation curve for all samples and subsetted by habitat type

#script assumes you have the following variables in working environment
#1) otu_table data frame/ matrix with otus as colnames and samples as rownames 

#2) a variable called per_sample_habitat with sample habitat per sample, with the same sample order as in otu_table.

#3) a list called unique_habitats_list containing unique habitat values that you want to subset the data by 
#this could be created by running unique_habitats_list=unique(per_sample_habitat)
#===========================================================================================================================================================================================
#load libraries
library(vegan)
library(RColorBrewer)
#===========================================================================================================================================================================================

#run species accumulation function on all data
specaccum.all<-specaccum(otu_table, method = "random", permutations = 10)

#make list to add  specaccum.all variable to and then subsequent species accumulation curves to based on subsets of data (subset by habitat) 
specaccum.list=list()
specaccum.list[["specaccum.all"]]=specaccum.all
for(habitat in unique_habitats_list){
  
  #make species accumulation curve for habitat subset and add it to list with suitable names
  specaccum.list[[paste0("specaccum.",habitat)]]=specaccum(subset(otu_table,per_sample_habitat==habitat), method="random", permutations = 10)
}


#plot species accumulation curves
par(mar=c(5,4,2,2))
par(mgp=c(3,0.5,0))


ylm=c(0,max(specaccum.all$richness)+max(specaccum.all$sd))

sapply(seq_along(specaccum.list), function(i) {
  if (i==1) { # If it's the first list element, use plot()
    with(specaccum.list[[i]], {
      plot(sites, richness, type='l',ylim=ylm, 
           xlab='Sites', ylab='Richness', las=1,cex=0.7)
      segments(seq_len(max(sites)), y0=richness - 2*sd, 
               y1=richness + 2*sd)
    })    
  } else {
    with(specaccum.list[[i]], { # for subsequent elements, use lines()
      lines(sites, richness, col=i)
      segments(seq_len(max(sites)), y0=richness - 2*sd, 
               y1=richness + 2*sd, col=i)
    })     
  }
})

legend('bottomright', c('All sites',unique_habitats_list), col=1:8, lty=1,lwd=2.5, 
       bty='n', inset=0.025,cex=1)

