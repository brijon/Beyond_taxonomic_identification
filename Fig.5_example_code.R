# Example script for Fig.5

# This figure demonstrates how HOF models built on a large scale dataset quantifying
# the responses of specific OTUs to soil pH ("train" dataset), can be used to predict community structure in
# an independent "test" dataset using only the pH data for the test dataset.

# The following datasets are required:

#script assumes you have the following variables in working environment

# train datasets: here the large scale Countryside survey data.

# train_otu : the OTU table for the large scale "train" dataset (in the paper the GB Countryside Survey data)
# train_pH : pH values for each sample in train_otu, used to buld HOF models predicting abundances of OTUs in train_otu

# train_model_features: a data frame containing the pH optimum of each otu as output from the HOF modelling. Includes two 
# columns: train_otu_name and train_OTU_pH_optimum. See XXXXXX for how this was created

#===========================================================================================================================================================================================

# test datasets : here an independent set of sample from the UGRASS soil survey

# test_pH: pH data for an independent dataset, which will be used to predict OTU abundances
# test_otu: an OTU table for the same samples as test_pH, which will be used to validate predictions

# train_test_otu_lookup: A data frame identifying matches between the train dataset OTUs and test dataset OTUs. 
# 				Obtained following blast searches of the representative sequences for the test dataset against the representative 
# 				sequences of the train datasets - as described in XXXXXXXXXX.


#===========================================================================================================================================================================================
#load in libraries
library(vegan)
library(labdsv)
library(eHOF)
#===========================================================================================================================================================================================

#Figure 5a
#Ordination plot of the test dataset showing pH associations in the data

#first generate a grouping factor defining broad pH groups
phgroup<-ifelse(test_pH<5.2,"Acid","Mid")
phgroup<-ifelse(test_pH$pH>7,"Neutral",env1$phgroup)

#plot ordination of test data by pH group
mod<-metaMDS(test_otu)
plot(mod, display="sites")
ordispider(mod, phgroup, label=T)

#===========================================================================================================================================================================================

#Figure 5b

# xy plot showing correspondence between species indicators of the three pH groups in the test datasets, and the
# pH optimum of these OTUs determined after sequence matching to OTUs in the train dataset.

#1. Determine OTU indicators of pH groups in the test data set

inds<-indval(test_otu[,colSums(test_otu)>0],phgroup)
indies<-data.frame("Acid_indval"=inds$indval$Acid,"Mid_indval"=inds$indval$Mid,
"Neutral_indval"=inds$indval$Neutral,"p"=inds$pval,test_OTU_ID=names(inds$pval))

#manipulate indies dataframe to categorise test data OTUs according to whether an indicator of acid/mid/neutral pH
indies<-indies[indies$p<0.05,]#remove non significant indicators
indies$max<-apply(indies[, 1:3], 1, max) #identify  max indval value across the three ph groups
indies$OTUclass<-colnames(indies)[apply(indies,1,which.max)] # create column containing the ph group with max indval
indies$OTUclass<-gsub('_indval', '', indies$OTUclass) # strip "_indval"
indies<-indies[indies$max>0.5,] #subset table to only include strong indicators (with indval value>0.5)


#2. Create lookup dataframe with predicted pH optimum of test OTUs 

#requires:
#train_model_features: Output from HOF models. Contains "train_OTU_ID", "train_pH_optimum" and derived "pH.Class" columns columns
#train_test_otu_lookup: Output from blasting test repseqs v train repseqs. 
				#contains "train_OTU_ID", test_OTU_ID, and "identity" (percentage match) columns


lookup<-merge(train_model_features,train_test_otu_lookup,by.x="train_OTU_ID", by.y="train_OTU_ID")

#3. Add pH optimum of matched OTU from train dataset to test indicators dataframe

indies$train_pH_class<-with(lookup,
                     pH.Class[match(indies$test_OTU_ID,
                                       test_OTU_ID)])
indies$identity<-with(lookup,
                     identity[match(indies$test_OTU_ID,
                                       test_OTU_ID)])
indies$ph_optimum<-with(lookup,
                     pH.Optimum.1[match(indies$test_OTU_ID,
                                      test_OTU_ID)])

indies$CS_class<-ifelse(indies$identity<97,"Weak match (<97%)",as.character(indies$CS_class))


#4. Plot using ggplot

# create dataframe for plotting in ggplot
 
dat<-data.frame("train_pH_optimum"=indies$ph_optimum,"train_pH_class"=indies$train_pH_class,"test_pH_class"=indies$OTUclass)
dat1<-dat[dat$train_pH_class=="Acid"|dat$train_pH_class=="Mid"|dat$train_pH_class=="Neutral",]#subset to only include classified OTUS
dat1<-dat1[rowSums(is.na(dat1)) != ncol(dat1), ]#remove other NA containing rows

# plot

ggplot() + 
geom_point(data=dat1,size=0.8,position = position_jitterdodge(jitter.width = 0.9, dodge.width = 0.4),aes(x=test_pH_class, y=train_pH_optimum, color=train_pH_class))+
ggtitle("b) Validating OTU pH predictions")+
theme(legend.position="top")+
labs(x = "Observed pH class", y="Predicted pH optimum\n  ",color="Predicted\npH class" )+
guides(colour = guide_legend(override.aes = list(size=2)))
theme_update(text = element_text(size=10),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
strip.background = element_blank()
)


#===========================================================================================================================================================================================

#Figure 5c

# xy plot showing for the test dataset, the fit between observed community structure (first axis NMDS score)
# and predicted community structure based on predicted abundances of the top 100 most abundant test data OTUS.
# The abundances of these OTUs are predicted from the test dataset sample pH, and HOF model for for the matched OTU 
# in the train dataset.

#1. Select the 100 most abundant OTUs in test dataset

abunds.test<-data.frame(sort(colSums(test_otu),decreasing=TRUE)[1:100])

#2. Determine matches in train dataset

lookup<-merge(train_model_features,train_test_otu_lookup,by.x="train_OTU_ID", by.y="train_OTU_ID")

abunds.train<-lookup[lookup$test_OTU_ID %in%row.names(abunds.test) ,]

train_select_otus<-abunds.train$train_OTU_ID #select only otus in train which are matched to top 100 otus in test

#3. Run HOF models on only selected OTUs in the train dataset

train_otu_subset<-train_otu[,match(train_select_otus, names(train_otu))] #subset train_otu to only include otus which match the top 100 otus in test data


mo<-HOF(train_otu_subset, train_pH,family=poisson,test="AIC",model=c("I","II","III","IV","V"))
mop<-pick.model(mo)

#4. Predict abundances of top 100 OTUs in test dataset based on test pH

# build table with parameters
mod.out<-data.frame(mop,abunds.train) # a dataframe contain picked model ID  and train OTU ID

# predict abundances of each OTU in the test dataset (only top 100 OTUs most abundant are predicted)

mylist <- list()
for (i in row.names(mod.out)){
 mylist[[i]]<-predict(mo[[i]],model=paste0(mod.out[i,]$mop),newdata=test_pH)
}

#unlist the list to dataframe
pr.df <- do.call("rbind",mylist)
row.names(pr.df)<-row.names(mod.out)#make sure row names are the otu names
colnames(pr.df)<-row.names(test_pH)#make sure col names are the sample names
pr.df<-t(pr.df)#transpose

#5. 
mod.pr<-metaMDS(pr.df)#an NMDS on the top 100 OTUs in test dataset, ABUNDANCES are PREDICTED from pH and the train HOF models
mod.obs<-metaMDS(test_otu)# an NMDS on the full observed train_otu dataset.

plot(mod.pr$points[,1],mod.obs$points[,1]) 


