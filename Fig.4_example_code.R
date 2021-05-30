# Example script for Fig.4

# Script produces a phylogenetic tree from OTU repseqs and visualised tree as a circle plot with pH classification and phyla annotations (using graphlan)

#script assumes you have the following variables in working environment
#1) model_features table containing two collumns OTU name and pH class based on HOF model optima (to run HOF models/ assign pH class see HOF_example_code.R)
#2) taxonomy table containing two collumns OTU name and Phyla
#3) sequencepath variable containing full path to fasta file (including fasta file name) of all OTU repseqs 
#4) variable called outdir containing path to directory to write alignment, tree and  graphlan inputs and outputs to
#note ensure that OTUS in taxonomy table and model_features table are in fasta file , otherwise when you run graphlan will likely result in errors
#===========================================================================================================================================================================================
#load in tidyr library
library(tidyr)
#===========================================================================================================================================================================================

#part a) linux commands to run sequence alignment and make tree 
#can ofcourse run these commands directly from command line but have included in R script for simiplicity
#these commands are assuming clustalo and fastree are in linux $PATH variable if they are not either add them to $PATH or edit command e.g change clustalo to /path/to/executable/clustalo

#align sequences
#define alnpath variable to save alignment to
alnpath=paste0(outdir,"/otus.aln")
#run alignment
system(paste0("clustalo -i ",sequencepath," -t DNA -o ",alnpath))
#or directly in command line
#clustalo -i <sequencepath> -t DNA -o <alnpath>

#generate phylogenetic tree
#define treepath variable to save tree to
treepath=paste0(outdir,"/otus.tre")
#run fasttree
system(paste0("fasttree -nt ",alnpath," > ",treepath))
#or directly in command line
#fasttree -nt <alnpath> > <treepath>

#===========================================================================================================================================================================================

#part b) Make mapping files
#should now have tree file we can use with graphlan
#now need to make mapping files to annotate tree with pH class and phyla

#first we will create a mapping file in graphlan format to annotate each OTU's pH class
#it should fit the following format:
#e.g OTU | "ring colour" | "1" | colour

#initially get dataframe with otus as first collumn then cols repeating "ring colour" and 1's (referring to ring 1) and then the pH class
ph_ring<-as.data.frame(cbind(as.character(model_features[,1]),rep('ring_color',nrow(model_features)),rep(1,nrow(model_features)),as.character(model_features[,2])),stringsAsFactors =FALSE)

#graphlan expects last col to be a colour not a character string
#lets assign each pH class to a colour
ph_ring$V4[ph_ring$V4=="Acid"]<-"#FF0000"
ph_ring$V4[ph_ring$V4=="Acid to Mid"]<-"#E67E22"
ph_ring$V4[ph_ring$V4=="Acid to Neutral"]<-"#F4D03F"
ph_ring$V4[ph_ring$V4=="Mid"]<-"#7CFC00"
ph_ring$V4[ph_ring$V4=="Mid to Neutral"]<-"#32CD32"
ph_ring$V4[ph_ring$V4=="Neutral"]<-"#20B2AA"
ph_ring$V4[ph_ring$V4=="None"]<-"#BFC9CA"

#write mapping file 
write.table(ph_ring,file=paste0(outdir,"/ph_anno.txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

#generate taxonomic anno mapping file 

#make colour pallete vector, these colours will be used to anotate otus belonging to dominant phyla (they will be used to assign node colour)
colourset=c("#CD6155","#5499C7","#F1C40F","#7DCEA0","#EB984E","#C39BD3","#9B59B6","#F5B041","#5993E5","#A6D785","#CD5555","#F87531","#FFF68F","#D3D3D3","#8B5742","#F6C9CC","#f28585","#f7be9e","#b79284","#e3ffd1","#bc6340","#d12121","#ffd39b","#a2cd5a","#cd950c","#ff7ff0","#a52a2a","#5f9ea0","#8a2ee2","#98f5ff","#8b7355","#f0f8ff","#7fffd4","#66cdaa","#8b8378","#8b7d6b","#f5f5dc","#ffebcd","#c1cdcd","#838b8b","#CD6155","#5499C7","#F1C40F","#7DCEA0","#EB984E","#C39BD3","#9B59B6","#F5B041","#5993E5","#A6D785","#CD5555","#F87531","#FFF68F","#D3D3D3","#8B5742","#F6C9CC","#f28585","#f7be9e","#b79284","#e3ffd1","#bc6340","#d12121","#ffd39b","#a2cd5a","#cd950c","#ff7ff0","#a52a2a","#5f9ea0","#8a2ee2","#98f5ff","#8b7355","#f0f8ff","#7fffd4","#66cdaa","#8b8378","#8b7d6b","#f5f5dc","#ffebcd","#c1cdcd","#838b8b","#CD6155","#5499C7","#F1C40F","#7DCEA0","#EB984E","#C39BD3","#9B59B6","#F5B041","#5993E5","#A6D785","#CD5555","#F87531","#FFF68F","#D3D3D3","#8B5742","#F6C9CC","#f28585","#f7be9e","#b79284","#e3ffd1","#bc6340","#d12121","#ffd39b","#a2cd5a","#cd950c","#ff7ff0","#a52a2a","#5f9ea0","#8a2ee2","#98f5ff","#8b7355","#f0f8ff","#7fffd4","#66cdaa","#8b8378","#8b7d6b","#f5f5dc","#ffebcd","#c1cdcd","#838b8b","#CD6155","#5499C7","#F1C40F","#7DCEA0","#EB984E","#C39BD3","#9B59B6","#F5B041","#5993E5","#A6D785","#CD5555","#F87531","#FFF68F","#D3D3D3","#8B5742","#F6C9CC","#f28585","#f7be9e","#b79284","#e3ffd1","#bc6340","#d12121","#ffd39b","#a2cd5a","#cd950c","#ff7ff0","#a52a2a","#5f9ea0","#8a2ee2","#98f5ff","#8b7355","#f0f8ff","#7fffd4","#66cdaa","#8b8378","#8b7d6b","#f5f5dc","#ffebcd","#c1cdcd","#838b8b")

#get most dominant phyla to annotate
Phyla_of_interest<-names(sort(table(taxonomy$Phylum),decreasing = T)[1:25])


#loop over phyla of interest
for (i in 1:length(Phyla_of_interest)){
  Phyla<-Phyla_of_interest[i]
  #assign colour to all otus belonging to that phyla
  taxonomy[which(taxonomy[,2]==Phyla),3]<-colourset[i]
  
}

#remove rows with no node colour assigned as they are OTU's belonging to rarer phyla
#these otus will be visualised in tree they just wont be assigned a node colour
taxonomy<-taxonomy[!is.na(taxonomy$V3),]

#now lets convert this to graphlan mapping file format
#it should fit the following format:
#e.g OTU | phyla_name | colour

taxonomy_anno<-cbind(taxonomy[,1],rep('clade_marker_color',nrow(taxonomy)),taxonomy[,3])
#write mapping file 
write.table(taxonomy_anno,file=paste0(outdir,"/taxonomy_anno.txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)


#get phyla colours for key at side of circle plot
Phyla_key<-as.data.frame(cbind(Phyla_of_interest,rep('clade_marker_color',length(Phyla_of_interest)),colourset[1:length(Phyla_of_interest)]),stringsAsFactors =FALSE)
#write mapping file 
write.table(Phyla_key,file=paste0(outdir,"/phyla_key.txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

#note all mapping files could be combined for example by using rbind before writing out
#you may also want to add additional aesthetic parameters in regards to node size, text size etc
#read graphlan's documentation for specific parameters https://github.com/biobakery/graphlan

#===========================================================================================================================================================================================

#c) run graphlan

#these commands are assuming graphlan is in linux $PATH variable if not add to linux $PATH variable or edit command e.g from graphlan_annotate.py to /path/to/executable/graphlan_annotate.py
system(paste0("graphlan_annotate.py --annot ",outdir,"/ph_anno.txt ",treepath))
system(paste0("graphlan_annotate.py --annot ",outdir,"/taxonomy_anno.txt ",treepath))
system(paste0("graphlan_annotate.py --annot ",outdir,"/phyla_key.txt ",treepath))
#or directly in command line
#graphlan_annotate.py --annot <annofile> <treepath>

system(paste0("graphlan.py --dpi 350 --size 20 ",treepath," ",outdir,"/circleplot.pdf"))
#or directly in command line
#graphlan.py --dpi 350 --size 20 <treepath> <pdfout>
