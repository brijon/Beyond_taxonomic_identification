#Example blastn command used to query independent datasets againts reference dataset to gain 'Query OTU percentage of hits' seen in Table.1
blastn -query <query_fasta> -db <database_of_reference_repseqs> -num_alignments 20 -evalue 0.001 -outfmt 7  
