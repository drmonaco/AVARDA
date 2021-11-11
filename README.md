# AVARDA

#AVARDA Script
### Example AVARDA run
Rscript Desktop/CDI_db/Software/AVARDA/AVARDA.R Desktop/CDI_db/ProcessedData/test1/VirscanLar_000/test1_VirscanLar_000_Hits.tsv 1 Desktop/CDI_db/Software/AVARDA/bin2/my_df.csv ~/Desktop/CDI_db/Software/AVARDA/bin2/total_probability_xr2.csv Desktop/CDI_db/Software/AVARDA/bin2/unique_probabilities3.csv Desktop/CDI_db/Software/AVARDA/bin2/VirScan_filtered_virus_blast_new.csv Desktop/CDI_db/ProcessedData/test1/VirscanLar_000/AVARDA/ test_1_ Desktop/CDI_db/Software/AVARDA/bin2/avarda_names.csv

### Template
Rscript Desktop/CDI_db/Software/AVARDA/AVARDA.R Input1 Input2 Input3 Input4 Input5 Input6 Input7 Input8 Input9

### Inputs 
Input1: path to file of data to be analyzed, rows are peptides, columns are samples, values are generally binary values (indicating hits)
Input2: thesholdiag value for input 1; set to 1 if using binary matrix
Input3: path to "dictionary" indicating all inter-peptide alignments (doesn't change between runs)
Input4: path to precomp table listing each viruses abundance in VirScan (does not change if using Virscan)
Input5: path to precomp table listing each viruses relative representation compared with other viruses (doesn't change) 
Input6: path to precomp matrix of Virscan peptide - virus alignments
Input7: Output path
Input8: Name of study 
Input9: Table converting between dataset names and avarda names

### Instructions
To run avarda place directory in preferred location and set input arguments as absolute paths to various above listed files.

If using names other than Larman Lab names replace the first column of avarda_names.csv with other peptide names.
