# POB_snMultiome_pipeline    
This snakemake pipeline was developed for processing snMultiome data and currently under active development. It has been developed solely on Biowulf.   

## 1.Overview    
The pipeline currently begins with tared cellranger output, completing per sample quality control, and multiple sample integration.  
![image](https://github.com/NCI-CCR-POB/POB_snMultiome_pipeline/assets/114779622/a60b67de-341c-47b0-83e5-7ac3025d3fad)  
**Overview of snMultiome analysis Pipeline**

## 2.Software/Package requirement  
The pipeline was built based on R4.3.2 and Python 3.9 
Please check if the R package in file `` or python packages in file `` have been installed. Please modify the path in `program.py` file.  
Please also check if the 


## 3.Usage  
  
  
## 4. References  
This repo was originally generated from the CCBR Nextflow Template. The template takes inspiration from nektool1 and the nf-core template. If you plan to contribute your pipeline to nf-core, don't use this template -- instead follow nf-core's instructions2.

## 5. Feedback
For comments/suggestions/advice please reach out to Ying Wu
