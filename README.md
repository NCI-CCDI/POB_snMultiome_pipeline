# POB_snMultiome_pipeline    
This snakemake pipeline was developed for processing snMultiome data and currently under active development. It has been developed solely on Biowulf.   

## 1.Overview    
The pipeline currently begins with tared cellranger output, completing per sample quality control, and multiple sample integration.  
![image](https://github.com/NCI-CCR-POB/POB_snMultiome_pipeline/assets/114779622/a60b67de-341c-47b0-83e5-7ac3025d3fad)  
**Overview of snMultiome analysis Pipeline**

## 2.Setup Dependencies
- 1) The R and python used in the pipeline are **R4.3.2** and **Python 3.9**     
- 2) Please check if the R package in file `workflow/snMultiome_split_v2/scripts/R/main/load_packages.R` or python packages in file `workflow/snMultiome_split_v2/scripts/python/check_py.py` have been installed.    
- 3) Please modify the path in `workflow/snMultiome_split_v2/config/program.py` file.     
- 4) Please check if [**reticulate**](https://rstudio.github.io/reticulate/) has been installed and modify the path in `workflow/snMultiome_split_v2/scripts/R/main/1.Multiome_preprocess_rawGEX_wSoupX_perSample.R` file. reticulate is the R interface to use Python script. Python based software [**scrublet**](https://github.com/swolock/scrublet) was used in R script for double removal.      

 
## 3.Usage
* 1) Preparing Input Manifest File in `assets/input_manifest_cellranger.csv`   
     “uniqueID”, “prefix”, "masterID",“subfolder” are required in the Column 1, 2,3 and last column. All the other fields between masterID and subfolder will be added into seurat object metadata.     
     example:   
     <img width="1328" alt="image" src="https://github.com/NCI-CCR-POB/POB_snMultiome_pipeline/assets/114779622/2fa836bc-4e58-418d-bfc8-d5c8abc8f1b5">      
     **subfolder**: is the folder that contains the tared cellranger output. The full path which will be used in the wrapper for looking for the tar files are project_rawdata_dir/subfolder. Please modify the wrapper in order to generate softlink of the tar files.   
     

           
* 2) Use the wrapper to generate softlink of the tar files, generate config files, copy the scripts to the analysis folder   
    **Please modify the source path in the wrapper!!!!!**   
    ```
    export PATH=$PATH:/path/to/snMultiomePipeline #this is the path which contain the run_snakemake_sc_v2.py  
    run_snakemake_sc_v2.py -p ccrccdi4 -a N -atacmin Y -mem regular -umapres 0.2 -parallel Y /data/CCRCCDI/rawdata/ccrccdi4 snmultiome hg38 #This step only prepare the analysis folder  
    run_snakemake_sc_v2.py --rerun #This step will start the pipeline   
    ```
    
    run_snakemake_sc_v2.py --help
    <img width="1455" alt="image" src="https://github.com/NCI-CCR-POB/POB_snMultiome_pipeline/assets/114779622/b3ce137b-5bf6-4deb-9ef2-e7c864c21064">   

        
    

* 3) Here is the detailed [**documentation**](https://github.com/NCI-CCR-POB/POB_snMultiome_pipeline/blob/main/snMultiome%20Pipeline%20Instruction_v2.docx)    

         
  
## 4. References    
* This Pipeline was built by [**Ying Wu**](ying.wu@nih.gov) based on vigenette from [**Seurat**](https://satijalab.org/seurat/) and [**Signac**](https://stuartlab.org/signac/)    
* Incorporate the snMultiome pipeline from POB staff scientist [**Xiyuan Zhang**](xiyuan.zhang@nih.gov)      
* Incorporate several functions from CCBR [**SINCLAIR**](https://github.com/CCBR/SINCLAIR) pipeline   

## 5. Feedback
For comments/suggestions/advice please reach out to [**Ying Wu**](ying.wu@nih.gov)  
