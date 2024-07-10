#!/usr/bin/env /usr/local/Anaconda/envs/py3.9/bin/python
import psutil, argparse, os, subprocess, re, sys, itertools, glob
import snakemake
import os
import pandas as pd
from argparse import ArgumentParser

#print("Snakemake version:", snakemake.__version__)
# get the path of the script: e.g. /is2/projects/CCR-SF/active/Software/scripts/bin/
bin_path = os.path.dirname(os.path.realpath(__file__))
# get the main script path: e.g. /is2/projects/CCR-SF/active/Software/scripts/
#script_path = os.path.dirname(bin_path)
script_path = bin_path
#print(script_path)

from argparse import ArgumentParser

#refpath = f"{bin_path}/currentsnake/"
snakepath = f"{bin_path}"
print(f"pipeline source path is: {snakepath}")

class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

def main(raw_args=None):
    parser = argparse.ArgumentParser(description="""Help to set up and run the Single Cell Pipelines""", formatter_class=SmartFormatter)
    parser.add_argument("cellranger_path", metavar="rawdata",
        nargs='?', action = "store", type=str,
        help="R|Path to cellranger output directory \n")
    parser.add_argument("pipeline", metavar="qc",
        nargs='?', action = "store", type=str.lower,
        help="""R|Pipeline to use:
    qc for QC analysis
    scrna for Single Cell RNASeq
    snmultiome for Single nuclei Multiome
    scrnapdx for Single Cell RNASeq with PDX samples """)
    parser.add_argument("refgenome", metavar="hg38",
        nargs='?', action = "store", type=str,
        help="""R|Reference genome for the run
    hg38 for human
    mm10 for mouse
    PDX for PDX sample""")
    parser.add_argument("-p", "--projectname", metavar="projectname",
        action = "store", type=str, help="Name of project. [Required]")
    parser.add_argument("-a", "--aggregate", metavar="N",
        action = "store", type=str, help="If aggregate needs to be done.Type N if no, Type Y if yes [Required]")
    parser.add_argument("-atacmin", "--ATACminimalfiltering", metavar="N",
        action = "store", type=str, help="If ATAC filtering thrshold should set to minimal.Type N if no, Type Y if yes")
    parser.add_argument("-mem", "--memory", metavar="regular",
        action = "store", type=str, help="Setup cluster.yaml file; If 1-6 samples, Type regular; 7-15 samples, Type medium; >15 samples, Type large; Or setup based on cell numbers per samples  [Required]")
    parser.add_argument("-umi_min", "--nCount_RNA_LL", metavar="500",
        action = "store", type=int, help="Default is 500, modify as needed")
    parser.add_argument("-umi_max", "--nCount_RNA_UL", metavar="3MADS",
        action = "store", type=str, help="Default is 3MADS, modify as needed, the value can be numberMADS, max+number, number(hard cutoff)")
    parser.add_argument("-gene_min", "--nFeature_RNA_LL", metavar="300",
        action = "store", type=int, help="Default is 300, modify as needed")
    parser.add_argument("-mt_max", "--mt_UL", metavar="10",
        action = "store", type=int, help="Default is 10, modify as needed")
    parser.add_argument("-afrag_min", "--nCount_ATAC_LL", metavar="500",
        action = "store", type=int, help="Default is 500, modify as needed")
    parser.add_argument("-apeak_min", "--nFeature_ATAC_LL", metavar="300",
        action = "store", type=int, help="Default is 300, modify as needed")
    parser.add_argument("-afrag_max", "--nCount_ATAC_UL", metavar="3MADS",
        action = "store", type=str, help="Default is 3MADS, modify as neededmodify as needed, the value can be numberMADS, max+number, number(hard cutoff)")
    parser.add_argument("-nucleosome", "--nucleosome_signal_UL", metavar="2",
        action = "store", type=int, help="Default is 2, modify as needed")
    parser.add_argument("-tss", "--TSSenrichment_LL", metavar="2",
        action = "store", type=int, help="Default is 2, modify as needed")
    parser.add_argument("-blacklist", "--blacklist_fraction_UL", metavar="0.05",
        action = "store", type=float, help="Default is 0.05, modify as needed")
    parser.add_argument("-afrip", "--pct_reads_in_peaks_LL", metavar="0.15",
        action = "store", type=float, help="Default is 0.15, modify as needed")
    parser.add_argument("-npcs", "--npcs_val", metavar="50",
        action = "store", type=int, help="Default is 5, modify as needed")
    parser.add_argument("-gkanchor", "--gexkanchor", metavar="5",
        action = "store", type=int, help="Default is 5, modify as needed")
    parser.add_argument("-gkfilter", "--gexkfilter", metavar="200",
        action = "store", type=int, help="Default is 200, modify as needed")
    parser.add_argument("-gkweight", "--gexkweight", metavar="100",
        action = "store", type=int, help="Default is 100, modify as needed")
    parser.add_argument("-akanchor", "--atackanchor", metavar="5",
        action = "store", type=int, help="Default is 5, modify as needed")
    parser.add_argument("-akfilter", "--atackfilter", metavar="200",
        action = "store", type=int, help="Default is 200, modify as needed")
    parser.add_argument("-akweight", "--atackweight", metavar="100",
        action = "store", type=int, help="Default is 100, modify as needed")
    parser.add_argument("-umapres", "--umapresolution", metavar="0.8",
        action = "store", type=float, help="Default is 0.8, modify as needed")
    parser.add_argument("-parallel", "--parallel", metavar="N",
        action = "store", type=str, help="Default is N. If Run the steps before merge parallelly, Type Y. ")
    parser.add_argument("-n", "--dryrun", action="store_true", help="dry run")
    parser.add_argument("--rerun", action="store_true",
        help="dry run then prompt for submitting jobs that need to be rerun")
    parser.add_argument("--unlock", action="store_true",
        help="unlock working directory")
    args = parser.parse_args(raw_args)
    atac_minimal = args.ATACminimalfiltering if args.ATACminimalfiltering else "N"
    nCount_RNA_LL = args.nCount_RNA_LL if args.nCount_RNA_LL else 500
    nCount_RNA_UL = args.nCount_RNA_UL if args.nCount_RNA_UL else "3MADS"
    nCount_ATAC_LL = args.nCount_ATAC_LL if args.nCount_ATAC_LL else 500
    nCount_ATAC_UL = args.nCount_ATAC_UL if args.nCount_ATAC_UL else "3MADS"
    nFeature_RNA_LL = args.nFeature_RNA_LL if args.nFeature_RNA_LL else 300
    nFeature_ATAC_LL = args.nFeature_ATAC_LL if args.nFeature_ATAC_LL else 300
    nucleosome_signal_UL = args.nucleosome_signal_UL if args.nucleosome_signal_UL else 2
    TSSenrichment_LL = args.TSSenrichment_LL if args.TSSenrichment_LL else 2
    pct_reads_in_peaks_LL = args.pct_reads_in_peaks_LL if args.pct_reads_in_peaks_LL else 0.15
    blacklist_fraction_UL = args.blacklist_fraction_UL if args.blacklist_fraction_UL else 0.05
    mt_UL = args.mt_UL  if args.mt_UL else 10
    npcs = args.npcs_val if args.npcs_val else 50
    GEXkanchor = args.gexkanchor if args.gexkanchor else 5
    GEXkfilter = args.gexkfilter if args.gexkfilter else 200
    GEXkweight = args.gexkweight if args.gexkweight else 100
    ATACkanchor = 5
    ATACkfilter = 200
    ATACkweight = 100
    UMAPresolution = 0.8
    if args.atackanchor:
        ATACkanchor = args.atackanchor
    if args.atackfilter:
        ATACkfilter = args.atackfilter
    if args.atackweight:
        ATACkweight = args.atackweight
    if args.umapresolution:
        UMAPresolution = args.umapresolution

    if args.dryrun:
        snakemake.snakemake('Snakefile', unlock=True)
        snakemake.snakemake('Snakefile', printshellcmds=True, dryrun=True)  
    elif args.rerun:
        snakemake.snakemake('Snakefile', unlock=True)
        #snakemake.snakemake('Snakefile', printshellcmds=True, dryrun=True, rerun_triggers='mtime')
        snakemake.snakemake('Snakefile', printshellcmds=True, dryrun=True, rerun_triggers='unchanged')
        print('CONFIG FILE:')
        with open('config/config.py') as f:
            print(f.read())
        submit = input('Submit now? (y/n): ')
        if submit.lower() == 'y':
            subprocess.check_output('sbatch submit.sh', shell=True)
            print("Submitted snakemake rerun in directory")
    elif args.unlock:
        snakemake.snakemake('Snakefile', unlock=True)
    elif args.cellranger_path == None or args.pipeline == None or args.refgenome == None:
        parser.print_help()
    elif args.pipeline not in ['qc', 'scrna', 'snmultiome', 'scrnapdx']:
        parser.print_help()
        print(args.pipeline)
        print("\nPlease enter a value that corresponds to a pipeline")
    else:
        print("pipeline",args.pipeline)
        project_name = args.projectname
        rawdata_dir = args.cellranger_path
        print("cellranger rawdata path :",rawdata_dir)
        analysisPath =  os.getcwd()
        oneup = os.path.dirname(analysisPath)
        print("analysis path:",analysisPath)
        print("oneup path:", oneup)
        out = "library.csv"
        projectPath = os.path.abspath(analysisPath)
        print("project path:",projectPath)
        manifest_file = f'{analysisPath}/assets/input_manifest_cellranger.csv'
        df = pd.read_csv(manifest_file)
        print(df)
        df['uniqueID']=df['uniqueID'].str.replace(' ', '')
        samples_in_metadata=set(df['uniqueID'].tolist())
        print(samples_in_metadata)
        log_dir = f'{analysisPath}/logs'
        stat_dir = f'{analysisPath}/stats'
        config_dir = f'{analysisPath}/config'
        #create logs and stats folders
        os.makedirs(log_dir) if not os.path.exists(log_dir) else print(f"Directory '{log_dir}' already exists.")
        os.makedirs(stat_dir) if not os.path.exists(stat_dir) else print(f"Directory '{stat_dir}' already exists.")
        os.makedirs(config_dir) if not os.path.exists(config_dir) else print(f"Directory '{config_dir}' already exists.")
        # Create softlink
        directory_name = f'{analysisPath}/00_FullCellrangerOutputs'
        if not os.path.exists(directory_name):
            os.makedirs(directory_name)
            #print(f"Directory '{directory_name}' created successfully.")
        else:
            print(f"Directory '{directory_name}' already exists.")
            print(f"!!!!!!!!!!!!!!!!!!!!!!!!!Please check the softlink inside {directory_name}")

        OUT = open(out, "w")
        OUT.write("Sample_ID,Sample,subfolder,path\n")
        file2 = open(f'{analysisPath}/Aggregate.csv',"w")
        file2.write("library_id,atac_fragments,per_barcode_metrics,gex_molecule_info\n")
        # Create softlink
        for i in samples_in_metadata:
        #    print(i)
            sub = df.loc[df['uniqueID'] == i, 'subfolder'].values[0]
        #    print (sub)
            rawdata_path = f'{rawdata_dir}'
            if sub != 0 :
                rawdata_path = f'{rawdata_dir}/{sub}'
            
            source_path = f'{rawdata_path}/{i}.tar'
            destination_path = f'{analysisPath}/00_FullCellrangerOutputs/{i}.tar'
            # Check if the link already exists
            if not os.path.exists(destination_path):

                # Create the symbolic link
                os.symlink(source_path, destination_path)
                #print(f"Symbolic link created from '{source_path}' to '{destination_path}'.")
            else:
                print(f"Symbolic link already exists at '{destination_path}'.")
            OUT.write(i + "," + i + "," +  str(sub) + "," +source_path + "\n")
            file2.write(i + "," + f'{analysisPath}/cellranger_output/{i}/outs/atac_fragments.tsv.gz,' + f'{analysisPath}/cellranger_output/{i}/outs/per_barcode_metrics.csv,'+f'{analysisPath}/cellranger_output/{i}/outs/gex_molecule_info.h5' + "\n")
        OUT.close()
        file2.close()        

        
        #additionalInfo = False
        #Copying Snakefile to be used and printing messages to user if needed        
        if args.pipeline == 'snmultiome':
                # create config.py file
            with open('config/config.py', 'w') as f:
                f.write('rawdata="%s"\n' % f'{analysisPath}/cellranger_output')
                f.write('analysis="%s"\n' % analysisPath)
                f.write('sampleInfo="%s"\n' % f'{analysisPath}/assets/input_manifest_cellranger.csv')
                f.write('species="%s"\n' % str(args.refgenome))
                f.write('project="%s"\n' % args.projectname)
                f.write('npcs="%s"\n' % npcs)
                f.write('pipeline="%s"\n' % args.pipeline)
                f.write('Aggregate="%s"\n' % args.aggregate)
                f.write('ATAC_minimal="%s"\n' % atac_minimal)
                f.write('nCount_RNA_LL="%s"\n' % nCount_RNA_LL)
                f.write('nCount_RNA_UL="%s"\n' % nCount_RNA_UL)
                f.write('nFeature_RNA_LL="%s"\n' % nFeature_RNA_LL)
                f.write('mt_UL="%s"\n' % mt_UL)
                f.write('nucleosome_signal_UL="%s"\n' % nucleosome_signal_UL)
                f.write('TSSenrichment_LL="%s"\n' % TSSenrichment_LL)
                if atac_minimal == "N":
                    f.write('nCount_ATAC_LL="%s"\n' % nCount_ATAC_LL)
                    f.write('nCount_ATAC_UL="%s"\n' % nCount_ATAC_UL)
                    f.write('nFeature_ATAC_LL="%s"\n' % nFeature_ATAC_LL)
                    f.write('pct_reads_in_peaks_LL="%s"\n' % pct_reads_in_peaks_LL)
                else:
                    f.write('nCount_ATAC_LL="%s"\n' % 0)
                    f.write('nCount_ATAC_UL="%s"\n' % "max+1000")
                    f.write('nFeature_ATAC_LL="%s"\n' % 0)
                    f.write('pct_reads_in_peaks_LL="%s"\n' % 0)
                f.write('blacklist_fraction_UL="%s"\n' % blacklist_fraction_UL)
                f.write('memory="%s"\n' % args.memory)
                f.write('gkanchor="%s"\n' % GEXkanchor)
                f.write('gkfilter="%s"\n' % GEXkfilter)
                f.write('gkweight="%s"\n' % GEXkweight)
                f.write('akanchor="%s"\n' % ATACkanchor)
                f.write('akfilter="%s"\n' % ATACkfilter)
                f.write('akweight="%s"\n' % ATACkweight)
                f.write('resolution="%s"\n' % UMAPresolution)

                #subprocess.check_output('cp %s/submit.sh submit.sh' % snakepath, shell=True)
            scriptpath = f'{snakepath}/snMultiome_v2'
            if args.parallel == "Y" :
                scriptpath = f'{snakepath}/snMultiome_split_v2'
            cluster_yaml =  f'cluster_{args.memory}.yaml'
            print (scriptpath,"--------------")
            print (cluster_yaml)
            subprocess.check_output('cp %s/submit.sh submit.sh' % scriptpath, shell=True)
            subprocess.check_output('cp %s/run_snakemake_sc_v2.py run_snakemake_sc.py' % snakepath, shell=True)
            subprocess.check_output('cp -r %s/scripts ./' % scriptpath, shell=True)
            #subprocess.check_output('cp %s/config/{cluster_yaml} config/cluster.yaml' % scriptpath, shell=True)
            subprocess.check_output(f'cp -r {scriptpath}/config/{cluster_yaml} config/cluster.yaml', shell=True)
            subprocess.check_output(f'cp -r {scriptpath}/config/program.py config/program.py', shell=True)
            if args.aggregate == "Y":
                subprocess.check_output('cp %s/Snakefile Snakefile' % scriptpath, shell=True)
            else:
                subprocess.check_output('cp %s/Snakefile_noAggr Snakefile' % scriptpath, shell=True)
             
        additionalInfo = True
        if additionalInfo:
            print('Exiting now. Fill in additional information in config file, then rerun with ./run_snakemake_sc.py --rerun')
        else:
            stdout = sys.stdout
            sys.stdout = open('dryrun_all_commands.txt', 'w')
            snakemake.snakemake('Snakefile', printshellcmds=True, dryrun=True)
            sys.stdout = stdout
            snakemake.snakemake('Snakefile', printshellcmds=True, dryrun=True)
            print('CONFIG FILE:')
            with open('config/config.py') as f:
                print(f.read())

            submit = input('Submit now? (y/n): ')
            while submit.lower() not in ['y', 'n']:
                submit = input('Please answer yes(y) or no(n): ')
            if submit.lower() == 'y':
                subprocess.check_output('sbatch submit.sh', shell=True)
                print("Submitted snakemake run for %s in %s" % (args.fastq_path, work_path))

if __name__ == '__main__':
    main()

