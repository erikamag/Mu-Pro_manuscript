### Nextflow RNA-Seq pipeline

[NF-Core RNASeq](https://github.com/nf-core/rnaseq)

- nfc/rnaseq is a bioinformatics analysis pipeline used for RNA sequencing data.
- The workflow processes raw data from FastQ inputs (FastQC, Trim Galore!), aligns the reads (STAR or HiSAT2), generates counts relative to genes (featureCounts, StringTie) or transcripts (Salmon, tximport) and performs extensive quality-control on the results (RSeQC, Qualimap, dupRadar, Preseq, edgeR, MultiQC). See the output documentation for more details of the results.
- The pipeline is built using Nextflow, a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

[nf-core/rnaseq usage docs](https://nf-co.re/rnaseq/3.4/usage)

1. **Copy RNAseq raw read .fq.gz files to one directory, where pipeline will be run**
2. **Create meta.tsv file for raw reads** 
		
		SampleID    Tissue    Genotype    Treatment    Replicate    source    paired    spots    avgLength    r0    r1    r2
		* spots: PF clusters from Illumina data, clusters passing filter or reads passing filter 
		* Treatment = (e.g. mutant vs. control)
			source: local = file path for local fq.gz files, sra = SRA submission number
		* paired = PE or SE
		* avgLength = SRA submission stats
		* r0: if sra, then r0 = SRA#, if local + SE, r0 = file path to SE fq.gz data, if local + PE, leave r0 blank
		* r1 = if local + PE - left or R1 reads
		* r2 = if local + PE - right or R2 reads

3. **Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) - MacOs Intall [Minconda](https://docs.conda.io/projects/conda/en/4.6.1/user-guide/install/macos.html)**

	    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
	    bash Miniconda3-latest-Linux-x86_64.sh 
	    
	Can also added [conda channels](http://bioconda.github.io/user/install.html#set-up-channels)
	After installing conda you will need to add the bioconda channel as well as the other channels bioconda depends on. It is important to add them in this order so that the priority is set correctly (that is, conda-forge is highest priority).
	The conda-forge channel contains many general-purpose packages not already found in the defaults channel.

        conda config --add channels defaults
        conda config --add channels r
        conda config --add channels bioconda
        conda config --add channels conda-forge
    
        cat .condarc #check ~/.condarc file to make sure the conda channels are in order that they were installed. with last one installed at top
        # channels:
          - conda-forge
          - bioconda
          - r
          - default

	[Update Miniconda](https://docs.anaconda.com/anaconda/install/update-version/) if already installed Already had miniconda3 installed in home directory, just need to update Miniconda

	    conda update --all # This updates all packages in the current environment to the latest version.
	    cd ~/miniconda3/envs
	    conda update -n myenv --all # if you have other environments that need to be updated 

4. **Install [Nextflow](https://github.com/nextflow-io/nextflow) (in a new conda environment)**

	    conda create -n nf python=3 #name of new environment = n, create conda environment nf.
	    conda activate nf
	    conda install nextflow # if this doesn't work can use -c option to get in correct environment # conda install -c bioconda nextflow
	    conda list
5. **Set up nextflow pipeline directory**

	Pipeline scripts/files in `miniconda/envs/nf/rnaseq` include:

	`genomes.yml  nextflow.config  meta.tsv  # meta.tsv file includes meta data on RNA-seq reads, see example file in rnaseq_input_scripts`

6. **Add these environmental variables (with necessary modification) to your `~/.bashrc` or `~/.bash_profile`**

        cd ~
        ls -la # will show hidden files
        vim .bashrc
        i # Insert within Vim
        esc # quote vim insert
        :wq # save and quit
    
        export NXF_HOME=/home/springer/magnu513/git/nf
        export NXF_CACHE=/scratch.global/magnu513/nf
        export NXF_EXECUTOR=slurm
        export NXF_CONDA_CACHEDIR=/home/springer/magnu513/miniconda3/envs
        export NXF_WORK=$NXF_CACHE/work
        export NXF_TEMP=$NXF_CACHE/tmp
        export NXF_SINGULARITY_CACHEDIR=$NXF_CACHE/singularity
        export NXF_OPTS='-Xms1g -Xmx10g'
    
    Log out and log in again (or run `source ~/.bashrc`) to make these variables into effect

7. **Edit `~/git/nf/configs/rnaseq.config`**
    
	    process.conda = "$NXF_CONDA_CACHEDIR/rnaseq"
	    
	    params {
	      qcdir = "$NXF_HOME/qc"
	      s3dir = "$NXF_HOME/s3"} 

8. **Create a conda environment named rnaseq to run the pipeline**

	This can be done either using an existing environment yml file:
	
	        conda env create -n rnaseq -f $NXF_HOME/configs/environments/rnaseq.yml
	        conda env list
	        # you should now see a new environment named "rnaseq"

	Or create an empty environment and then manual install the packages:
	
	        conda create -n rnaseq
	        conda env list
	        # you should now see a new environment named "rnaseq"
	        
	        conda activate rnaseq
	        conda install sra-tools pigz fastqc trim-galore bwa hisat2 star picard samtools bcftools bedtools bamtools pysam sambamba preseq alfred bioawk biopython deeptools jcvi kallisto minimap2 numpy pandas perl plotly pyfaidx pybigwig ucsc-bedgraphtobigwig subread r-base r-tidyverse r-glue r-argparse r-purrr r-readr r-readxl rseqc salmon stringtie multiqc

    If run into package inconsistencies

        conda update --all

    then for specific packages with inconsistencies.

        conda install {package_name}

	In order to create this conda environment, because it takes a long time use an [interactive job](https://www.msi.umn.edu/content/interactive-queue-use-qsub) [interactive job](https://www.msi.umn.edu/news/mesabi-interactive-queue-replacing-labqi-cluster)
	
	    srun -N 1 -n 1 -c 1 -t 180 --mem=20gb -p interactive --pty bash     # this gives you 3 hours (-t 180) interactive session, -N: number of nodes, -n: number of tasks per node, -c: cpus-pre-task--mem: memory required per node # once session starts input code
	    cd ~git/ #whatever directory contains nf OR refer to above $NXF_HOME from .bashrc file
	    conda env create -n rnaseq -f nf/configs/environments/rnaseq.yml #create new environment - lists all software dependencies need. # this process might take some time.  # rnaseq.yml file contains all the packages needed for the nextflow pipeline, conda will use this file to install all packages. # config.yml files are YAML formatted to contain configuration options
	    conda env list # you should now see a new environment named "rnaseq" # should see nf and rnaseq
	    conda activate rnaseq
	    conda list # all packages needed should be in this rnaseq env

9. **Make necessary changes to pipeline input files:
`~/git/nf/`**
	- `meta.tsv` which contains paths to your fastq sequences
	- `nextflow.config` 
	- `genomes.yml` # genome index configuration, make sure genome.yml file contains all files for genome mapping to

	`nextflow.config` **EDIT**
	
		launchDir = "/home/springer/magnu513/git/nf/rn21a/rnaseq" **EDIT** #wherever the `meta.tsv` file is
		
		workDir = "$NXF_CACHE/work/test"
			
		includeConfig "$NXF_HOME/configs/nextflow.config"
		includeConfig "$NXF_HOME/configs/fastq.config"
		includeConfig "$NXF_HOME/configs/rnaseq.config"

		params {
		  genome = 'Zmays_W22' # check to make sure genomes.yml file has genome want
		  name = 'rn21a' # name of project
		  design = 'rn21a.tsv' # meta.tsv
		  // output locations
		  outdir = "./raw"
		  qcdir = "./qc"
		  s3dir = "./s3"
		  // sequence source: "local", "sra", "s3" or "mixed"
		  source = 'sra' 
		  // paired-end or single-end: "SE", "PE" or "mixed"
			  paired = 'PE' 
			  // strand-specific RNA-Seq? :"no", "forward" or "reverse"
			  stranded = 'no'
			  interleaved = false
			  save_fastq = false
			  save_trimmed = false
			  aligner = "hisat2" // one of "hisat2" or "star"
			  saveBAM = true  #if want to look at raw reads in IGV
			  skip_preseq = true
			  run_salmon = false
			  run_stringtie = false
			  ase = false
			  ril = false
			  cage = false
			  // send email?
			  email = false
			  email_on_fail = 
			}

10. **Finally, inside directory where `meta.tsv` is `~git/nf/rn21a/rnaseq` we are set to run the test pipeline using existing genome database** 
    
	Submit a interactive job or a job- `interactive.sh`
	[Interactive Job-Slurm](https://www.msi.umn.edu/sites/default/files/Slurm_Workshop.pdf)
		
		#SBATCH -N 1
		#SBATCH -n 1
		#SBATCH -c 1
		#SBATCH --mem=150g
		#SBATCH -p interactive
	
	    cd /home/springer/magnu513/git/nf/rn21a/rnaseq
	    
	    conda activate nf
	    
	    nextflow run $NXF_HOME/rnaseq -params-file genomes.yml -profile mangi
	
	Submit interactive job with
	
	    srun --pty bash
	    
 Or submit job with the following shell script
 
		#!/bin/bash -l
		#SBATCH --time=90:00:00
		#SBATCH --ntasks=8
		#SBATCH --mem=245g
		#SBATCH --tmp=200g
		#SBATCH --job-name=rn21a
		
		cd ~/git/nf/rn21a/rnaseq
		
		conda activate nf
		
		nextflow run $NXF_HOME/rnaseq -params-file genomes.yml -profile mangi

**Results** results should be in `nextflow.config` outdir, `~/raw`

`~/git/nf/rn21a/rnaseq`

    rn21a.tsv  genomes.yml  nextflow.config  qc  raw  reads.tsv  reads.xlsx  results  s3
    
    ../qc/Zmays_W22/rn21a/ #contains 00.meta.tsv  00.raw.rds # read into .R
    ../s3/Zmays_W22/ #contains rn21a.html #qc of data in html format
    ../raw/ # contains .bam files and all data processed
        00.meta.tsv  02_trim_galore  20_bam       26_qualimap  31_featureCounts       40_multiqc  pipeline_info 01_fastqc    11_hisat2       26_dupradar  26_rseqc     33_sample_correlation  50_final  # 11_hisat2 directory should be in the outdir if BAM files saved = true`

 **Other pipeline codes**

    conda deactivate # to deactivate current env
    
    -resume # resume to pipeline, add onto submission run code


**Errors when running pipeline**
    
    cd /scratch.global/magnu513/nf/work/test/da/7d659c3acf48bd20cac75eff28593d
    less .command.err
    
**Error 1:** Caused by: java.io.IOException: No space left on device.

This error is saying either you ran out of space on the scratch folder, the server node, or your home directory

Looking at the `command.run.txt` file, the script is requesting just under 16GB of RAM. Slurm defaults to a temporary directory allocation of half of the requested RAM, so it would be just under 8GB. If you can add in a directive like `#SBATCH --tmp=24gb` then that might address the issue. Also, I think increasing the amount of memory requested for the qualimap step wouldn’t hurt the job, either.

 `less command.run.txt`
		
		#!/bin/bash
		#SBATCH -D /scratch.global/magnu513/nf/work/test/da/7d659c3acf48bd20cac75eff28593d
		#SBATCH -J rnaseq:qmap.s06
		#SBATCH -o /scratch.global/magnu513/nf/work/test/da/7d659c3acf48bd20cac75eff28593d/.command.log
		#SBATCH --no-requeue
		#SBATCH -t 25:00:00
		#SBATCH --mem 15360M
		#SBATCH -p amdsmall
		#NEXTFLOW TASK: rnaseq:qmap (s06)
		set -e
		set -u
		NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
		NXF_ENTRY=${1:-nxf_main}
		
		Error executing process > 'rnaseq:qmap (C07)'
		
		Caused by:
		  Process `rnaseq:qmap (C07)` terminated for an unknown reason -- Likely it has been terminated by the external system
		
		Command executed:
		
		  unset DISPLAY
		  qualimap --java-mem-size=155G rnaseq -p non-strand-specific -pe \
		    -bam C07.bam -gtf 10.gtf -outdir C07
		
		Command exit status:
		
		Command output:
		  (empty)
		
		Work dir:
		  /scratch.global/magnu513/nf/work/test/77/4735bf6184a00bc89e254921c63585
		
Fix to the error - Add the following section to `~/nextflow.config` (no need to change anything under `~git/nf/`) to increase memory for the qualimap step. Add this section after the params {} section**

    # have increased GB to 200GB for check_max and --tmp
    process {
      withName:qmap {
        cpus = { check_max( 1, 'cpus' ) }
        memory = { check_max( 24.GB + 5.GB * task.attempt, 'memory' ) }
        time = { check_max( 8.h * task.attempt, 'time' ) }
        clusterOptions = '--tmp=24gb'
      }
    }

**Error 2:** When running qmap step and trying to `-resume` same job**

When you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

`~/git/nf/rn21a/rnaseq/rn21a.nextflow.sh`

    nextflow run $NXF_HOME/rnaseq -params-file genomes.yml -profile mangi -resume # added the -resume option to script and increased --tmp memory.

**Error 3:** following error when attempting to run `-resume` on same job because .nextflow/cache/[directory]/db is locked

[lock on nextflow session](https://github.com/nextflow-io/nextflow/issues/1300)
	
	`N E X T F L O W  ~  version 20.10.0
	Launching /home/springer/magnu513/git/nf/rnaseq/main.nf [nasty_hawking] - revision: 78f0a943da
	Unable to acquire lock on session with ID 726edee9-cbd7-44c2-b13d-61a3154a2f76`

Common reasons of this error are:
 - You are trying to resume the execution of an already running pipeline
 - A previous execution was abruptly interrupted leaving the session open

You can check what process is holding the lock file by using the following command:

	lsof /panfs/roc/groups/15/springer/magnu513/git/nf/rn21a/rnaseq/.nextflow/cache/726edee9-cbd7-44c2-b13d-61a3154a2f76/db/LOCK

 In order to view this directory with a LOCK ` ~/git/nf/rn21a/rnaseq`

    conda activate nf
    cd /.nextflow/cache/726edee9-cbd7-44c2-b13d-61a3154a2f76/db/
    ls # if see LOCK file, remove it
    rm -rf LOCK

In order to resubmit the job properly with `-resume` after removing LOCK file in `.nextflow/cache` directory
[Nextflow `-resume`](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html)

You can use the resume command with either the session ID or the run name to recover a specific execution. For example:

    nextflow run rnaseq-nf -resume mighty_boyd
or

    nextflow run rnaseq-nf -resume 4dc656d2-c410-44c8-bc32-7dd0ea87bebf

`~/git/nf/rn21a/rnaseq/rn21a.nextflow.sh` # found which [run name] was running and failed and added line the `-resume` command to script after adjusting qmap memory in `nexflow.config`

    cd ~/git/nf/rn21a/rnaseq/
    head slurm-1044918.out
    
    N E X T F L O W  ~  version 20.10.0
    Launching `/home/springer/magnu513/git/nf/rnaseq/main.nf` [wise_bassi] - revision: 78f0a943da
    Both GTF and GFF have been provided: Using GTF as priority.
    ----------------------------------------------------
      orionzhou  rnaseq v0.9.1
      Nextflow RNA-Seq analysis pipeline
      https://github.com/orionzhou/nf
    ----------------------------------------------------

added run name to `rn21a.nextflow.sh`

    nextflow run $NXF_HOME/rnaseq -params-file genomes.yml -profile mangi -resume wise_bassi
    
    sbatch rn21a.nextflow.sh # resubmitted script


**If there are package dependency issues when trying to update conda environment**

    srun -N 1 -n 1 -c 1 -t 180 --mem=20gb -p interactive --pty bash # interactive job    # this gives you 3 hours (-t 180) interactive session, -N: number of nodes, -n: number of tasks per node, -c: cpus-pre-task--mem: memory required per node # once session starts input code

    conda env update -n rnaseq -f $NXF_HOME/configs/environments/rnaseq.yml # ran this instead in an interactive job and it worked to solve the software dependencies
    
	ERROR: package dependencies. Ran into an error where there was a package that needed to be downgraded
	pkg_resources.ContextualVersionConflict: (decorator 5.1.0 (`/panfs/roc/groups/15/springer/magnu513/miniconda3/envs/rnaseq/lib/python3.7/site-packages`), Requirement.parse('decorator<5,>=4.3'), {'networkx'})
	went to working directory `/scratch.global/magnu513/nf/work/rnaseq/rn20a_b/23/06289d12bcd90991f901b332672c4e``

	less .command.log

needed to find the exact package the `rnaseq.yml` file was wanting for decorator. wanted decorator 4.4.2 versus decorator=>5. Conda would not downgrade. 

So found exact label for package needed on conda 
[decorator label](https://anaconda.org/conda-forge/decorator)

	conda install -c conda-forge/label/cf202003 decorator
	conda list #still had decorator 5.1.0 installed
	cd ~/miniconda3/envs/rnaseq/lib/python3.7/site-packages
	ls #see both decorator=5.1.0 and decorator=4.4.2
	remove -r decorator-5.1.0.dist-info 