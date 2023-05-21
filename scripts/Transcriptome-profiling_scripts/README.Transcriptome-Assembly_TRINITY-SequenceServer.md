### Transcriptome-Assembly with TRINITY

**Transcriptome Assembly Workflow: **

1. Assemble + clean raw reads from fastq file (.fastq.gz” or “.fq.gz”) with Trinity 
2. Create Trinity.blastdb from Trinity.fasta files with Sequence Server
3. Build GMAP index 
4. Align Trinity output to the W22 genome using GMAP 

Note: Did not find mapping assembled transcripts to the genome as informative as identifying transcripts with local BLAST database

**1. De novo assembly with Trinity**
[Documentation](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity)

Assemble + clean raw reads from fastq file (.fastq.gz” or “.fq.gz”) with Trinity 

1. Make sure all .fq.gz read files are in one directory
	
	R1 and R2 files for each sample because there are two reads from a single read pair with paired-end sequencing 

2. Prepare tab-delimited samples.txt files for Trinity to recognize multiple reps/TF mutant

 create a `samples.txt` file that describes the data. This way you do no need to supply the command line with each rep by PE read.   the `—samples_file argument` can call this file
    
	    --samples_file <string>       tab-delimited text file indicating biological replicate relationships.
	                            cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
	                            cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
	                            cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
	                            cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq
         
                            
	I only want to include 1 condition or one TF allele x 3 reps in my samples.txt files. Trinity will assemble all samples included in the samples.txt file into one de novo assembly. A potential reason why multiple conditions are used for one assembly includes: a genome assembly where multiple tissues were sampled to represent one genome. 
	
	[Trinity generates a single assembly from all of the reads](https://github.com/trinityrnaseq/trinityrnaseq/issues/320)
	
	For SE only include one .fq file/replicate

3. Create a Trinity .sh script**

`trinity.sh` example

    module load trinity
    module load trimmomatic/0.33

    cd /home/springer/magnu513/rn20h/rn20h_00_fastq/

    Trinity --seqType fq --max_memory 150G --samples_file 24.baf6021m1.samples.txt --CPU 6 --trimmomatic --min_contig_length 200 --output trinity.out.baf6021m1

Check to make sure there are no wierd^M characters before running shell script

    cat -te filename

There are several ways to translate a file between ASCII CR+LF (DOS/Windows) and LF (Unix) newlines. The following are different options to remove, convert or translate the ^M characters:
The simplest solution, use the dos2unix command (sometimes named fromdos, d2u or unix2dos):

    dos2unix filename
    
    sbatch *.sh #submit shell script
    

**2. Sequence Server. Make local BLAST database**

[Use sequence server in local terminal](https://sequenceserver.com/#installation)

Sequence Server Installation Requirements
Linux or Mac and Ruby (≥ 1.8.7; preferably ≥ 2.0)

NCBI BLAST+ (2.2.31+) is interactively downloaded if absent

One of SequenceServer's dependencies (the module used to parse BLAST's XML output) compiles some C code as part of the installation process. This means that the standard Unix build tools (e.g., gcc, make) are required to install SequenceServer. On a Mac, this means having Xcode and CLI tools for Xcode installed. On Ubuntu and other Debian-based Linux systems, you would have to install the ruby-dev and build-essential packages in addition to ruby.

[Install Xcode Command Line Tools via terminal](https://wilsonmar.github.io/xcode/#InstallIDE)

	xcode-select --install
	
Install and update sequence server 

    sudo gem install sequenceserver # this will ask you if you want to download NCBI+, click Enter and it will install 

Create a folder or directory on your local computer to store your TRINITY.fasta file(s) you want to create a BLAST database from

go into that directory via terminal

	cd ~/genome/ # this should be where your TRINITY.fasta file is
	
Configure and run sequence server 

    sequenceserver # this should create a blastdb in the directory you are in and a web page should pop up
  
If sequence server does not know where to create the database use the following option 

    sequenceserver -d [file path to .fasta file(s)] # if all blastdb are made, then sequence server will open properly without '-d' option  
        sequenceserver -d ~/genome/

    sequenceserver -m -d [file path to .fasta file(s] # updates blast data base with -m and shows sequenceserver where to find the data base

    sequenceserver -m # if added new .fasta files to folder and need to update data base

#### Did not use GMAP to identify Mu-promoter transcripts. Only Sequence Server BLAST database. See above.

**3. Align Trinity-reconstructed transcripts to the genome using GMAP**

Build a gmap index for the genome, required prior to aligning Trinity assembled transcripts, k-mer of 15 (default): 

    gmap_build -D /scratch.global/erika/rn20yh/rn20yh_fq -d w22.gmap ~/genome/Zmays_W22/10.fasta

**4. Running GMAP**

[gmapREADME](http://research-pub.gene.com/gmap/src/README)
[manpagesREADME](http://manpages.ubuntu.com/manpages/trusty/man1/gmap.1.html)
[manpages Ubuntu gmap_build](http://manpages.ubuntu.com/manpages/bionic/man1/gmap_build.1.html)
[msi gmap](https://www.msi.umn.edu/sw/gmap)

    module load samtools/1.9
    gmap -n 0 -D /home/springer/magnu513/rn20h/rn20h_00_fastq -d w22.gmap trinity.out.baf6021m1/Trinity.fasta -f gff3_gene > gmap.baf6021m1.gff3


Indexing a genome sorted BAM file allows one to quickly extract alignments overlapping particular genomic regions

    samtools index trinity.gmap.baf6021m1.bam

convert .bam file to .sam file for viewing - .bam file is compressed

    module load samtools
    samtools view -h trinity.gmap.baf6021m1.bam > trinity.gmap.baf6021m1.sam
    

