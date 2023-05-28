### Transcript abundance (CPM per Fragment) calculated for shared genomic sequence regions between wild-type W22 and mutant transcript assemblies. 

1. Identify transcribed regions of homology between control W22 wild-type and mutant assembled transcripts. Gene TSS partial and Mu TSS transcripts for each mutant and wild-type were aligned in Benchling using the annotated gene gDNA sequence as the template. The annotated gene cDNA transcript was aligned along with the transcript assemblies to analyze transcript splicing. Transcript `Zm_T001`were used as the canonical annotated transcript for most genes analyzed. For some genes, based on coverage of RNA-seq reads aligned to the reference genome, other isoforms of the gene `Zm_T002` or `Zm_T003` may have represented the wild-type transcript generated in the tissue sampled. 

	**Note:** Template is gDNA and if gDNA has -1 orientation reverse complement in benchling and input start index with annotated start
2. Regions of shared sequence between the wild-type transcript assembly and each mutant transcript assembly (gene TSS parital or Mu TSS) were identified. These regions included only gene sequence that was annotated in the canonical gene cDNA. Genomic coordinates of these shared sequence regions were recorded by using the gene gDNA sequence coordinates (template for the alignment).
3. Create a BED file of sequence regions shared by between W22 wild-type and mutant transccript assemblies for each gene by each transcript - gene TSS partial or Mu TSS. Used the gene.gff3 file as a template. In the original `gene.gff3` file the CDS overlaps with the exon (5'UTR or 3'UTR). Make sure genomic coordinates of shared regions for each gene do not overlap to avoid counting reads twice.

	**Example**
	
	For this HSF13, Zm00004b000433, Mu TSS transcript
	
	Orginal `gene.gff3` from `Zm-W22-REFERENCE-NRGENE-2.0_Zm00004b.1.gff3`
	
		chr1	ENSEMBLPEP	gene	13646839	13649038	.	+	.	ID=Zm00004b000433;Name=Zm00004b000433;biotype=protein_coding
		chr1	ENSEMBLPEP	mRNA	13646839	13649038	.	+	.	ID=Zm00004b000433_T001;Parent=Zm00004b000433;Name=Zm00004b000433_T001;biotype=protein_coding
		chr1	ENSEMBLPEP	intron	13647178	13647920	.	+	.	Parent=Zm00004b000433_T001;Name=intron.12580
		chr1	ENSEMBLPEP	exon	13646839	13647177	.	+	.	Parent=Zm00004b000433_T001;Name=exon.12581
		chr1	ENSEMBLPEP	exon	13647921	13649038	.	+	.	Parent=Zm00004b000433_T001;Name=exon.12582
		chr1	ENSEMBLPEP	CDS		13646839	13647177	.	+	0	Parent=Zm00004b000433_T001;Name=CDS.12583
		chr1	ENSEMBLPEP	CDS		13647921	13648769	.	+	0	Parent=Zm00004b000433_T001;Name=CDS.12584
	
	Formatted into BED format with shared sequence coordinates only
	
		chr01	13646839	13647177	exon.12581_CDS.12583	.	+
		chr01	13648770	13649038	exon.12582	.	+
		chr01	13647921	13648769	CDS.12584	.	+
		
	See `Table_S5_CPM-Fragment_raw.xlsx` in github repository for all alleles 
	
	**Optional:** Use [bedtools sort](https://bedtools.readthedocs.io/en/latest/content/tools/sort.html) to sort `gene-shared-sequence-regions.bed` file
	
		bedtools sort -i e2f13_Zm00004b023063_ct_CDS.bed >  #sorted gene.bed file	e2f13_Zm00004b023063_ct_CDS_sorted.bed
	
4. Convert RNA-seq BAM files from all W22 WT control and mutant biological replicates into BED files with [bedtools bamtobed](https://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html)
	
		bedtools bamtobed -i reads.bam > read.bed #bam to bed file conversion
		bedtools bamtobed -i E38_2.bam > ~/git/nf/rn21a/rnaseq/raw/21_bed/E38_2.bed
		
5. 	Use [bedtools intersect](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) to intersect the read aligned BED files (converted from BAM) with BED file of shared sequence regions.  

		bedtools intersect -a gene-shared-sequence-regions-sorted.bed -b reads-from-bam.bed -s -c > read_counts.bed # need to force strandedness with -s and -c counts		
**Note:** may need to reduce the reads-from-bam.bed file to fewer chromosomes because bedtools has a hard time counting all the reads at once 

		awk '$1 == "chr03" || $1 == "chr04" || $1 == "chr05"' C01.bam.bed > C01_ct.bed #to reduce the file size of the reads-from-bam.bed for W22 control samples used to compare multiple genes on multiple chromosomes, can pipe in an OR command || 
		
		bedtools intersect -a e2f13_Zm00004b023063_ct_CDS_sorted.bed -b C14_W04.bed -s -c > e2f13_C14_W04_count.bed 
		
	**Note:** only counting reads for one gene at a time, computational power to count reads from all genes in genome may not work with current pipeline (using bedtools). could think about feature counts or realign reads to new .gff3 file.
	
6. Calculate CPM for each read count/CDS or exon by using a scalar factor of sum raw read counts per library or biological replicate. 
7. Calculate and average CPM/Fragment for all biological replicates/Genotype for each shared CDS or exon
Step 6 and 7 completed in .R

```{r}library(dplyr)library(reshape2)library(tidyverse)library(tidyr)library(readr)

x = readRDS(file.path(dird, "rn21a.raw.rds")) tm = data.frame(x[["tm"]])th = data.frame(x[["th"]])

aa = left_join(tm, th)ab = aa %>% select(1:3,11:17) # filter columns. now I want to group by and sum how many reads there were in each library or replicateac = ab %>% group_by(SampleID, Tissue, Genotype, Treatment, Replicate) %>%  dplyr::summarise(ReadCount_sum = sum(ReadCount)) #calculate the read sum of each biological replicate for each sample (mutant or control). 

# edgeR calcNormFactors()# This function computes scaling factors to convert observed library sizes into effective library sizes. The effective library sizes for use in downstream analysis are lib.size * norm.factors where lib.size contains the original library sizes and norm.factors is the vector of scaling factors computed by this function.

# let edgeR calculate the normalization factors by creating matrix of ReadCounts. gid = rownames x SampleID = colnames

m1 = ab %>% select(1:3)
m1a = ab %>% pivot_wider(id_cols = gid, names_from = SampleID, values_from = ReadCount)m1a = column_to_rownames(m1a, var = "gid")m1b = as.matrix(m1a) 

# Considering that the formula is CPM = ((counts on the features) / library size) X 1,000,000
# The effective library size is then the original library size multiplied by the scaling factor.# count / (library size * normalization factor)# Then multiply that by a million to get CPM.

library(edgeR)en = calcNormFactors(m1b)
ens = colSums(m1b)
na = as.data.frame(en) %>% rownames_to_column(var = "SampleID") %>% dplyr::rename("norm.factors" = "en")nc = as.data.frame(ens) %>% rownames_to_column(var = "SampleID") %>% dplyr::rename("lib.size" = "ens")
a1 = left_join(na,nc)
a2 = a1 %>% mutate(lib.effect.size = lib.size*norm.factors) #library effect size calculated

# now caluclate CPM/Fragment by taking (Read Counts/Library Effect size) * 1000000 
a2 = read_tsv(file.path(dirr, "rn21a_libeffectsize.tsv")) # Library Effect size filemi = read_tsv(file.path(dirr, "cpm_cds.tsv")) # Read count per CDS shared regions filemh = left_join(mi, a2) #Join these filesmj = mh %>% mutate(cpm_cid = (count/(lib.size*norm.factors))*1000000) # calculate CPM/Fragment#now calculate the average CPM/Fragment for each shared CDS among biological replicates/Genotypeae = mj %>% separate_rows(allele, transcript, homology, sep = ",") af = ae %>% group_by(Genotype, cid) %>% mutate(avg_cpm_cid = mean(cpm_cid)) #group by Genotype and by CDS ID (cid)# now add column of mutant x control to pivot wider ai = ah %>% mutate(treatment = if_else(Genotype == "W22", "W22", "mutant"))aii = ai %>% select(-Genotype)#sum by transcript and treatment before pivoting wider aij = aii %>% group_by(tf,allele,transcript,homology,treatment) %>%  mutate(sum_avg_cpm_cid = sum(avg_cpm_cid)) #average CPM/CDS or Fragment for all shared Fragments/Transcript comparison (e.g. W22 x mutant Mu TSS)aij = ungroup(aij)aik = aij %>% group_by(tf,Tissue,allele,transcript,homology,treatment,chr,srd, sum_avg_cpm_cid) %>%   dplyr::summarise(cid = str_c(cid, collapse = ","), start = str_c(start, collapse = ","), end = str_c(end, collapse = ","),avg_cpm_cid = str_c(avg_cpm_cid, collapse = ","))la = ungroup(aik) %>% select(tf,Tissue,allele,transcript,homology,treatment,sum_avg_cpm_cid)aj = la %>% pivot_wider(names_from = treatment, values_from = sum_avg_cpm_cid)
```
See `Table_S5_CPM-Fragment_raw.xlsx` in github repository. Second Sheet includes Library Effect Size calculated