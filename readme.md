# This script is an overview of the quality control (QC) processes that were done for the ReD-Lat Exomes.
*Developed by Juliana Acosta-Uribe October 2022*

1.  **Set up environment and data**\
1.a Set up R\
1.b Set up data

2.  **General report of data quality and statistics previous to QC**\
2.a Autosomes quality report\
            - i. Extract autosome\
            - ii. Calculate quality metrics\
            - iii. Plot quality metrics\
            - iv. Create a record of your different QC metrics per sample (optional)\
2.b X and Y quality report\
            - i. Extract sex chromosomes\
            - ii. Get quality metrics for X and Y\
            - iii. Plot quality metrics for X\
            - iv. Plot quality metrics for Y\
            - v. Check chromosomal sex in plink

3: **Genotype quality control**\
3.a Autosomes\
3.b X Chromosome\
3.c Y Chromosome\




## 1. Set up environment and data for QC

### 1.a Set up R

It may be necessary to start R before, as you need to have knitr to run this markdown properly

```{r setup, echo = TRUE}
## Install and load required packages
#install.packages("psych")
library(psych)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("knitr")
library(knitr)
#install.packages("dplyr")
library(dplyr)

## Set your working directory
knitr::opts_chunk$set(root.dir = "/Users/acostauribe/ReDLat/JulianaQC/", echo = TRUE,
                      engine.path = list(plink = '~/bin/plink',
                                         king = '~/bin/king'))
## Set up the prefix for our dataset
PREFIX = 'ReDLat'
```

### 1.b Set up data

We will use **vcftools** and **bcftools** to clean our "starting dataset".

Since these are exomes, we first need to extract the 'exome sequencing targets' from the raw dataset. Targets are a '.bed' file downloaded from IDT website for xGen™ Exome Hybridization Panel v2 hg38 [the sequencing platform that was used]. This .bed file has three columns: *chr*, *start bp*, *end bp*. Additionally, we will remove individuals that are not part of the ReDLat project Exomes, but are in the same file.

**Note**: VCFtools output can be compressed directly into a .gz file by adding `–stdout | gzip -c > newname.gz`, but it wont generate a .log file.

Then, we will use bcftools to perform Left-alignment and normalization of indels, as well as a check to verify that the REF/ALT alleles are correct. However, it will NOT fix strand issues in your VCF (make sure VCF is properly aligned).


```{bash extract-target-align-reference-normalize, eval=FALSE, include=FALSE}
PREFIX=ReDLat
RAW_dataset='ADFTD.vqsr.snp.indel'
EXOME_targets='xgen-exome-research-panel-v2-targets-hg38.bed'
RAW_dataset_exclude_samples='RAW_samples_to_exclude.txt'
fasta_file='hg38.fa.gz'

# I. Extract targets 
vcftools --gzvcf $RAW_dataset.vcf.gz --bed $EXOME_targets --remove $RAW_dataset_exclude_samples --recode --recode-INFO-all --out $RAW_dataset.exome-targets

mv $RAW_dataset.exome-targets.recode.vcf $RAW_dataset.exome-targets.vcf
bgzip $RAW_dataset.exome-targets.vcf #Files for bcftools should be compressed
tabix -p vcf $RAW_dataset.exome-targets.vcf.gz #File should be indexed with Tabix

# II. Check reference allele and normalize INDELs
bcftools norm --check-ref ws --fasta-ref $fasta_file --output-type --output-type z $RAW_dataset.exome-targets.vcf.gz > $RAW_dataset.exome-targets.fasta-checked.vcf.gz
#--check-ref warn (w), exclude (x), or set/fix (s)
#--output-type compressed VCF (z) 
#fasta file was downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/
#The index file fasta.fai was created using http://www.htslib.org/doc/samtools-faidx.html

cp $RAW_dataset.exome-targets.fasta-checked.vcf.gz PREFIX.vcf.gz

# III. Save original files in a new directory
mkdir Raw_data
mv $RAW_dataset.* ./1.Raw_data
mv $RAW_dataset_exclude_samples ./1.Raw_data
mv $EXOME_targets ./1.Raw_data
```

## 2: General report of data quality and statistics previous to QC

For the following steps we are going to split the data into autosomes and the sex chromosomes. The reason we do this is because quality metrics like *Depth* or *Missingness* are very different between autosomes (where all individuals are diploid) and X and Y (where individuals can be diploid or haploid)

### 2.a Autosomes quality report
#### i. Extract autosomes
```{bash retain-autosomes, eval=FALSE, include=FALSE}
PREFIX=ReDLat

## Split into autosomes and sex chromosomes
## It is important to know how the chromosomes numbers are encoded in your vcf, e.g. they can be '1' or 'chr1'
vcftools --gzvcf $PREFIX.vcf.gz --chr chr1 --chr chr2 --chr chr3 --chr chr4 --chr chr5 --chr chr6 --chr chr7 --chr chr8 --chr chr9 --chr chr10 --chr chr11 --chr chr12 --chr chr13 --chr chr14 --chr chr15 --chr chr16 --chr chr17 --chr chr18 --chr chr19 --chr chr20 --chr chr21 --chr chr22 --recode --recode-INFO-all --out $PREFIX.autosomes

mv $PREFIX.autosomes.recode.vcf $PREFIX.autosomes.vcf
bgzip $PREFIX.autosomes.vcf
```

#### ii. Calculate quality metrics

```{bash autosome-preQC-stats, eval=FALSE, include=FALSE}
PREFIX=ReDLat

## I. Data description 
bcftools stats $PREFIX.autosomes.vcf.gz > $PREFIX.autosomes.vchk
# General distribution of depth, missingness, heterozygosity
vcftools --gzvcf $PREFIX.autosomes.vcf.gz --FILTER-summary --out $PREFIX.autosomes
# Generates a summary of the number of SNPs and Ts/Tv ratio for each FILTER category. 
# The output file has the suffix ".FILTER.summary"

## II. Individual missingness
vcftools --gzvcf $PREFIX.autosomes.vcf.gz --missing-indv --out $PREFIX.autosomes
# Generates a file reporting the missingness on a per-individual basis. 
# The output file has the suffix ".imiss".
# Individuals whose missingness is >10% should be added to 'flagged_samples' file

## III. Individual depth 
vcftools --gzvcf $PREFIX.autosomes.vcf.gz --depth --out $PREFIX.autosomes
# Generates a file containing the mean depth per individual. 
# This file has the suffix ".idepth".
# Individuals whose mean depth is <20 should be added to 'flagged_samples' file

## VI. Individual heterozygosity 
vcftools --gzvcf $PREFIX.autosomes.vcf.gz --het --out $PREFIX.autosomes
# Inbreeding coefficient, F, is estimated for each individual using a method of moments. 
# The resulting file has the suffix ".het"
# Individuals whose heterozygosity deviated more than 3 SD from the main should be added to 'flagged_samples' file

## V. Site missingness 
vcftools --gzvcf $PREFIX.autosomes.vcf.gz --missing-site --out $PREFIX.autosomes
# Generates a file reporting the missingness on a per-site basis. 
# The file has the suffix ".lmiss".

## VI. Site depth
vcftools --gzvcf $PREFIX.autosomes.vcf.gz --site-mean-depth --out $PREFIX.autosomes
# Generates a file containing the mean depth per site across all individuals. 
# This output file has the suffix ".ldepth.mean"
```

**Checkpoint:** At the end of the "General report" you should have at least the 7 following files: <file>.vchk, <file>.FILTER.summary, <file>.lmiss, <file>.ldepth, <file>.imiss, <file>.het and <file>.idepth

#### iii. Plot quality metrics

```{r autosome-preQC-stats-plots, echo=TRUE, fig.keep='all'}
# I. Sample based Metrics:
## 1. Missingness per sample
imiss = read.delim((paste0(PREFIX,".autosomes.imiss")), header = T, sep = "")

### Get some stats of missingness rate per sample (given as a fraction)
imiss_F_MISS = describe(imiss$F_MISS) #describe function from the "Psych" package gives basic stats
rownames(imiss_F_MISS) = c("sample_missingness_F")
### Individuals whose missingness is >10% should be identified
high_missingness = filter(imiss, F_MISS>0.1)
### Save it. It may come useful for debugging
write.table(high_missingness, 
            "samples_high_missingness_raw.txt",
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')

### Missingness counts per sample (given as an integer where N = variants)
imiss_N_MISS = describe(imiss$N_MISS)
rownames(imiss_N_MISS) = c("sample_missingness_N")

missingness_sample = bind_rows(imiss_F_MISS,
                               imiss_N_MISS)
                               
### Plot Missingness rate per sample
hist(imiss$F_MISS,
     xlab="Missingness rate",
     ylab="Samples", 
     main="Missingness rate per sample in - Autosomes preQC", 
     col="paleturquoise3")
```
![Autosomes/](https://github.com/acostauribe/exome_qc_tutorial/blob/main/Autosomes/missingness_sample_raw.png)            
```
## 2. Mean depth per sample
idepth = read.delim((paste0(PREFIX,".autosomes.idepth")), header = T, sep = "")

### Get some site depth stats:
depth_sample = describe(idepth$MEAN_DEPTH)
rownames(depth_sample) = c("sample_depth")
### Individuals whose mean depth is <20 should be identified
low_mean_depth = filter(idepth, MEAN_DEPTH<20)
write.table(low_mean_depth, 
            "samples_low_mean_depth_raw.txt",
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')

### Plot depth per sample
hist(idepth$MEAN_DEPTH,
     xlab="Mean Depth ",
     ylab="Samples", 
     main="Mean Depth per sample - Autosomes preQC", 
     col="paleturquoise3",
     breaks=50)
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/Autosomes/mean_depth_sample_raw.png)
            
```
## 3. Individual heterozygosity 

het = read.delim((paste0(PREFIX,".autosomes.het")), header = T, sep = "")

### Get some site depth stats:
heterozygosity_sample = describe(het$F)
rownames(heterozygosity_sample) = c("sample_heterozygosity_F")
### Identify the values for +3 and -3 standard deviations
heterozygosity_low_limit = mean(het$F)-(3*(sd(het$F)))
heterozygosity_high_limit = mean(het$F)+(3*(sd(het$F)))

### Plot heterozygosity per sample

hist(het$F,  
     freq=TRUE, 
     xlab="Heterozygosity F coefficient",  
     ylab="Samples", 
     main="Heterozygosity rate per sample - Autosomes preQC",
     col="paleturquoise3",
     breaks=50)
abline(v = (heterozygosity_low_limit), col="red")
abline(v = (heterozygosity_high_limit), col="red")
abline(v = (mean(het$F)), col="blue") 
legend("topleft",
       c("+/-3 SD","mean"),
       col=c("red","blue"),
       pch=16)

### Individuals whose heterozygosity deviated more than 3 SD from the mean should be identified
het_outlier_low = filter(het, F<heterozygosity_low_limit)
het_outlier_high = filter(het, F>heterozygosity_high_limit)
het_outlier_both = bind_rows(het_outlier_low,
                             het_outlier_high)
write.table(het_outlier_both, 
            "samples_heterozygosity_outliers_raw.txt",
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')
```     
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/Autosomes/samples_heterozygosity_outliers_raw.png)            
```
## 4. Generate a file with the descriptive statistics per sample
stats_sample = bind_rows(missingness_sample,
                         depth_sample,
                         heterozygosity_sample)
write.table(stats_sample,
            "stats_sample_autosomes.txt", 
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')
print(stats_sample)

## 5. Create a dataframe  with the samples that failed the autosome quality thresholds
high_missingness_id = select(high_missingness, INDV)
low_mean_depth_id = select(low_mean_depth, INDV)
het_outlier_low_id = select(het_outlier_low, INDV)
het_outlier_high_id = select(het_outlier_high, INDV)
flagged_samples_autosomes = bind_rows(high_missingness_id,
                                      low_mean_depth_id,
                                      het_outlier_low_id, 
                                      het_outlier_high_id)
# II. Site based metrics:

## 1. Missingness per site
lmiss = read.delim((paste0(PREFIX,".autosomes.lmiss")), header = T, sep = "")

### Get missingness rate per site stats (given as a fraction)
lmiss_F_MISS = describe(lmiss$F_MISS)
rownames(lmiss_F_MISS) = c("site_missingness_F")
### Missingness counts (given as an integer)
lmiss_N_MISS = describe(lmiss$N_MISS)
rownames(lmiss_N_MISS) = c("site_missingness_N")
missingness_site = bind_rows(lmiss_F_MISS,
                             lmiss_N_MISS)
### Plot Missingness per site
hist(lmiss$F_MISS,
        xlab="Missingness rate",
        ylab="Number of sites", 
        main="Missingness rate per site - Autosomes preQC", 
        col="paleturquoise3",
        breaks=50)
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/Autosomes/missingness_site_raw.png)         
```        
boxplot(lmiss$F_MISS,
        ylab="Missingness rate",
        xlab="Raw dataset", 
        main="Missingness rate per site - Autosomes preQC", 
        col="paleturquoise3")
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/Autosomes/missingness_site_boxplot_raw.png)                     
```            

## 2. Mean depth per site
ldepth.mean = read.delim((paste0(PREFIX,".autosomes.ldepth.mean")), header = T, sep = "")

### Get basic stats
depth_site = describe(ldepth.mean$MEAN_DEPTH)
rownames(depth_site) = c("site_mean_depth")

### Plot
hist(ldepth.mean$MEAN_DEPTH,
     xlab="Mean depth",
     ylab="Sites", 
     main="Mean depth per site - Autosomes preQC", 
     col="paleturquoise3")
            
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/Autosomes/mean_depth_site_raw.png)       
```            
            
### Take a zoom at the lower end
boxplot(ldepth.mean$MEAN_DEPTH,
        ylab="Mean depth per variant",
        xlab="Raw dataset", 
        main="Mean depth per site in Autosomes - preQC - 80X and lower",
        col = c("paleturquoise3"),
        ylim = c(0, 80))
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/Autosomes/mean_depth_site_raw_autosomes_80andlower.png)      
```
## 3. Generate a file with the descriptive statistics per site
stats_sites = bind_rows(missingness_site,
                         depth_site)
write.table(stats_sites,
            "stats_sites_autosomes.txt", 
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')
print(stats_sites)

## 4. VQSR (Variant Quality Score Recalibration) FILTER
filter_raw = read.delim((paste0(PREFIX,".autosomes.FILTER.summary")), header = T, sep = "")
print(filter_raw)

#Plot it
filter_raw %>%
  filter(!is.na(N_VARIANTS)) %>%
  arrange(N_VARIANTS) %>%
  mutate(FILTER=factor(FILTER, FILTER)) %>%
  ggplot( aes(x=FILTER, y=N_VARIANTS, label = N_VARIANTS) ) +
    geom_segment( aes(x=FILTER ,xend=FILTER, y=0, yend=N_VARIANTS), color="grey") +
    geom_point(size=3, color="#69b3a2") +
    geom_text(vjust=-1, size = 3) +
    coord_flip() +
    theme_minimal() +
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position="none") +
    scale_y_continuous(name = "Number of variants") +
    labs(title = "VQSR in Autosomes - preQC")
ggsave("VQSR-autosomes-preQC.png")  

```
#![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/Autosomes/VQSR-autosomes-preQC.png)        

Taking a close look to the descriptive statistics in the `stats_sample_autosomes.txt` and `stats_sites_autosomes.txt` files can help decide which are the best thresholds for QC.

#### iv. Create a record of your different QC metrics per sample (optional)

This step is absolutely optional, but I like to generate a dataframe where I annotate all samples with their different QC metrics. It is very useful when you come back and check the reason why a sample was dropped.

```{r qc-metrics-dataframe, echo=TRUE}
## I. Take the information of individuals from the .idepth file. 
## We are changing the names of the columns to specify these values come from raw data
sample_metrics = rename(idepth, n_sites_raw= N_SITES, mean_depth_raw = MEAN_DEPTH)

## II. Using the "match" function, we will create a new column in 'sample_metrics' with the missingness per sample.
## Technically, imiss and idepth should have the same samples in the same order, but using the "match" function will be useful when we start dropping samples
sample_metrics$missingness_raw = imiss$F_MISS[match(sample_metrics$INDV, imiss$INDV)]

## III. Using the "match" function, we will create a new column in 'sample_metrics' with the heterozygosity per sample 
sample_metrics$heterozygosity_F_raw = het$F[match(sample_metrics$INDV, het$INDV)]

## IV. Save as a file (optional)
write.table(sample_metrics,
            "sample_metrics_autosomes_preQC.txt", 
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE)
```

We will add additional columns in the X and Y quality report before saving `sample_metrics` as a file

### 2.b X and Y quality report

Since this data set is an Exome, we wanted to explore the distribution of markers prior to the quality control and determination of chromosomal sex. You will need a dataframe/table with the reported phenotypic sex of the included samples `sample_sex.txt` sex should be encoded: 1=male, 2=female

#### i. Extract sex chromosomes

Before extracting the sex chromosomes check the file to see how are these encoded, a good idea is to do `grep X` or `grep x`. For example, in the ReDLat file, these are encoded as 'chrX' and 'chrY'
```{bash retain-sex-chromosomes, eval=FALSE, include=FALSE}
PREFIX=ReDLat
for chr in X Y
do
vcftools --gzvcf $PREFIX.vcf.gz --chr chr$chr --recode --recode-INFO-all --out $PREFIX.$chr
mv $PREFIX.$chr.recode.vcf $PREFIX.$chr.vcf
bgzip $PREFIX.$chr.vcf
done
```

#### ii. Get quality metrics for X and Y

It is useful to know the number of expected males and females

```{r count-males-females, echo=TRUE}
## Import a dataframe with disclosed sex of samples
## sex should be encoded 1=male, 2=female
sample_sex = read.delim(file = "sample_sex.txt", header=T)

## Count number of females
females = sum(sample_sex$sex == 2) 
print(paste0("Number of reported females in dataset: ", females))

## Count number of males
males = sum(sample_sex$sex == 1) 
print(paste0("Number of reported males in dataset: ", males))

## Annotate the sample_metrics dataframe 
sample_metrics$reported_sex = sample_sex$sex[match(sample_metrics$INDV, sample_sex$sample)]
```

Calculate quality statistics

```{bash XY-preQC-check, eval=FALSE, include=FALSE}
PREFIX=ReDLat

# Make a list of only males 
awk '{ if($2 == 1) {print $1}}' sample_sex.txt > male.samples.txt
# Make a list of only females 
awk '{ if($2 == 2) {print $1}}' sample_sex.txt > female.samples.txt

# These metrics are discussed in section 2.a
for chr in X Y 
do
bcftools stats $PREFIX.$chr.vcf.gz > $PREFIX.$chr.vchk
vcftools --gzvcf $PREFIX.$chr.vcf.gz --FILTER-summary --out $PREFIX.$chr
vcftools --gzvcf $PREFIX.$chr.vcf.gz --missing-indv --out $PREFIX.$chr
vcftools --gzvcf $PREFIX.$chr.vcf.gz --depth --out $PREFIX.$chr
  # For site based metrics we will split the analyses into male and female
  for sex in female male
  do
  vcftools --gzvcf $PREFIX.$chr.vcf.gz --keep $sex.samples.txt --missing-site --out $PREFIX.$chr.$sex
  vcftools --gzvcf $PREFIX.$chr.vcf.gz --keep $sex.samples.txt --site-mean-depth --out $PREFIX.$chr.$sex
  done
done

# You should end with 14 files. For each chromosome ( X and Y) you get:
## I. Sample based:
### 1. missingness per sample (.imiss)
### 2. depth per sample (.idepth)
## II. Site based
### 1. Missingness per site 
#### a. Females (female.lmiss)
#### b. Males (male.lmiss)
### 2. Mean depth per site
#### a. Females (female.ldepth.mean)
#### b. Males (female.ldepth.mean)
### 3. VQSR quality 
## III. General stats (.vchk)
```

#### iii. Plot quality metrics for X

Plot quality metrics of all samples, females and males. This script will also generate a txt file with the descriptive statistics

```{r X-preQC-check-plots, echo=TRUE}
# I.  Sample based:
## 1. Missingness per sample
imiss_X = read.delim((paste0(PREFIX,".X.imiss")), header = T, sep = "")
imiss_X$sex = sample_sex$sex[match(imiss_X$INDV,sample_sex$sample)]
imiss_X_female = filter(imiss_X, sex == 2)
imiss_X_male = filter(imiss_X, sex == 1)

## Get some site missingness stats:
## Missingness rate (given as a fraction)
imiss_X_F_MISS_all = describe(imiss_X$F_MISS)
rownames(imiss_X_F_MISS_all) = c("sample_F_missingness_X_all")
imiss_X_F_MISS_female = describe(imiss_X_female$F_MISS)
rownames(imiss_X_F_MISS_female) = c("sample_F_missingness_X_female")
imiss_X_F_MISS_male = describe(imiss_X_male$F_MISS)
rownames(imiss_X_F_MISS_male) = c("sample_F_missingness_X_male")

## Missingness counts (given as an integer where N = samples x 2)
imiss_X_N_MISS_all = describe(imiss_X$N_MISS)
rownames(imiss_X_N_MISS_all) = c("sample_N_missingness_X_all")
imiss_X_N_MISS_female = describe(imiss_X_female$N_MISS)
rownames(imiss_X_N_MISS_female) = c("sample_N_missingness_X_female")
imiss_X_N_MISS_male = describe(imiss_X_male$N_MISS)
rownames(imiss_X_N_MISS_male) = c("sample_N_missingness_X_male")

missingness_sample_X_chromosome = bind_rows(imiss_X_F_MISS_all,
                                          imiss_X_F_MISS_female,
                                          imiss_X_F_MISS_male,
                                          imiss_X_N_MISS_all,
                                          imiss_X_N_MISS_female,
                                          imiss_X_N_MISS_male)

## Plot Missingness per sample
## Missingness counts helps to see distribution
boxplot(imiss_X$F_MISS, imiss_X_female$F_MISS, imiss_X_male$F_MISS,
        ylab="Missingness",
        xlab="Raw dataset", 
        main="Missingness rate per sample in X chromosome - preQC",
        sub = "Distribution according to disclosed sex",
        col = c("lavender","paleturquoise3", "lightgoldenrod1"),
        names = c("all", "Female", "Male"))
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/chrX/missingness_sample_raw_Xchr.png)
```  
## 2. Mean depth per sample
idepth_X = read.delim((paste0(PREFIX,".X.idepth")), header = T, sep = "")
idepth_X$sex = sample_sex$sex[match(idepth_X$INDV,sample_sex$sample)]
idepth_X_female = filter(idepth_X, sex == 2)
idepth_X_male = filter(idepth_X, sex == 1)

## Get some site depth stats:
## Mean depth per sample
idepth_X_stats_all = describe(idepth_X$MEAN_DEPTH)
rownames(idepth_X_stats_all) = c("sample_depth_X_all")
idepth_X_stats_female = describe(idepth_X_female$MEAN_DEPTH)
rownames(idepth_X_stats_female) = c("sample_depth_X_female")
idepth_X_stats_male = describe(idepth_X_male$MEAN_DEPTH)
rownames(idepth_X_stats_male) = c("sample_depth_X_male")

depth_sample_X_chromosome = bind_rows(idepth_X_stats_all,
                                      idepth_X_stats_female,
                                      idepth_X_stats_male)

## Plot depth per sample
boxplot(idepth_X$MEAN_DEPTH, idepth_X_female$MEAN_DEPTH, idepth_X_male$MEAN_DEPTH,
        ylab="Mean sample depth",
        xlab="Raw dataset", 
        main="Mean sample depth in X chromosome - preQC",
        sub = "Distribution according to disclosed sex",
        col = c("lavender","paleturquoise3", "lightgoldenrod1"),
        names = c("all", "Female", "Male"))
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/chrX/mean_depth_sample_raw_chrX.png)
```
## 3. Generate a file with the descriptive statistics per sample
stats_sample_X_chromosome = bind_rows(missingness_sample_X_chromosome,
                                      depth_sample_X_chromosome)

write.table(stats_sample_X_chromosome,
            "stats_sample_X_chromosome.txt", 
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')

print(stats_sample_X_chromosome)

# II. Site based:
### 1. Missingness per site

### a. Females (female.lmiss)
lmiss_X_female = read.csv((paste0(PREFIX,".X.female.lmiss")), header = T, sep = "")
### b. Males (male.lmiss)
lmiss_X_male = read.csv((paste0(PREFIX,".X.male.lmiss")), header = T, sep = "")

## Get some site missingness stats:
## Missingness rate (given as a fraction)
lmiss_X_female_F_MISS = describe(lmiss_X_female$F_MISS)
rownames(lmiss_X_female_F_MISS) = c("site_F_missingness_X_female")
lmiss_X_male_F_MISS = describe(lmiss_X_male$F_MISS)
rownames(lmiss_X_male_F_MISS) = c("site_F_missingness_X_male")

## Missingness counts (given as an integer
lmiss_X_female_N_MISS = describe(lmiss_X_female$N_MISS)
rownames(lmiss_X_female_N_MISS) = c("site_N_missingness_X_female")
lmiss_X_male_N_MISS = describe(lmiss_X_male$N_MISS)
rownames(lmiss_X_male_N_MISS) = c("site_N_missingness_X_male")

missingness_site_X_chromosome = bind_rows(lmiss_X_female_F_MISS,
                                     lmiss_X_male_F_MISS,
                                     lmiss_X_female_N_MISS,
                                     lmiss_X_male_N_MISS)

## Plot Missingness per site
## Missingness counts helps to see distribution
boxplot(lmiss_X_female$N_MISS, lmiss_X_male$N_MISS,
        ylab="Missingness counts",
        xlab="Raw dataset", 
        main="Missingness counts per site in X chromosome - preQC",
        sub = "Distribution according to disclosed sex",
        col = c("paleturquoise3", "lightgoldenrod1"),
        names = c("Female", "Male"))
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/chrX/missingness-counts_site_boxplot_raw_chrX.png)
```            
## Take a closer look at the missingness RATE per site
boxplot(lmiss_X_female$F_MISS, lmiss_X_male$F_MISS,
        ylab="Missingness rate",
        xlab="Raw dataset", 
        main="Missingness rate per site in X chromosome - preQC",
        sub = "Distribution according to disclosed sex",
        col = c("paleturquoise3", "lightgoldenrod1"),
        names = c("Female", "Male"),
        ylim = c(0, 0.05))
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/chrX/missingness-rate_site_boxplot_raw_chrX.png)
```            
### 2. Mean depth per site
#### a. Females (female.ldepth.mean)
ldepth_X_female = read.csv((paste0(PREFIX,".X.female.ldepth.mean")), header = T, sep = "")
#### b. Males (female.ldepth.mean)
ldepth_X_male = read.csv((paste0(PREFIX,".X.male.ldepth.mean")), header = T, sep = "")

## Get some stats
ldepth_X_female_depth = describe(ldepth_X_female$MEAN_DEPTH)
rownames(ldepth_X_female_depth ) = c("site_mean-depth_X_female")
ldepth_X_male_depth = describe(ldepth_X_male$MEAN_DEPTH)
rownames(ldepth_X_male_depth ) = c("site_mean-depth_X_male")

mean_site_depth_X_chromosome = bind_rows(ldepth_X_female_depth,
                                         ldepth_X_male_depth)

## Plot as boxplots
boxplot(ldepth_X_female$MEAN_DEPTH, ldepth_X_male$MEAN_DEPTH,
        ylab="Mean depth per site",
        xlab="Raw dataset", 
        main="Mean depth per site in X chromosome - preQC",
        sub = "Distribution according to disclosed sex",
        col = c("paleturquoise3", "lightgoldenrod1"),
        names = c("Female", "Male"))
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/chrX/mean_depth_sample_raw_chrX.png)
```            
## Take a zoom at the lower end
boxplot(ldepth_X_female$MEAN_DEPTH, ldepth_X_male$MEAN_DEPTH,
        ylab="Mean depth per variant",
        xlab="Raw dataset", 
        main="Mean depth per site in X chromosome - preQC - 80X and lower",
        sub = "Distribution according to disclosed sex",
        col = c("paleturquoise3", "lightgoldenrod1"),
        names = c("Female", "Male"),
        ylim = c(0, 80))
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/chrX/mean_depth_sample_raw_80xandlower_chrX.png) 
```            
## 3. Generate a file with the descriptive statistics per site

stats_site_X_chromosome = bind_rows(missingness_site_X_chromosome,
                                    mean_site_depth_X_chromosome)
write.table(stats_site_X_chromosome,
            "stats_site_X_chromosome.txt", 
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')
print(stats_site_X_chromosome)

## 4. VQSR (Variant Quality Score Recalibration) FILTER
filter_X_raw = read.delim((paste0(PREFIX,".X.FILTER.summary")), header = T, sep = "")
print(filter_X_raw)

#Plot it
filter_X_raw %>%
  filter(!is.na(N_VARIANTS)) %>%
  arrange(N_VARIANTS) %>%
  mutate(FILTER=factor(FILTER, FILTER)) %>%
  ggplot( aes(x=FILTER, y=N_VARIANTS, label = N_VARIANTS) ) +
    geom_segment( aes(x=FILTER ,xend=FILTER, y=0, yend=N_VARIANTS), color="grey") +
    geom_point(size=3, color="#69b3a2") +
    geom_text(vjust=-1, size = 3) +
    coord_flip() +
    theme_minimal() +
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position="none") +
    scale_y_continuous(name = "Number of variants") +
    labs(title = "VQSR in X chromosome - preQC")
ggsave("VQSR-X-preQC.png")  

```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/chrX/VQSR-X-preQC.png)            
```
In the ReDLat dataset you can observe that all samples (regardless chromosomal sex) have a mean sample depth above 20X and the missingness rate per sample for X is below 1%.

#### iii. Plot quality metrics for Y

```{r Y-preQC-check-plots, echo=TRUE}
#### These metrics use almost the same code as X-preQC-check-plots

# I.  Sample based:
## 1. Missingness per sample
imiss_Y = read.delim((paste0(PREFIX,".Y.imiss")), header = T, sep = "")
imiss_Y$sex = sample_sex$sex[match(imiss_Y$INDV,sample_sex$sample)]
imiss_Y_female = filter(imiss_Y, sex == 2)
imiss_Y_male = filter(imiss_Y, sex == 1)

## Get some site missingness stats:
## Missingness rate (given as a fraction)
imiss_Y_F_MISS_female = describe(imiss_Y_female$F_MISS)
rownames(imiss_Y_F_MISS_female) = c("sample_F_missingness_Y_female")
imiss_Y_F_MISS_male = describe(imiss_Y_male$F_MISS)
rownames(imiss_Y_F_MISS_male) = c("sample_F_missingness_Y_male")

## Missingness counts (given as an integer)
imiss_Y_N_MISS_female = describe(imiss_Y_female$N_MISS)
rownames(imiss_Y_N_MISS_female) = c("sample_N_missingness_Y_female")
imiss_Y_N_MISS_male = describe(imiss_Y_male$N_MISS)
rownames(imiss_Y_N_MISS_male) = c("sample_N_missingness_Y_male")

missingness_sample_Y_chromosome = bind_rows(imiss_Y_F_MISS_female,
                                          imiss_Y_F_MISS_male,
                                          imiss_Y_N_MISS_female,
                                          imiss_Y_N_MISS_male)

## Plot Missingness per sample
## Missingness counts helps to see distribution
boxplot(imiss_Y_female$F_MISS, imiss_Y_male$F_MISS,
        ylab="Missingness",
        xlab="Raw dataset", 
        main="Missingness rate per sample in Y chromosome - preQC",
        sub = "Distribution according to disclosed sex",
        col = c("paleturquoise3", "lightgoldenrod1"),
        names = c("Female", "Male"))
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/chrY/missingness_sample_raw_Ychr.png)
```
## 2. Mean depth per sample
idepth_Y = read.delim((paste0(PREFIX,".Y.idepth")), header = T, sep = "")
idepth_Y$sex = sample_sex$sex[match(idepth_Y$INDV,sample_sex$sample)]
idepth_Y_female = filter(idepth_Y, sex == 2)
idepth_Y_male = filter(idepth_Y, sex == 1)

## Get some site depth stats:
## Mean depth per sample
idepth_Y_stats_female = describe(idepth_Y_female$MEAN_DEPTH)
rownames(idepth_Y_stats_female) = c("sample_depth_Y_female")
idepth_Y_stats_male = describe(idepth_Y_male$MEAN_DEPTH)
rownames(idepth_Y_stats_male) = c("sample_depth_Y_male")

depth_sample_Y_chromosome = bind_rows(idepth_Y_stats_female,
                                      idepth_Y_stats_male)

## Plot depth per sample
boxplot(idepth_Y_female$MEAN_DEPTH, idepth_Y_male$MEAN_DEPTH,
        ylab="Mean sample depth",
        xlab="Raw dataset", 
        main="Mean sample depth in Y chromosome - preQC",
        sub = "Distribution according to disclosed sex",
        col = c("paleturquoise3", "lightgoldenrod1"),
        names = c("Female", "Male"))
```
![image](  https://github.com/acostauribe/exome_qc_tutorial/blob/main/chrY/mean_depth_sample_raw_chrY.png)         
```            
## 3. Generate a file with the descriptive statistics per sample
stats_sample_Y_chromosome = bind_rows(missingness_sample_Y_chromosome,
                                      depth_sample_Y_chromosome)
write.table(stats_sample_Y_chromosome,
            "stats_sample_Y_chromosome.txt", 
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')
print(stats_sample_Y_chromosome)

# II. Site based:
### 1. Missingness per site

### a. Females (female.lmiss)
lmiss_Y_female = read.csv((paste0(PREFIX,".Y.female.lmiss")), header = T, sep = "")
### b. Males (male.lmiss)
lmiss_Y_male = read.csv((paste0(PREFIX,".Y.male.lmiss")), header = T, sep = "")

## Get some site missingness stats:
## Missingness rate (given as a fraction of samples missing the variant)
lmiss_Y_female_F_MISS = describe(lmiss_Y_female$F_MISS)
rownames(lmiss_Y_female_F_MISS) = c("site_F_missingness_Y_female")
lmiss_Y_male_F_MISS = describe(lmiss_Y_male$F_MISS)
rownames(lmiss_Y_male_F_MISS) = c("site_F_missingness_Y_male")

## Missingness counts (given as an integer where N = samples x 2)
lmiss_Y_female_N_MISS = describe(lmiss_Y_female$N_MISS)
rownames(lmiss_Y_female_N_MISS) = c("site_N_missingness_Y_female")
lmiss_Y_male_N_MISS = describe(lmiss_Y_male$N_MISS)
rownames(lmiss_Y_male_N_MISS) = c("site_N_missingness_Y_male")
missingness_site_Y_chromosome = bind_rows(lmiss_Y_female_F_MISS,
                                     lmiss_Y_male_F_MISS,
                                     lmiss_Y_female_N_MISS,
                                     lmiss_Y_male_N_MISS)

## Plot Missingness per site
## Missingness counts helps to see distribution
boxplot(lmiss_Y_female$N_MISS, lmiss_Y_male$N_MISS,
        ylab="Missingness counts",
        xlab="Raw dataset", 
        main="Missingness counts per site in Y chromosome - preQC",
        sub = "Distribution according to disclosed sex",
        col = c("paleturquoise3", "lightgoldenrod1"),
        names = c("Female", "Male"))
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/chrY/missingness-counts_site_boxplot_raw_chrY.png)
```           
## Take a look at missingness rate
boxplot(lmiss_Y_female$F_MISS, lmiss_Y_male$F_MISS,
        ylab="Missingness rate",
        xlab="Raw dataset", 
        main="Missingness rate per site in Y chromosome - preQC",
        sub = "Distribution according to disclosed sex",
        col = c("paleturquoise3", "lightgoldenrod1"),
        names = c("Female", "Male"))
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/chrY/missingness-rate_site_boxplot_raw_chrY.png)
```
### 2. Mean depth per site
#### a. Females (female.ldepth.mean)
ldepth_Y_female = read.csv((paste0(PREFIX,".Y.female.ldepth.mean")), header = T, sep = "")
#### b. Males (female.ldepth.mean)
ldepth_Y_male = read.csv((paste0(PREFIX,".Y.male.ldepth.mean")), header = T, sep = "")

## Get some stats
ldepth_Y_female_depth = describe(ldepth_Y_female$MEAN_DEPTH)
rownames(ldepth_Y_female_depth ) = c("site_mean-depth_Y_female")
ldepth_Y_male_depth = describe(ldepth_Y_male$MEAN_DEPTH)
rownames(ldepth_Y_male_depth ) = c("site_mean-depth_Y_male")

mean_site_depth_Y_chromosome = bind_rows(ldepth_Y_female_depth,
                                         ldepth_Y_male_depth)

## Plot as boxplots
boxplot(ldepth_Y_female$MEAN_DEPTH, ldepth_Y_male$MEAN_DEPTH,
        ylab="Mean depth per site",
        xlab="Raw dataset", 
        main="Mean depth per site in Y chromosome - preQC",
        sub = "Distribution according to disclosed sex",
        col = c("paleturquoise3", "lightgoldenrod1"),
        names = c("Female", "Male"))
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/chrY/mean_depth_site_raw_chrY.png)         
```            
## 3. Generate a file with the descriptive statistics per sample
stats_site_Y_chromosome = bind_rows(missingness_site_Y_chromosome,
                                      mean_site_depth_Y_chromosome)
write.table(stats_site_Y_chromosome,
            "stats_site_Y_chromosome.txt", 
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')
print(stats_sample_Y_chromosome)

## 4. VQSR (Variant Quality Score Recalibration) FILTER
filter_Y_raw = read.delim((paste0(PREFIX,".Y.FILTER.summary")), header = T, sep = "")
print(filter_Y_raw)

#Plot it
filter_Y_raw %>%
  filter(!is.na(N_VARIANTS)) %>%
  arrange(N_VARIANTS) %>%
  mutate(FILTER=factor(FILTER, FILTER)) %>%
  ggplot( aes(x=FILTER, y=N_VARIANTS, label = N_VARIANTS) ) +
    geom_segment( aes(x=FILTER ,xend=FILTER, y=0, yend=N_VARIANTS), color="grey") +
    geom_point(size=3, color="#69b3a2") +
    geom_text(vjust=-1, size = 3) +
    coord_flip() +
    theme_minimal() +
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position="none") +
    scale_y_continuous(name = "Number of variants") +
    labs(title = "VQSR in Y chromosome - preQC")
ggsave("VQSR-Y-preQC.png")  
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/chrY/VQSR-Y-preQC.png)
```
Y chromosome sample statistics can also give us a hint of chromosomal vs disclosed/phenotypic sex discordances. Notice how there are only 83 markers in the y chromosome exome!

#### v. Check chromosomal sex in plink

Since this is exome data, and the exomic portion of the Y chromosome is so small, we will can only use plink to check for sex using the X chromosome.

```{bash check-sex-plink, eval=FALSE, include=FALSE}
PREFIX=ReDLat
sample_sex=sample_sex.txt

## I. Import files into plink and add sex information
# more info in https://www.cog-genomics.org/plink/1.9/input#vcf
plink --vcf $PREFIX.X.vcf.gz --keep-allele-order --vcf-half-call h --double-id --make-bed --out $PREFIX.X.plink

awk '{print $1, $1, $2}' $sample_sex > sex_plink.txt

plink --bfile $PREFIX.X.plink --update-sex sex_plink.txt --make-bed --out $PREFIX.X.plink.sex

## II.Remove X pseudoautosomal region (if your data is hg19, use 'hg19')
plink --bfile $PREFIX.X.plink.sex --split-x hg38 --make-bed --out $PREFIX.X.plink.sex.split-x

## III. Check if variants have an id in the second column of the bim file.
# Assign IDs to variants if they dont have it (maximum lenght of variant IDs is 20, which means longer insertions will be left without a variant ID)
plink --bfile $PREFIX.X.plink.sex.split-x --set-missing-var-ids '@:#' --new-id-max-allele-len 20 --make-bed --out $PREFIX.X.plink.sex.split-x.id

## IV. Prune for Linkage Disequilibrium 
#(make sure that the variants have an ID in the bim file)
plink  --bfile $PREFIX.X.plink.sex.split-x.id --indep-pairphase 20000 2000 0.5
# this produces plink.prune.in and plink.prune.out

## V. Retain independent markers (in linkage equilibrium)
plink  --bfile $PREFIX.X.plink.sex.split-x.id --extract plink.prune.in --make-bed --out $PREFIX.X.plink.sex.split-x.id.LD

## VI. Check sex 
plink --bfile $PREFIX.X.plink.sex.split-x.id.LD --check-sex 0.3 0.7 --out $PREFIX.X.plink.sex.split-x.id.LD.Xsex
# more info on https://www.cog-genomics.org/plink/1.9/basic_stats#check_sex

## VII. Save the files in a directory of their own
#mkdir 3.Check_sex
#mv $PREFIX.* ./3.Check_sex/
#mv ./3.Check_sex/$PREFIX.vcf.gz ./
```

Take a look at the results

```{r plot-check-sex-plink, eval=FALSE, include=FALSE}
#install.packages("ggbeeswarm")
library(ggbeeswarm)
#install.packages("scales")
library(scales)

## Load X chromosome F coefficients calculated by Plink
X_file = read.delim(file = paste0(PREFIX,".X.plink.sex.split-x.id.LD.Xsex.sexcheck"), 
                    header = TRUE, 
                    sep = '')
## Plot these coefficients comparing males vs. females
## Roughly you expect Female to have an F coefficient < 0.2-0.3 and males have an F coefficient > 0.7-0.8
ggplot(X_file, 
       aes(x = factor(PEDSEX,
                      labels = c("Male", "Female")),
           y = F,
           color = PEDSEX)) +
    geom_quasirandom(alpha = 0.7,
                   size = 1.5) + 
    labs(title = "Chromosomal sex assignement in samples",
         x = "Disclosed sex",
         y = "F coefficient X chromosome") +
  theme_minimal() +
  theme(legend.position = "none") +
  #geom_hline(aes(yintercept = 0.3)) + 
  #geom_hline(aes(yintercept = 0.7)) + 
  ggsave("Chromosomal sex assignement in samples.png", width = 8, height = 5) 
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/chrX/Chromosomal%20sex%20assignement%20in%20samples.png)

We can add the missingness of Y chromosome as another column to select the samples that will be taken out. Individuals that fail sex check should be added to 'flagged_samples' file

```{r determine-sex-check-fails}
X_file = read.delim(file = paste0(PREFIX,".X.plink.sex.split-x.id.LD.Xsex.sexcheck"), 
                    header = TRUE, 
                    sep = '')

## Create a new column in 'sample_metrics' with the X chromosome F coefficients
sample_metrics$X_Fcoeff = X_file$F[match(sample_metrics$INDV, X_file$IID)]
sample_metrics$X_STATUS = X_file$STATUS[match(sample_metrics$INDV, X_file$IID)]
sample_metrics$Y_F_MISS = imiss_Y$F_MISS[match(sample_metrics$INDV, imiss_Y$INDV)]

## Generate a file
write.table(sample_metrics,
            "sample_metrics_autosomes_sex-check.preQC.txt", 
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE,
            sep = '\t')

## I am taking a combination of X F coefficient and missingness in Y chromosome to filter samples that DEFINITELY fail check sex
sex_fail = filter(sample_metrics, 
                        (reported_sex == "1" & X_Fcoeff < 0.3 & Y_F_MISS > 0.1 |
                        reported_sex == "2" & X_Fcoeff > 0.7 & Y_F_MISS < 0.7)) 
 
## Samples that were flagged by plink but dont meet the Y chromosome condition are labeled as Warning
sex_warning = filter(sample_metrics, 
                        (reported_sex == "1" & X_Fcoeff < 0.7 & Y_F_MISS < 0.1 |
                        reported_sex == "2" & X_Fcoeff > 0.3 & Y_F_MISS > 0.7)) 

## Add the sex-fail samples to the flagged sample list
sex_fail_id = select(sex_fail, INDV)

flagged_samples_autosomes_XY = bind_rows(flagged_samples_autosomes,
                                         sex_fail_id)
write.table(flagged_samples_autosomes_XY,
            "flagged_samples_preQC.txt", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)
```
All the samples of `sex_fail` and `sex_warning` should be double checked with the submitting center to determine if there was a sex labeling error.


## 3: Genotype quality control
### 3.a Autosomes 

Retain genotypes with Depth >= 20 & Quality >= 20 and check how this affected the missingness rates
```{bash autosomes-genotype-QC, eval=FALSE, include=FALSE}
PREFIX=ReDLat
DP=20
GQ=20

## I. Retain genotypes with Depth >= 20 & Quality >= 20
# Other genotypes will be set to missing
vcftools --gzvcf $PREFIX.autosomes.vcf.gz --minDP $DP --minGQ $GQ --recode --recode-INFO-all --out $PREFIX.autosomes.DP$DP.GQ$GQ
#--minDP <float> Includes only genotypes greater than or equal to the "--minDP"value. DP <30 are set as missing
#This option requires that the "DP" FORMAT tag is specified for all sites.
#--minGQ <float> Excludes all genotypes with a quality below the threshold specified. 
#This option requires that the "GQ" FORMAT tag is specified for all sites.

mv $PREFIX.autosomes.DP$DP.GQ$GQ.recode.vcf $PREFIX.autosomes.DP$DP.GQ$GQ.vcf
bgzip $PREFIX.autosomes.DP$DP.GQ$GQ.vcf

## II. Check how this affected the missingness rates
# Per site
vcftools --gzvcf $PREFIX.autosomes.DP$DP.GQ$GQ.vcf.gz --missing-site --out $PREFIX.autosomes.DP$DP.GQ$GQ
## Per sample
vcftools --gzvcf $PREFIX.autosomes.DP$DP.GQ$GQ.vcf.gz --missing-indv --out $PREFIX.autosomes.DP$DP.GQ$GQ
```

Plot missingness rates per site and per sample comparing raw dataset and after genotype quality control.
Declare the same GQ and DP values.
Notice that since we are retaining genotypes with Depth >= 20, low depth is not a problem any more
            
```{r compare-missingness-rates-autosomes}
DP_autosomes = 20
GQ_autosomes = 20
# We had previously generated an .lmiss and .imiss files of raw dataset
## 1. Missingness per site
lmiss_GT = read.delim((paste0(PREFIX,".autosomes.DP",DP_autosomes,".GQ",GQ_autosomes,".lmiss")), header = T, sep = "")
lmiss_GT_F_MISS = describe(lmiss$F_MISS)
rownames(lmiss_GT_F_MISS) = c("site_missingness_F_GT")

boxplot(lmiss$F_MISS, lmiss_GT$F_MISS,
        ylab ="Missingness rate",
        names = c("Raw dataset", "After genotype QC"),
        main = "Missingness rate per site - Autosomes genotype QC", 
        col = c("paleturquoise3", "lavender"))
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/Autosomes-GT/missingness-rate_site_boxplot_GT-QC.png)            
```
boxplot(lmiss$F_MISS, lmiss_GT$F_MISS,
        ylab ="Missingness rate",
        names = c("Raw dataset", "After genotype QC"),
        main = "Missingness rate per site - Autosomes genotype QC <0.01", 
        col = c("paleturquoise3", "lavender"),
        ylim = c(0,0.01))
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/Autosomes-GT/missingness-rate_site_boxplot_GT-QC_lowend.png)            
```
imiss_GT = read.delim((paste0(PREFIX,".autosomes.DP",DP_autosomes,".GQ",GQ_autosomes,".imiss")), header = T, sep = "")
imiss_GT_F_MISS = describe(imiss$F_MISS) 
rownames(imiss_GT_F_MISS) = c("sample_missingness_F_GT")

boxplot(imiss$F_MISS, imiss_GT$F_MISS,
        ylab = "Missingness rate",
        names = c("Raw dataset", "After genotype QC"),
        main = "Missingness rate per sample - Autosomes genotype QC", 
        col = c("paleturquoise3", "lavender"))
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/Autosomes-GT/missingness-rate_sample_boxplot_GT-QC.png)            
```
stats_missingness_autosomeGT = bind_rows(lmiss_GT_F_MISS,
                                         imiss_GT_F_MISS)
write.table(stats_missingness_autosomeGT,
            "stats_missingness_autosomeGT.txt", 
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')
print(stats_sample_Y_chromosome)

high_missingness_GT = filter(imiss_GT, F_MISS>0.1)

# Annotate the sample_metrics file with the new Missingness values
sample_metrics$missingness_GT = imiss_GT$F_MISS[match(sample_metrics$INDV, imiss$INDV)]
```

### 3.b X Chromosome

When we checked the missingness rates and depth rates in x chromosomes, the mean values were similar to those of autosomes (both for male and female samples). We tried two settings. 1. DP10 GQ20 2. DP20 GQ20
```{bash X-genotype-QC, eval=FALSE, include=FALSE}
PREFIX=ReDLat
DP=10
GQ=20

## I. Retain genotypes with Depth >= 10 & Quality >= 20
vcftools --gzvcf $PREFIX.X.vcf.gz --minDP $DP --minGQ $GQ --recode --recode-INFO-all --out $PREFIX.X.DP$DP.GQ$GQ
mv $PREFIX.X.DP$DP.GQ$GQ.recode.vcf $PREFIX.X.DP$DP.GQ$GQ.vcf
bgzip $PREFIX.X.DP$DP.GQ$GQ.vcf

## II. Check how this affected the missingness rates
# Per site
vcftools --gzvcf $PREFIX.X.DP$DP.GQ$GQ.vcf.gz --missing-site --out $PREFIX.X.DP$DP.GQ$GQ
## Per sample
vcftools --gzvcf $PREFIX.X.DP$DP.GQ$GQ.vcf.gz --missing-indv --out $PREFIX.X.DP$DP.GQ$GQ
```
Make some comparative plots
```{r compare-missingness-rates-X-GT}
#Declare the same GQ and DP values as in the filtering
DP_1 = 10
DP_2 = 20
GQ_X = 20

# We had previously generated an X.lmiss and X.imiss files of raw dataset
## I. Missingness per sample
imiss_GT_X1 = read.delim((paste0(PREFIX,".X.DP",DP_1,".GQ",GQ_X,".imiss")), header = T, sep = "")
imiss_GT_X2 = read.delim((paste0(PREFIX,".X.DP",DP_2,".GQ",GQ_X,".imiss")), header = T, sep = "")
                    
# Get stats
imiss_GT_X1_F_MISS = describe(imiss_GT_X1$F_MISS) 
rownames(imiss_GT_X1_F_MISS) = c(paste0("sample_chrX_missingness_F_GT.DP",DP_1))
imiss_GT_X2_F_MISS = describe(imiss_GT_X2$F_MISS) 
rownames(imiss_GT_X2_F_MISS) = c(paste0("sample_chrX_missingness_F_GT.DP",DP_2))

# Plot
boxplot(imiss$F_MISS, imiss_GT_X1$F_MISS,imiss_GT_X2$F_MISS,
        ylab = "Missingness rate",
        names = c("Raw dataset", paste0("Genotype QC DP=",DP_1), paste0("Genotype QC DP=",DP_2)),
        main = "Missingness rate per sample - Chromosome X genotype QC", 
        col = c("paleturquoise3", "lavender", "lightgoldenrod1"))
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/chrX-GT/missingness-rate_sample_boxplot_XGT-comparison.png)
```                    
lmiss_GT_X1 = read.delim((paste0(PREFIX,".X.DP",DP_1,".GQ",GQ_X,".lmiss")), header = T, sep = "")
lmiss_GT_X2 = read.delim((paste0(PREFIX,".X.DP",DP_2,".GQ",GQ_X,".lmiss")), header = T, sep = "")
# Get stats
lmiss_GT_X1_F_MISS = describe(lmiss_GT_X1$F_MISS) 
rownames(lmiss_GT_X1_F_MISS) = c(paste0("site_chrX_missingness_F_GT.DP",DP_1))
lmiss_GT_X2_F_MISS = describe(lmiss_GT_X2$F_MISS) 
rownames(lmiss_GT_X2_F_MISS) = c(paste0("site_chrX_missingness_F_GT.DP",DP_2))


boxplot(lmiss$F_MISS, lmiss_GT_X1$F_MISS, lmiss_GT_X2$F_MISS,
        ylab ="Missingness rate",
        names = c("Raw dataset", paste0("Genotype QC DP=",DP_1), paste0("Genotype QC DP=",DP_2)),
        main = "Missingness rate per site - Chromosome X genotype QC", 
        col = c("paleturquoise3", "lavender", "lightgoldenrod1"))
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/chrX-GT/missingness-rate_site_boxplot_XGT-comparison.png)
```

boxplot(lmiss$F_MISS, lmiss_GT_X1$F_MISS, lmiss_GT_X2$F_MISS,
        ylab ="Missingness rate",
        names = c("Raw dataset", paste0("Genotype QC DP=",DP_1), paste0("Genotype QC DP=",DP_2)),
        main = "Missingness rate per site - Chromosome X genotype QC <0.05", 
        col = c("paleturquoise3", "lavender", "lightgoldenrod1"),
        ylim = c(0,0.05))
```
![image](https://github.com/acostauribe/exome_qc_tutorial/blob/main/chrX-GT/missingness-rate_site_boxplot_XGT-comparison_lowend.png)
```
# Annotate the sample_metrics file with the new Missingness values
sample_metrics$missingness_Xchrom_GT_DP1 = imiss_GT_X1$F_MISS[match(sample_metrics$INDV, imiss_GT_X1$INDV)]
names(sample_metrics)[names(sample_metrics) == 'missingness_Xchrom_GT_DP1'] = paste0('missingness_Xchrom_GT_DP', DP_1)
sample_metrics$missingness_Xchrom_GT_DP2 = imiss_GT_X2$F_MISS[match(sample_metrics$INDV, imiss_GT_X2$INDV)]
names(sample_metrics)[names(sample_metrics) == 'missingness_Xchrom_GT_DP2'] = paste0('missingness_Xchrom_GT_DP', DP_2)
write.table(sample_metrics, "sample_metrics.txt")
```
