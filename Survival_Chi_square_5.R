##########################################################################################################
############################# Plotting lymphoma personalis data KM curves ################################
##########################################################################################################
############################### Following code adapted from komen2 project ###############################
###################### Griffith Lab git repo /git/komen2/Rscripts/SurvivalAnalysis.R #####################
##########################################################################################################
library(survival)
library(ggplot2)
library(multtest)
library(openxlsx)
library(plyr)

## Define output directory
outdir2 <-"/Users/fgomez/Felicia/lymphoma/ibrutinib_Personalis/survival_out"
setwd("/Users/fgomez/Felicia/lymphoma/ibrutinib_Personalis/survival_out")

##########################################################################################################
##### Read in files
##########################################################################################################
## Complete list of personalis variants -- Using this to get the full list of samples and gene names prior to filtering
personalis <- read.csv(file='~/Box Sync/LymphomaProject/P2C ibrutinib study data/Survival_analysis/test_input/all_personalis_variants.csv', header=T, stringsAsFactors=FALSE)

#pers_round2 <- read.table(file='~/Box\ Sync/LymphomaProject/P2C\ ibrutinib\ study\ data/data/personalis_filtered_recurrent.tsv', header=T, sep="\t", stringsAsFactors=FALSE)

## Readcount data#### 
merged_bams <- read.table(file='/Users/fgomez/Felicia/lymphoma/ibrutinib_Personalis/bam_readcount/all_samples_readcount_4.txt',header=T, fill=T, stringsAsFactors = FALSE)

### Manually reviwed readcount data-- NB these are sites that remained after filtering described below####
manual_review <- read.delim(file='~/Box\ Sync/LymphomaProject/P2C\ ibrutinib\ study\ data/data/bed_files/all_sample_manual_review_results.tsv',sep='\t', fill=T, header=T, stringsAsFactors=FALSE)

## ExAC data####
exac_data <- read.table(file='/Users/fgomez/Felicia/lymphoma/ibrutinib_Personalis/annotation/exac.pass_round_3.tsv',sep='\t',header=T, fill=T, stringsAsFactors = FALSE)

###MGI Annotations###
annotation_data <- read.table(file='/Users/fgomez/Felicia/lymphoma/ibrutinib_Personalis/annotation/personalis_unique_indel_fix_annotated_c.tsv',sep='\t',header=T, fill=T, stringsAsFactors = FALSE)

### List of variants with MGI annotations, ExAC allele frequencies and readcounts
#vars_mgi <- read.table(file='~/Box Sync/LymphomaProject/P2C ibrutinib study data/data/personalis_unique_indel_fix_annotated_c.tsv', sep="\t", header=T, stringsAsFactors=FALSE)

## List of variants that passed the TCGA breast cancer panel of normals
pers_pass <- read.table(file='~/Box Sync/LymphomaProject/P2C ibrutinib study data/data/pon_personalis_variants.pass.tsv', header=T, sep="\t", stringsAsFactors=FALSE)
## List of variants that failed the TCGA breast cancer panel of normals
pers_fail <- read.table(file='~/Box Sync/LymphomaProject/P2C ibrutinib study data/data/pon_personalis_variants.fail.tsv', header=T, sep="\t", stringsAsFactors=FALSE)

### Clinical data
clinicaldata <- read.xlsx("~/Box Sync/LymphomaProject/P2C ibrutinib study data/Mayo_Clinic/MC1282_timetoevent_021616.xlsx",colNames=T,rows=4:44)

##########################################################################################################
##### Input data cleanup
##########################################################################################################
##### Clean up personalis variant list
## Create new dataframe for manipulation
vars_pers <- personalis
## Add "PH" prefix to patient IDs in variant files
vars_pers$sample <- gsub("(\\S+)","PH\\1",vars_pers$sample)
## Harmonize column names
colnames(vars_pers)[colnames(vars_pers) == "Symbol"] <- "gene_name"
colnames(vars_pers)[colnames(vars_pers) == "Amino_Acid_Change"] <- "amino_acid_change"
## Create a copy of all personalis variants
vars_pers_all <- vars_pers

##### Merge MGI Annotations 
merged_bams_exac <- merge(merged_bams,exac_data,by=c("chromosome_name","start","stop","reference","variant"))
merged_bams_exac_annotation <- merge(merged_bams_exac,annotation_data,by=c("chromosome_name","start","stop","reference","variant"))
write.table(merged_bams_exac_annotation,"/Users/fgomez/Felicia/lymphoma/ibrutinib_Personalis/data/merged_bams_exac_annotation_2.tsv",sep="\t",col.names=T,quote=F, row.names=F)
write.table(merged_bams_exac_annotation,"/Users/fgomez/Box\ Sync/LymphomaProject/P2C\ ibrutinib\ study\ data/data/merged_bams_exac_annotation.tsv",sep="\t",col.names=T,quote=F, row.names=F)

vars_mgi <- merged_bams_exac_annotation
vars_mgi_2 <-merged_bams_exac_annotation

##### Clean up MGI file
## Add "PH" prefix to patient IDs in variant files
vars_mgi$sample <- gsub("(\\S+)","PH\\1",vars_mgi$sample)

#### Clean up clinical data column names
## Remove paranthetical information
colnames(clinicaldata) <- gsub("(\\S+)\\.\\((\\S+)", "\\1", colnames(clinicaldata))
## Rename sample column
colnames(clinicaldata) <- gsub("DCNTR_ID", "sample", colnames(clinicaldata))
## Rename overall survival columns and make things uniformly lower case
colnames(clinicaldata) <- tolower(gsub("FU|fu", "os", colnames(clinicaldata)))

#### Get the list of passed or failed variants from the panel of normals
## Add a new column to be merged in and filtered on
pers_pass$PON_results <- "PASS"
pers_fail$PON_results <- "FAIL" 
## Combine list of variants
pers_pon <- rbind(pers_pass,pers_fail)
## Add the panel of normals annotation
vars_mgi <- merge(vars_mgi, pers_pon)

##########################################################################################################
##### Filter sites
##########################################################################################################
setwd("/Users/fgomez/Box\ Sync/LymphomaProject/P2C\ ibrutinib\ study\ data/Survival_analysis/output_3")
outdir2 <- "/Users/fgomez/Box\ Sync/LymphomaProject/P2C\ ibrutinib\ study\ data/Survival_analysis/output_3"
do.call(file.remove, list(list.files("/Users/fgomez/Box\ Sync/LymphomaProject/P2C\ ibrutinib\ study\ data/Survival_analysis/output_2", full.names = TRUE)))

##### Create a stats table to follow filtering progression
stats <- matrix(c("Total Variants",nrow(vars_pers), 0, 0, nrow(vars_mgi), 0, 0), nrow=1, ncol=7)
colnames(stats) <- c("source","remaining_personalis","filtered_personalis","total_filtered_personalis","remaining_MGI","filtered_MGI","total_filtered_MGI")

##### Filter by coverage metrics
## Exclude variants with <100x coverage or >5000x coverage
vars_mgi <- vars_mgi[which(vars_mgi$coverage > 100 & vars_mgi$coverage < 5000),]
stats <- rbind(stats,c("Exclude <100x coverage or >5000x coverage", nrow(vars_pers), (as.numeric(stats[nrow(stats),2]) - nrow(vars_pers)), (as.numeric(stats[1,2]) - nrow(vars_pers)), nrow(vars_mgi), as.numeric(stats[nrow(stats),5]) - nrow(vars_mgi), (as.numeric(stats[1,5]) - nrow(vars_mgi))))

pdf(file="coverage_full_data.pdf")
hist(vars_mgi_2$coverage,xlim=c(0,30000), breaks=400, main= "Personalis Coverage - All data", xlab="total read count")
dev.off()

pdf(file="coverage_limited_data.pdf")
hist(vars_mgi_2$coverage,xlim = c(0,5000),breaks=400, main = "Personalis Coverage - Limited",xlab="total read count")
dev.off()

meancoverage <- mean(vars_mgi_2$coverage)

## Exclude variants with <5 reads of support
vars_mgi <- vars_mgi[which(vars_mgi$tumor_var_count >= 5),]
stats <- rbind(stats,c("Exclude <5 reads of support", nrow(vars_pers), (as.numeric(stats[nrow(stats),2]) - nrow(vars_pers)), (as.numeric(stats[1,2]) - nrow(vars_pers)), nrow(vars_mgi), as.numeric(stats[nrow(stats),5]) - nrow(vars_mgi), (as.numeric(stats[1,5]) - nrow(vars_mgi))))
## Exclude variants with VAF <2.5
vars_mgi <- vars_mgi[which(vars_mgi$tumor_VAF >= 2.5),]
stats <- rbind(stats,c("Exclude VAF <2.5", nrow(vars_pers), (as.numeric(stats[nrow(stats),2]) - nrow(vars_pers)), (as.numeric(stats[1,2]) - nrow(vars_pers)), nrow(vars_mgi), as.numeric(stats[nrow(stats),5]) - nrow(vars_mgi), (as.numeric(stats[1,5]) - nrow(vars_mgi))))

##### Filter by population frequencies
## Exclude variants that are not unannotated or <0.001 allele freq in 1000 genomes, ESP, UK10K and Korean Genomes
vars_pers <- vars_pers[which(vars_pers$X1000genomes_Total <= 0.001 | is.na(vars_pers$X1000genomes_Total)),] 
vars_pers <- vars_pers[which(vars_pers$ESP5400_Total <= 0.001 | is.na(vars_pers$ESP5400_Total)),] 
vars_pers <- vars_pers[which(vars_pers$UK10K_TWINS_Freq <= 0.001 | is.na(vars_pers$UK10K_TWINS_Freq)),] 
vars_pers <- vars_pers[which(vars_pers$KoreanGenomes_Freq <= 0.001 | is.na(vars_pers$KoreanGenomes_Freq)),] 
## Exclude variants with >0.001 ExAC allele frequency
vars_mgi <- vars_mgi[which(vars_mgi$ExAC_adj_AF <= 0.001 | is.na(vars_mgi$ExAC_adj_AF)),]
stats <- rbind(stats,c("Exclude >0.001 population frequency", nrow(vars_pers), (as.numeric(stats[nrow(stats),2]) - nrow(vars_pers)), (as.numeric(stats[1,2]) - nrow(vars_pers)), nrow(vars_mgi), as.numeric(stats[nrow(stats),5]) - nrow(vars_mgi), (as.numeric(stats[1,5]) - nrow(vars_mgi))))

##### Keep only nonsynonymous coding mutations for each annotation file
eff_type = c("NON_SYNONYMOUS_CODING","STOP_GAINED","CODON_INSERTION","EXON","FRAME_SHIFT","SPLICE_SITE_ACCEPTOR","SPLICE_SITE_DONOR","CODON_DELETION","CODON_CHANGE_PLUS_CODON_DELETION","START_GAINED","START_LOST","NON_SYNONYMOUS_START","STOP_LOST","CODON_CHANGE_PLUS_CODON_INSERTION")
vars_pers <- vars_pers[vars_pers$Effect %in% eff_type, ]
mut_type = c("frame_shift_del","frame_shift_ins","in_frame_del","missense","nonsense","nonstop","splice_site","splice_site_del","splice_site_ins","in_frame_ins","rna")
vars_mgi <- vars_mgi[vars_mgi$trv_type %in% mut_type, ]
stats <- rbind(stats,c("Restrict to exonic and splice site", nrow(vars_pers), (as.numeric(stats[nrow(stats),2]) - nrow(vars_pers)), (as.numeric(stats[1,2]) - nrow(vars_pers)), nrow(vars_mgi), as.numeric(stats[nrow(stats),5]) - nrow(vars_mgi), (as.numeric(stats[1,5]) - nrow(vars_mgi))))

##### Exclude variants that failed the panel of normals
vars_mgi <- vars_mgi[which(vars_mgi$PON_results == "PASS"),]
stats <- rbind(stats,c("Exclude recurrent in TCGA BRCA panel of normals", nrow(vars_pers), (as.numeric(stats[nrow(stats),2]) - nrow(vars_pers)), (as.numeric(stats[1,2]) - nrow(vars_pers)), nrow(vars_mgi), as.numeric(stats[nrow(stats),5]) - nrow(vars_mgi), (as.numeric(stats[1,5]) - nrow(vars_mgi))))

##### Exclude relapse samples
relapse <- c("1966RR","1990RR","2060RR","PH1966RR","PH1990RR","PH2060RR")
vars_mgi <- vars_mgi[!(vars_mgi$sample %in% relapse),]
vars_pers_all <- vars_pers_all[!(vars_pers_all$sample %in% relapse),]
vars_pers <- vars_pers[!(vars_pers$sample %in% relapse),]
## Add figures to stats table
stats <- rbind(stats,c("Excluding relapse samples", nrow(vars_pers), (as.numeric(stats[nrow(stats),2]) - nrow(vars_pers)), (as.numeric(stats[1,2]) - nrow(vars_pers)), nrow(vars_mgi), as.numeric(stats[nrow(stats),5]) - nrow(vars_mgi), (as.numeric(stats[1,5]) - nrow(vars_mgi))))

##### Add a statistic for variants in genes mutated in 3 or more samples
rec_pers <- ddply(vars_pers, "gene_name", nrow)
rec_mgi <- ddply(vars_mgi, "gene_name", nrow)
stats <- rbind(stats,c("Restrict to genes mutated in >2 samples", nrow(vars_pers[vars_pers$gene_name %in% rec_pers[rec_pers$V1 >2,"gene_name"],]), (as.numeric(stats[nrow(stats),2]) - nrow(vars_pers[vars_pers$gene_name %in% rec_pers[rec_pers$V1 >2,"gene_name"],])), (as.numeric(stats[1,2]) - nrow(vars_pers[vars_pers$gene_name %in% rec_pers[rec_pers$V1 >2,"gene_name"],])), nrow(vars_mgi[vars_mgi$gene_name %in% rec_mgi[rec_mgi$V1 >2,"gene_name"],]), as.numeric(stats[nrow(stats),5]) - nrow(vars_mgi[vars_mgi$gene_name %in% rec_mgi[rec_mgi$V1 >2,"gene_name"],]), (as.numeric(stats[1,5]) - nrow(vars_mgi[vars_mgi$gene_name %in% rec_mgi[rec_mgi$V1 >2,"gene_name"],]))))

##Number of Recurrent genes##
number_of_genes_pers <-nrow(rec_pers[which(rec_pers$V1 >2),])
number_of_genes_mgi <-nrow(rec_mgi[which(rec_mgi$V1 >2),])

vars_mgi<- vars_mgi[vars_mgi$gene_name %in% rec_mgi[rec_mgi$V1 >2,"gene_name"],]
vars_pers <-vars_pers[vars_pers$gene_name %in% rec_pers[rec_pers$V1 >2,"gene_name"],]

###Filter manual review 
manual_review$call_upper <-toupper(manual_review$call)
manual_review$call <-NULL 
manual_review$call <-manual_review$call_upper
manual_review$call_upper <-NULL 
somatic_manual_review <- manual_review[which(manual_review$call=="S"),]
#test_merge <-merge(x,manual_review,by=c("chromosome_name","start","stop", "reference","variant","sample"))
#vars_pers_test <-merge(y,somatic_manual_review,by=c("chromosome_name","start","stop", "reference","variant","sample"), all.y=TRUE)
vars_mgi <- merge(vars_mgi,somatic_manual_review,by=c("chromosome_name","start","stop", "reference","variant","sample"))
stats <- rbind(stats,c("Exclude sites that failed manual review", nrow(vars_pers), (as.numeric(stats[nrow(stats),2]) - nrow(vars_pers)), (as.numeric(stats[1,2]) - nrow(vars_pers)), nrow(vars_mgi), as.numeric(stats[nrow(stats),5]) - nrow(vars_mgi), (as.numeric(stats[1,5]) - nrow(vars_mgi))))

##### Add a statistic for variants in genes mutated in 3 or more samples -AGAIN After manual review filter#####
rec_pers_b <- ddply(vars_pers, "gene_name", nrow)
rec_mgi_b <- ddply(vars_mgi, "gene_name", nrow)
stats <- rbind(stats,c("Restrict to genes mutated in >2 samples following manual review", nrow(vars_pers[vars_pers$gene_name %in% rec_pers_b[rec_pers_b$V1 >2,"gene_name"],]), (as.numeric(stats[nrow(stats),2]) - nrow(vars_pers[vars_pers$gene_name %in% rec_pers_b[rec_pers_b$V1 >2,"gene_name"],])), (as.numeric(stats[1,2]) - nrow(vars_pers[vars_pers$gene_name %in% rec_pers_b[rec_pers_b$V1 >2,"gene_name"],])), nrow(vars_mgi[vars_mgi$gene_name %in% rec_mgi_b[rec_mgi_b$V1 >2,"gene_name"],]), as.numeric(stats[nrow(stats),5]) - nrow(vars_mgi[vars_mgi$gene_name %in% rec_mgi_b[rec_mgi_b$V1 >2,"gene_name"],]), (as.numeric(stats[1,5]) - nrow(vars_mgi[vars_mgi$gene_name %in% rec_mgi_b[rec_mgi_b$V1 >2,"gene_name"],]))))
##Number of Recurrent genes##
number_of_genes_pers_b <-nrow(rec_pers_b[which(rec_pers_b$V1 >2),])
number_of_genes_mgi_b <-nrow(rec_mgi_b[which(rec_mgi_b$V1 >2),])
stats <-rbind(stats,c("Number of recurrent genes",nrow(rec_pers_b[which(rec_pers_b$V1 >2),]), "-","-",nrow(rec_mgi_b[which(rec_mgi_b$V1 >2),]),"-","-"))

##### Write several summary files
#rec_pers_2 <-as.data.frame(rec_pers_b[which(rec_pers_b$V1 >2),])
rec_pers_df <-as.data.frame(rec_pers_b)
colnames(rec_pers_df) <-c("gene_name","number_of_mutants")

rec_mgi_df <-as.data.frame(rec_mgi_b)
colnames(rec_mgi_df) <-c("gene_name","number_of_mutants")
vars_mgi_out <-as.data.frame(vars_mgi)

write.table(stats, 'filtering_stats.tsv', quote=F, row.names=F, sep="\t")
write.table(rec_pers_df, 'personalis_recurrence_data.tsv',row.names=F,  quote=F,sep="\t")
write.table(rec_mgi_df, 'mgi_recurrence_data.tsv', row.names=F, quote=F,  sep="\t")
write.table(vars_mgi_out, 'mgi_annotation_recurrence_filtered_manual_review_filtered.tsv', col.names= T,row.names=F, quote=F,  sep="\t")

##### Restrict to only necessary columns
## Restrict the MGI annotated variant list to only necessary columns 
vars_mgi <- vars_mgi[,c("sample","gene_name","trv_type","c_position","amino_acid_change")]
## Restrict the personalis variant list to only necessary columns
vars_pers <- vars_pers[,c("sample","gene_name","amino_acid_change","Effect")]
## Restrict the personalis variant list to only necessary columns
vars_pers_all <- vars_pers_all[,c("sample","gene_name","amino_acid_change","Effect")]

##########################################################################################################
##### Recode clinical data event columns
##########################################################################################################
## Recode PFS event status from 1=censored, 2=event ==> 0=censored, 1=event
clinicaldata$pfs_stat_recoded <- clinicaldata$pfs_stat
clinicaldata$pfs_stat_recoded <- gsub("1","0", clinicaldata$pfs_stat_recoded)
clinicaldata$pfs_stat_recoded <- gsub("2","1", clinicaldata$pfs_stat_recoded)

## Recode overall event status from Alive, Dead ==> 0=Aive, 1=Dead
clinicaldata$os_stat_recoded <- clinicaldata$os_stat
clinicaldata$os_stat_recoded <- gsub("Alive","0", clinicaldata$os_stat_recoded)
clinicaldata$os_stat_recoded <- gsub("Dead","1", clinicaldata$os_stat_recoded)

## Recode duration of response event status from 1=censored, 2=event  ==> 0=censored, 1=event
clinicaldata$dur_stat_recoded <- clinicaldata$dur_stat
clinicaldata$dur_stat_recoded <- gsub("1","0", clinicaldata$dur_stat_recoded)
clinicaldata$dur_stat_recoded <- gsub("2","1", clinicaldata$dur_stat_recoded)
##########################################################################################################
##### Create mutation matricies
##########################################################################################################
## Create a list genes from the original list of personalis data, exclude blank gene symbol   
#genes <- unique(vars_mgi[vars_mgi$gene_name != "","gene_name"])

## Create a list of genes from the receurrent genes in the MGI list 
vars_mgi <-vars_mgi[vars_mgi$gene_name %in% rec_mgi_b[rec_mgi_b$V1 >2,"gene_name"],]
#genes <- unique(vars_mgi[vars_mgi$gene_name != "","gene_name"])
genes <- unique(vars_mgi$gene_name)

define_sample <-function(clin_dat,outcome){
  samples <- na.omit(clin_dat[, c("sample", outcome)])$sample
  clin_dat <-clin_dat[clin_dat$sample %in% samples,]
  return(clin_dat)
  }
clinicaldata_s <- define_sample(clinicaldata,"pfs_stat_recoded")
#clinicaldata_s <- define_sample(clinicaldata,"dur_stat_recoded")

## Create a list of samples
#samples <- as.factor(unique(clinicaldata[which(!(is.na(clinicaldata$dur_stat_recoded))),"sample"]))
samples <- as.factor(unique(clinicaldata_s[,"sample"]))
## Determine minimum number of mutations to plot
min_mut <- 3
## Set a raw p-value cutoff for shortened output file
rawpcutoff <- 0.05

## Create a matrix of mutation status by gene (All personalis data)
mutations_pers_all=matrix(data=0,nrow=length(samples),ncol=length(genes), dimnames=list(samples,genes))
for (i in 1:length(genes)){
  gene=genes[i]
  gene_data=vars_pers_all[which(vars_pers_all[,"gene_name"]==gene),]
  mut_samples=unique(gene_data[,"sample"])
  mut_status=data.frame(samples,"0")
  mut_status=apply(mut_status,2,as.vector)
  rownames(mut_status)=samples
  colnames(mut_status)[2]="mut_status"
  mut_status[which(mut_status[,"samples"]%in%mut_samples),"mut_status"]=1
  mutations_pers_all[samples,gene]=mut_status[samples,"mut_status"]
}

## Create a matrix of mutation status by gene (filtered personalis data)
mutations_pers=matrix(data=0,nrow=length(samples),ncol=length(genes), dimnames=list(samples,genes))
for (i in 1:length(genes)){
  gene=genes[i]
  gene_data=vars_pers[which(vars_pers[,"gene_name"]==gene),]
  mut_samples=unique(gene_data[,"sample"])
  mut_status=data.frame(samples,"0")
  mut_status=apply(mut_status,2,as.vector)
  rownames(mut_status)=samples
  colnames(mut_status)[2]="mut_status"
  mut_status[which(mut_status[,"samples"]%in%mut_samples),"mut_status"]=1
  mutations_pers[samples,gene]=mut_status[samples,"mut_status"]
}

# Create a matrix of mutation status by gene (mgi annotated data)
mutations_mgi=matrix(data=0,nrow=length(samples),ncol=length(genes), dimnames=list(samples,genes))
for (i in 1:length(genes)){
  gene=genes[i]
  gene_data=vars_mgi[which(vars_mgi[,"gene_name"]==gene),]
  mut_samples=unique(gene_data[,"sample"])
  mut_status=data.frame(samples,"0")
  mut_status=apply(mut_status,2,as.vector)
  rownames(mut_status)=samples
  colnames(mut_status)[2]="mut_status"
  mut_status[which(mut_status[,"samples"]%in%mut_samples),"mut_status"]=1
  mutations_mgi[samples,gene]=mut_status[samples,"mut_status"]
}


#########################################################################################################
##### Set up parameters for survival plotting function calls
##########################################################################################################
## Name the test comparisons
test_types=c("PFS_vs_filtered_MGI_list")
## Identify columns with event status
eventvars=c("pfs_stat_recoded")
## Identify columns with time to event/censor for patient
timevars=c("pfs_mos")
## Identify the value that corresponds to an event
eventvals=c("1")
## List the input maticies for each comparison
mutationdatalist=list(mutations_mgi)

## Name the test comparisons
#test_types=c("DurRes_vs_filtered_MGI_list")
## Identify columns with event status
#eventvars=c("dur_stat_recoded")
## Identify columns with time to event/censor for patient
#timevars=c("dur_mos")
## Identify the value that corresponds to an event
#eventvals=c("1")
## List the input maticies for each comparison
#mutationdatalist=list(mutations_mgi)


## Name the test comparisons
#test_types=c("PFS_vs_full_personalis_list","PFS_vs_filtered_personalis_list","PFS_vs_filtered_MGI_list")
## Identify columns with event status
#eventvars=c("pfs_stat_recoded","pfs_stat_recoded","pfs_stat_recoded")
## Identify columns with time to event/censor for patient
#timevars=c("pfs_mos","pfs_mos","pfs_mos")
## Identify the value that corresponds to an event
#eventvals=c("1","1","1")
## List the input maticies for each comparison
#mutationdatalist=list(mutations_pers_all,mutations_pers,mutations_mgi)

# ## Name the test comparisons
# test_types=c("PFS_vs_filtered_personalis_list","PFS_vs_filtered_MGI_list")
# ## Identify columns with event status
# eventvars=c("pfs_stat_recoded","pfs_stat_recoded")
# ## Identify columns with time to event/censor for patient
# timevars=c("pfs_mos","pfs_mos")
# ## Identify the value that corresponds to an event
# eventvals=c("1","1")
# ## List the input maticies for each comparison
# mutationdatalist=list(mutations_pers,mutations_mgi)


#Define function to perform survival analysis and plot KM curve
#gene="FBXO11"
#samples=samples
#mutations=mutations_mgi
#min_mutations=min_mut
#clinicaldata=clinicaldata
#timevar="pfs_mos"
#eventvar="pfs_stat_recoded"
#eventval="1"
#survtype="PFS_vs_filtered_MGI_list"
#resultdir=paste(outdir2,"PFS_vs_filtered_MGI_list_test",sep="/")

plotKM=function(gene,samples,mutations,min_mutations,clinicaldata,timevar,eventvar,eventval,survtype,resultdir){
  #Create vector to store results and initiate as NA
  results=rep("NA",5)
  names(results)=c("n_wt","n_mt","n_events","p","HR")
  
  print(paste("processing ",gene, sep=""))
  mut_status=mutations[samples,gene]
  
  #skip genes with too few mutations
  if (sum(as.numeric(mut_status))<min_mutations){
    return(results)
  }
  
  #skip genes with all samples mutated
  if (sum(as.numeric(mut_status))==length(samples)){
    return(results)
  }
  
  #Create dataframe for survival analysis
  mutations_mgi_2 <-as.data.frame(mutations_mgi)
  mutations_mgi_2$sample <-rownames(mutations_mgi_2)
  #surv_data=data.frame(clinicaldata[clinicaldata$sample %in% samples,timevar],clinicaldata[clinicaldata$sample %in% samples,eventvar],mut_status)
  surv_data <- merge(mutations_mgi_2[,c("sample",gene)], clinicaldata[,c("sample",timevar,eventvar)])
  colnames(surv_data)=c("sample","mut_status","time","event")
  
  #Count number of events
  n_events=length(which(clinicaldata[clinicaldata$sample %in% samples,eventvar]==eventval))
  
  #Exclude patients with survival event missing
  surv_data=surv_data[which(!is.na(surv_data[,"event"])),]
  
  #Create a survival object
  surv_data.surv = with(surv_data, Surv(time, event==eventval))
  
  #Calculate p-value
  survdifftest=survdiff(surv_data.surv ~ mut_status, data = surv_data)
  survpvalue = 1 - pchisq(survdifftest$chisq, length(survdifftest$n) - 1)
  survpvalue = format(as.numeric(survpvalue), digits=3)
  krfit.by_RFgroup = survfit(surv_data.surv ~ mut_status, data = surv_data)
  
  #Calculate a hazard ratio
  #http://stats.stackexchange.com/questions/4528/what-is-the-difference-between-the-coef-and-expcoef-output-of-coxph-in-r
  fit1 = coxph(Surv(time, event==1) ~ mut_status, data = surv_data)
  hazard_ratio=format(exp(fit1$coefficients), digits=3)
  
  #Grab Confidence Intervals from coxph
  fit_sum <- summary(fit1, conf.int=.95)
  hr_int <- fit_sum$conf.int[,3:4]
  hr_low <- format(hr_int[1], digits=3)
  hr_high <- format(hr_int[2], digits=3)
  
  #Plot survival curve
  setwd(resultdir)
  pdf(file=paste(gene,"_",survtype,".pdf",sep=""))
  par(mar=c(5,5,4,2) + 0.1)
  #colors = rainbow(5)
  colors = c("red","blue")
  title=paste(survtype," ",gene," mut status", sep="")
  plot(krfit.by_RFgroup, col = colors, xlab = "Time (Months)", ylab = "Survival", main=title, cex.axis=1.5, cex.lab=1.6)
  groups=sort(unique(surv_data[,"mut_status"])) #returns unique factor levels sorted alphabetically
  names(colors)=groups
  groups_custom=c("0","1")
  colors_custom=colors[groups_custom]
  group_sizes_custom=table(surv_data[,"mut_status"])[groups_custom]
  groups_custom=c("WT","MT") #Reset names for readability
  legend_text=c(paste(groups_custom, " ", "(", group_sizes_custom, ")", sep=""),paste("p =", survpvalue, sep=" "),paste("HR =", hazard_ratio, "(95% CI: ",hr_low," - ",hr_high,")", sep=" "))
  legend(x = "topright", legend = legend_text, col = c(colors_custom,"white","white"), lty = "solid", bty="n", cex=1.5)
  dev.off()
  
  results[c("n_wt","n_mt","n_events","p","HR","CI_low","CI_high")]=c(group_sizes_custom,n_events,survpvalue,hazard_ratio,hr_low,hr_high)
  return(results)
}

#Create folders to store each set of KM plots
sapply(test_types,dir.create,showWarnings = FALSE)

#Create matrix to store survival results
tests <- as.vector(t(outer(test_types, genes, paste, sep="_")))
surv_results=matrix(data=NA,nrow=length(tests),ncol=10,dimnames=list(tests,c("test","test_type","gene","n_wt","n_mt","n_events","p","HR","CI_low","CI_high")))

for (i in 1:length(test_types)){
  for (j in 1:length(genes)){
    test_type=test_types[i]
    gene=genes[j]
    test=paste(test_type,gene,sep="_")
    resultdir=paste(outdir2,test_type,sep="/")
    result=plotKM(gene=gene,samples=samples,mutations=mutationdatalist[[i]],min_mutations=min_mut,clinicaldata=clinicaldata_s,
                  timevar=timevars[i],eventvar=eventvars[i],eventval=eventvals[i],survtype=test_type,resultdir=resultdir)
    
    #Save results in dataframe
    
    surv_results[test,c("test","test_type","gene","n_wt","n_mt","n_events","p","HR","CI_low","CI_high")]=c(test,test_type,gene,result[c("n_wt","n_mt","n_events","p","HR","CI_low","CI_high")])
  }
}

#Perform multiple testing corrections
#Correct p-values
pvalues=as.numeric(surv_results[,"p"])
pvalues_adj=mt.rawp2adjp(pvalues, proc=c("BH"), na.rm=T)
pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
surv_results=cbind(surv_results, pvalues_adj_orig_order) 


#Save survival results
setwd(outdir2)
test_split <-strsplit(test_types,"_")
test <-test_split[[1]][1]

write.table(surv_results, file=paste0("SurvivalPvaluesTable_",test,".tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
#Output results with raw p-value < defined above
write.table(surv_results[which(as.numeric(surv_results[,"rawp"])<rawpcutoff),], file=paste0("SurvivalPvaluesTable_pvaluecutoff_",test,".tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

####################################################################################################################################
###### Do Chi square analysis
###################################################################################################################################

mutations_mgi_3 <-(mutations_mgi)
#mutations_mgi_3$sample <-rownames(mutations_mgi)

#clinicaldata_2 <- merge(clinicaldata,mutations_mgi_3,by.x ="sample",by.y = "row.names")
clinicaldata_2 <-clinicaldata[,1:12]

marker_names <- colnames(mutations_mgi_3)
marker_data <-mutations_mgi_3

prognostic_names=c("bestresp")
prognostic_data=clinicaldata_2[,c("bestresp","sample")]

prognostic_data <- na.omit(prognostic_data)
marker_data <-marker_data[prognostic_data$sample,]


rec_marker_data <- colSums(as.numeric(as.array(marker_data)))

marker_vs_prog_var=vector(length=length(marker_names)*length(prognostic_names))
n=0
for (a in 1:length(marker_names)){
  for (b in 1:length(prognostic_names)){
    n=n+1
    marker_vs_prog_var[n]=paste(marker_names[a],"_x_",prognostic_names[b],sep="")
  }
}

#Create array to store results for prognostic contingency table stats
table_results = array(0, dimnames = list(marker_vs_prog_var, c("pvalue", "test","mutations")), dim=c(length(marker_vs_prog_var),3))

#Loop through all markers for categorical tests
n=0
for (a in 1:length(marker_names)){
  #loop through all categorical prognostic variables and perform contingency table statistics
  for (b in 1:length(prognostic_names)){
    n=n+1
    
    #Skip any marker/prognostic variable combinations which do not have at least two factors when cases are excluded (i.e. due to missing values)
    
    if (nlevels(as.factor(as.vector(marker_data[,a][!is.na(prognostic_data[,b])])))<2){
      table_results[n,"pvalue"]=NA
      table_results[n,"test"]=NA
      next
    }
    if (nlevels(as.factor(as.vector(prognostic_data[,b][!is.na(marker_data[,a])])))<2){
      table_results[n,"pvalue"]=NA
      table_results[n,"test"]=NA
      next
    }
    #Calculate Pearson Chi-square for all markers vs diagnosis
    chisquare_result=chisq.test(x=marker_data[,a], y=prognostic_data[,b])
    observed=chisquare_result$observed
    expected=chisquare_result$expected
    #If one or more cells in the table have expected count less than or equal to 5 use a Fisher test instead of Pearson Chi-square
    if (length(expected[expected<=5])>0){
      fisher_test_result=fisher.test(x=marker_data[,a], y=prognostic_data[,b], alternative = "two.sided")
      table_results[n,"pvalue"]=fisher_test_result$p.value
      table_results[n,"test"]="Fisher"
    }else{
      #Enter results into array for print out
      table_results[n,"pvalue"]=chisquare_result$p.value
      table_results[n,"test"]="Pearson_Chi"
    }
  }
  table_results[n,"mutations"]=sum(as.numeric(marker_data[,a]))
}

#Remove rows where no p-value could be calculated.  This screws up multtest even though it shouldn't matter
table_results=table_results[!is.na(table_results[,1]),]

#Correct p-values
pvalues=as.numeric(table_results[,"pvalue"])
pvalues_adj=mt.rawp2adjp(pvalues, proc=c("Bonferroni","BH"))
pvalues_adj_orig_order=pvalues_adj$adjp[order(pvalues_adj$index),]
table_results=cbind(table_results, pvalues_adj_orig_order[,2:3])
table_results_2 <- as.data.frame(table_results)
table_results_2$test <-rownames(table_results_2)
rownames(table_results_2) <-c()
table_results_2 <- table_results_2[c(2,3,1,4,5)]


#Create folder to stoer Chi Square Results
dir.create("/Users/fgomez/Box\ Sync/LymphomaProject/P2C\ ibrutinib\ study\ data/Survival_analysis/output_3/ChiSquare_results",showWarnings = FALSE)
outdir3 <- "/Users/fgomez/Box\ Sync/LymphomaProject/P2C\ ibrutinib\ study\ data/Survival_analysis/output_3/ChiSquare_results"
setwd(outdir3)
write.table(table_results_2, file="ChiSquare_pvals.tsv", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE,append=FALSE)

