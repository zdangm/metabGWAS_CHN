####metabolite data ####
#ml load GCC/8.2.0 OpenMPI/3.1.4 pandas/0.24.2 R/3.6.0 
#ml GCC/10.2.0  OpenMPI/4.0.5 R/4.0.5
library(readxl)
library(data.table)
library(dplyr)
LHmetab<- fread("/path/data/metabolites/liuheng/metab/data/MWdata.txt", header=T,check.names=F,sep=" ")
metabolite_list_totest<-unique(colnames(LHmetab[,1:1912]))
  post_outlier_removal_mean_sd<-data.frame(
    metabolite=character(),
    mean=numeric(),
    sd=numeric(),
    standardized_mean=numeric(),
    standardized_sd=numeric())
  
  for (metabo in metabolite_list_totest) {
    dat_sub<-LHmetab %>% select(BLOODID,all_of(metabo))
    dat_sub$log_metabo<-apply(dat_sub, 1, function(x) log(as.numeric(x[2])))##natural log transform the data
    ##remove the measurements that out 3 SD 
    mean_metabo<-mean(dat_sub$log_metabo,na.rm = TRUE)
    sd_metabo<-sd(dat_sub$log_metabo,na.rm = TRUE)
    dat_sub_filter<-dat_sub %>% filter(abs((log_metabo-mean_metabo)/sd_metabo)<=3)
    ##standardize data
    mean_metabo_post<-mean(dat_sub_filter$log_metabo, na.rm = TRUE)
    sd_metabo_post<-sd(dat_sub_filter$log_metabo, na.rm = TRUE)
    dat_sub_filter$std_metabo<-apply(dat_sub_filter, 1, function(x) (as.numeric(x[3])-mean_metabo_post)/sd_metabo_post)
    phenoCol<-dat_sub_filter
    colnames(phenoCol)<-c("IID","metab","log_metab","Standardized_metabolites_level")
    standardized_mean<-mean(phenoCol$Standardized_metabolites_level, na.rm = T)
    standardized_sd<-sd(phenoCol$Standardized_metabolites_level, na.rm = T)
    tmp_dat<-data.frame(metabo, mean_metabo_post, sd_metabo_post, standardized_mean, standardized_sd)
    colnames(tmp_dat)<-c("metabolite", "mean","sd", "standardized_mean", "standardized_sd")
    post_outlier_removal_mean_sd<-rbind(post_outlier_removal_mean_sd,tmp_dat)
    sdphenoCol<-dat_sub_filter[,c(1,1,4)]
    write.table(sdphenoCol, paste0("/path/data/metabolites/liuheng/metab/metab_for_GWAS/sdmetab/sdpheno_",metabo,"_LH.txt"),sep="\t",quote = FALSE,row.names = FALSE, col.names = FALSE)

 }


 #### Covariate #### 

##grm
gcta64 --bfile /path/data/metabolites/liuheng/genome/liuheng_1159 --autosome --make-grm --out /path/data/metabolites/liuheng/genome/liuheng_1159_grm

##pc
gcta64 --grm  /path/data/metabolites/liuheng/genome/liuheng_1159_grm --pca 10 --thread-num 10 --out /path/data/metabolites/liuheng/genome/liuheng_1159_PC

##age
awk '{gsub(/[0-9]+_/, "", $2); print}' liuheng_1159_covar.txt > liuheng.txt
awk '{ $1 = $2; print }' liuheng.txt > liuheng2.txt


####GWAS analysis####
##demo
./gcta64 --mlma \
--bfile /path/data/metabolites/liuheng/genome/liuheng_1159 \
--grm /path/data/metabolites/liuheng/genome/grm_pc/liuheng_1159_grm \
--pheno /path/data/metabolites/liuheng/metab/metab_for_GWAS/sdmetab/sdpheno_1_LH.txt \
--covar /path/data/metabolites/liuheng/cov/liuheng_1159_covar.txt \
--out testmlm \
--thread-num 10

####sbatch.sh####
#!/bin/bash
#SBATCH --mail-user=  #############################
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint=haswell|skylake #sandybridge|
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=mlma
#SBATCH --array=1-1912

data_dir="/path/data/metabolites/liuheng/metab/metab_for_GWAS/sdmetab/"
files=($data_dir/*)
current_file=${files[$SLURM_ARRAY_TASK_ID-1]}
filename=$(basename -- "$current_file")
filename_no_ext="${filename%.*}"
# Run GCTA analysis 
cd /path/data/metabolites/liuheng/genome/gcta_1.94.0beta
./gcta64 --mlma \
--bfile /path/data/metabolites/liuheng/genome/liuheng_1159 \
--grm /path/data/metabolites/liuheng/genome/grm_pc/liuheng_1159_grm \
--pheno "$current_file" \
--qcovar /path/data/metabolites/liuheng/cov/liuheng_1159_qcovar.txt \
--covar /path/data/metabolites/liuheng/cov/liuheng_1159_covar.txt \
--out "/path/data/metabolites/liuheng/mlmaGWAS/GWAS/${filename_no_ext}"


