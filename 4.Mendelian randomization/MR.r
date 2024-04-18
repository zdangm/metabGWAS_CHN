args <- as.numeric(commandArgs(trailingOnly = TRUE))
#set up sub jobs
run_id<-1
run_list<-list()
for (outcome_i in c(30)){  #outcome id 
  for (array_i in 1:11){  #11 subjobs per job   only need about 20 min
    run_list[[run_id]]<-c(outcome_i, array_i)
    run_id=run_id+1
  }
}

run_i<-run_list[[args]]
print(run_i)

outcome_i = run_i[1]  #for outcome id
args = run_i[2]  #for array id

library(data.table)
library(foreach)
library(dplyr)
library(TwoSampleMR)
library(doParallel)
library(foreach)
library(ieugwasr)
library(mr.raps)
library(stringr)

####data####
exposure_folder <-"/path/metaGWAS/All_metab/P1e5_metab"
exposure_files <- list.files(exposure_folder) #, pattern = "*.xxx"

outcome_folder <- "/path/metaGWAS/gwas_summary_stat/bbj" #all
outcome_files <- list.files(outcome_folder, pattern = ".txt")
outcome_file <- outcome_files[outcome_i] 
outcome_name = sub('....$','',outcome_file)

#make dir for the outcome
output_dir = paste0('/path/metaGWAS/MR/result/',outcome_name,'/')
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

exposure_select <- c("SNP", "b", "se", "A1", "A2", "Freq", "p") #mlma
outcome_select <- c("SNP","b","se","A1","A2","freq","p")
exposure_colnames <- c("exp_SNP", "exp_beta", "exp_se", "exp_effect_allele", "exp_other_allele", "exp_eaf", "exp_pval")
outcome_colnames <- c("out_SNP", "out_beta", "out_se", "out_effect_allele", "out_other_allele", "out_eaf", "out_pval")

#### Parameters of MR####
p_select <- 1e-5
clump_kb <- 1000
clump_r2 <- 0.01

#check exist
exist_file = dir(output_dir)
output_files = paste0('result_file_array_',seq(1,length(exposure_files)),'.csv')
need_to_run = setdiff(output_files, exist_file)
need_to_run = as.numeric(sapply(need_to_run, function(x) sub('....$','',strsplit(x,"_")[[1]][4])))
exposure_files_original = exposure_files
exposure_files = exposure_files[need_to_run]


#id
i_start = (args-1)*600+1   #11 subjobs
i_end = min(c(args*600,length(exposure_files)))


#registerDoParallel(cores = detectCores())
registerDoParallel(cores = 90)

foreach(i = i_start:i_end) %dopar% {
  exposure_file <- file.path(exposure_folder, exposure_files[i])
  Metab <- fread(exposure_file)
  colnames(Metab) = c("SNP", "A1", "A2", "Freq", "b", "se",  "p", "N")
  exposure_raw_data <- as.data.frame(Metab)
  exposure_raw_data <- exposure_raw_data[, exposure_select]
  colnames(exposure_raw_data) <- exposure_colnames
  exposure_filter_pvalue <- exposure_raw_data[exposure_raw_data$exp_pval < p_select,]
  
  outcome_path <- file.path(outcome_folder, outcome_file)
  outcome <- fread(outcome_path)
  outcome_raw_data <- as.data.frame(outcome)
  outcome_raw_data <- outcome_raw_data[, outcome_select]
  colnames(outcome_raw_data) <- outcome_colnames
  
  #merge
  merge_exposure_and_outcome <- merge(exposure_filter_pvalue, outcome_raw_data, by.x = "exp_SNP", by.y = "out_SNP")
  df <- merge_exposure_and_outcome
  
  outcome_phenotype <- sub("\\.txt$", "", outcome_file)
  exposure_phenotype<- sub("\\.filtered.txt$", "",  exposure_files[i])
  print(exposure_phenotype)
  print(outcome_phenotype)
  nrow(df)
  
  # MR
  ans_NA = list()
  ans_NA$exposure = exposure_phenotype
  ans_NA$outcome = outcome_phenotype
  ans_NA$nsnp = NA
  ans_NA$mr_rap_beta = NA
  ans_NA$mr_rap_se = NA
  ans_NA$mr_rap_p=NA
  ans_NA$ivw_beta = NA
  ans_NA$ivw_se = NA
  ans_NA$ivw_p = NA
  ans_NA$egger_beta = NA
  ans_NA$egger_se = NA
  ans_NA$egger_p = NA
  ans_NA$wme_beta = NA
  ans_NA$wme_se = NA
  ans_NA$wme_p = NA
  ans_NA$wm_beta = NA
  ans_NA$wm_se = NA
  ans_NA$wm_p = NA
  ans_NA$sm_beta = NA
  ans_NA$sm_se = NA
  ans_NA$sm_p = NA
  ans_NA$egger_p_intercept = NA
  ans_NA$Q_pval = NA
  
  
  # 检查列是否为空
  if (nrow(df) >= 1) {
    exposure_df <- data.frame(
      SNP = df$exp_SNP,
      beta = df$exp_beta,
      se = df$exp_se,
      effect_allele = df$exp_effect_allele,
      other_allele=df$exp_other_allele,
      eaf=df$exp_eaf,
      pval=df$exp_pval,
      Phenotype=exposure_phenotype
    )
    
    outcome_df <- data.frame(
      SNP = df$exp_SNP,
      beta = df$out_beta,
      se = df$out_se,
      effect_allele = df$out_effect_allele,
      other_allele=df$out_other_allele,
      eaf=df$out_eaf,
      pval=df$out_pval,
      Phenotype=outcome_phenotype
    )
    
    
    exposure_dat <- format_data(exposure_df, type="exposure")
    outcome_dat <- format_data(outcome_df, type="outcome")
    print("the number of snps after mergeing exposue_df and outcome_df:")
    nrow(exposure_dat)
    nrow(outcome_dat)
    exposure_dat <- exposure_dat[complete.cases(exposure_dat[,c('effect_allele.exposure','other_allele.exposure')]),]
    outcome_dat <- outcome_dat[complete.cases(outcome_dat[,c('effect_allele.outcome','other_allele.outcome')]),]
    
    # harmonise
    dat <- harmonise_data(
      exposure_dat = exposure_dat,
      outcome_dat = outcome_dat,
      action=2
    )
    print("the number of snps after harmonise:")
    nrow(dat)
    length(which(dat$mr_keep=="TRUE"))
    dat <-dat[dat$mr_keep=="TRUE",]
    length_of_dat <- nrow(dat)
    print("okdata")
    nrow(dat)
    
    if (length_of_dat > 1) {
      # clump
      #if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
      #remotes::install_github("MRCIEU/genetics.binaRies")
      #local
      dat_clump <- dat %>% dplyr::select(rsid=SNP, pval=pval.exposure) %>%
        ieugwasr::ld_clump(., plink_bin="/data/tools/plink_v1.90b7.2/plink", #local plink
                           bfile="/path/metaGWAS/1kg/EAS",
                           pop='EAS',
                           clump_kb=clump_kb, clump_r2=clump_r2,)%>%{.$rsid}
      dat <- subset(dat, SNP %in% dat_clump)
      
      print("the number of snps after clumping:")
      nrow(dat)
      
      
      if (nrow(dat)>=1){
        #run MR
        ans_regular = mr(dat)
        ans_egger = mr_egger_regression(b_exp = dat$beta.exposure,
                                        se_exp = dat$se.exposure,
                                        b_out = dat$beta.outcome,
                                        se_out = dat$se.outcome)
        
        ans_rap=mr.raps(dat$beta.exposure, dat$beta.outcome,dat$se.exposure, dat$se.outcome)
        #summarize results
        ans = list()
        ans$exposure = as.character(exposure_phenotype)
        ans$outcome = as.character(outcome_phenotype)
        ans$nsnp = ans_regular[1,'nsnp']
        ans$mr_rap_beta = ans_rap$beta.hat
        ans$mr_rap_se = ans_rap$beta.se
        ans$mr_rap_p = ans_rap$beta.p.value
        
        
        row_index_ivw = which(ans_regular[,'method'] == 'Inverse variance weighted')
        ans$ivw_beta = ans_regular[row_index_ivw,'b']
        ans$ivw_se = ans_regular[row_index_ivw,'se']
        ans$ivw_p = ans_regular[row_index_ivw,'pval']
        
        row_index_egger = which(ans_regular[,'method'] == 'MR Egger')
        ans$egger_beta = ans_regular[row_index_egger,'b']
        ans$egger_se = ans_regular[row_index_egger,'se']
        ans$egger_p = ans_regular[row_index_egger,'pval']
        ans$egger_p_intercept = ans_egger$pval_i
        ans$Q_pval = ans_egger$Q_pval
        
        
        row_index_wme = which(ans_regular[,'method'] == 'Weighted median')
        ans$wme_beta = ans_regular[row_index_wme,'b']
        ans$wme_se = ans_regular[row_index_wme,'se']
        ans$wme_p = ans_regular[row_index_wme,'pval']
        
        row_index_sm = which(ans_regular[,'method'] == 'Simple mode')
        ans$sm_beta = ans_regular[row_index_sm,'b']
        ans$sm_se = ans_regular[row_index_sm,'se']
        ans$sm_p = ans_regular[row_index_sm,'pval']
        
        row_index_wm = which(ans_regular[,'method'] == 'Weighted mode')
        ans$wm_beta = ans_regular[row_index_wm,'b']
        ans$wm_se = ans_regular[row_index_wm,'se']
        ans$wm_p = ans_regular[row_index_wm,'pval']
        
      }else{
        ans = ans_NA
      }
      
    } else {
      ans = ans_NA
      # print("Skipping clump because the number of SNPs in dat is <= 1")
    }
    
  } else {
    ans = ans_NA
    # cat("Skipping file with less than 1 row of data:", exposure_file, "\n")
  }
  
  
  for(j in 1:length(ans)){
    if(length(ans[[j]]) == 0){
      ans[[j]] = NA
    }
  }
  
  ans = as.data.frame(ans)
  
  i_original = which(exposure_files_original == exposure_files[i])
  
  write.csv(ans, paste0(output_dir,'/result_file_array_',i_original,'.csv'), row.names = FALSE)
  
  
}
stopImplicitCluster()
