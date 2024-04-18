#### GCTA Heritability ####
#demo
cd /path/metabolites/liuheng/genome/gcta_1.94.0beta
./gcta64 \
--grm /path/metabolites/liuheng/genome/grm_pc/liuheng_1159_grm \
--reml \
--pheno /path/metabolites/liuheng/metab/metab_for_GWAS/sdmetab/sdpheno_1_ZMSC.txt \
--qcovar /path/metabolites/liuheng/cov/liuheng_1159_qcovar.txt \
--covar /path/metabolites/liuheng/cov/liuheng_1159_covar.txt \
--thread-num 10 \
--out test

# run.sh
#!/bin/bash
#SBATCH --mail-user=  #############################
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint=haswell|skylake #sandybridge|
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=h2
##SBATCH --array=1-19
#SBATCH --array=1

cd /path/metabolites/liuheng/genome/gcta_1.94.0beta
data_folder="/path/metabolites/liuheng/metab/metab_for_GWAS/sdmetab"
base_command="./gcta64 \
--grm /path/metabolites/liuheng/genome/liuheng_1159_grm \
--reml \
--covar /path/metabolites/liuheng/cov/liuheng_1159_covar.txt \
--thread-num 10 \
--qcovar /path/metabolites/liuheng/cov/liuheng_1159_qcovar.txt \
"


file_ranges=("1-100" "101-200" "201-300" "301-400" "401-500" "501-600" "601-700" "701-800" "801-900" "901-1000" "1001-1100" "1101-1200" "1201-1300" "1301-1400" "1401-1500" "1501-1600" "1601-1700" "1701-1800" "1801-1912")

current_range=${file_ranges[${SLURM_ARRAY_TASK_ID}-1]}

start=$(echo $current_range | cut -d'-' -f1)
end=$(echo $current_range | cut -d'-' -f2)

for ((i=start; i<=end; i++)); do

  phenotype_file="$data_folder/sdpheno_${i}_ZMSC.txt"

  output_file="/path/metabolites/liuheng/GWAS_result/heritability/h2/sdpheno_${i}_gwa"

  command="$base_command --pheno $phenotype_file --out $output_file"
  echo "Running command: $command"
  $command
  echo "Phenotype $i processed."
done


#combine files

ml GCC/10.2.0  OpenMPI/4.0.5 R/4.0.5
library(data.table)
file_list <- list.files("/path/metabolites/liuheng/GWAS_result/heritability/h2/hsq", full.names = TRUE)
summary_data <- data.table(filename = character(), h2 = numeric(), h2_SE = numeric(), P = numeric(), n = numeric())
for (file_path in file_list) {
  a <- fread(file_path, fill = TRUE,header=F,sep="\t")
  
 h2 <- as.numeric(a[5, 2])
  h2_SE <- as.numeric(a[5, 3])
  P <- as.numeric(a[10, 2])
  n <- as.numeric(a[11, 2])
  filename <- basename(file_path)
  summary_data <- rbind(summary_data, data.table(filename = filename, h2 = h2, h2_SE = h2_SE, P = P, n = n))
}
fwrite(summary_data, "ZMSC_h.csv")
