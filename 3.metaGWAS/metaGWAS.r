
# =====================get txt for metal =====================
library(data.table)
index_table <- fread("path/meta/meta_index/match.txt", header = TRUE, sep = "\t")

for (i in 1:nrow(index_table)) {
  index <- index_table$Index[i]
  ZMSCID <- index_table$ID_ZMSC[i]
  THSBCID <- index_table$ID_THSBC[i]
  WeID <- index_table$We_ID[i]

  txt_content <- paste(
    "AVERAGEFREQ ON",
    "MINMAXFREQ ON",
    "SCHEME   STDERR",
    "STDERR   se",
    "MARKER   SNP",
    "ALLELE   A1 A2",
    "FREQ     Freq",
    "EFFECT   b",
    "PVAL     p",
    "",
    paste("PROCESS /path/meta/ZMSC/GWAS/summary/sdpheno_", ZMSCID, "_ZMSC.mlma", sep = ""),

    "AVERAGEFREQ ON",
    "MINMAXFREQ ON",
    "SCHEME   STDERR",
    "STDERR   se",
    "MARKER   SNP",
    "ALLELE   A1 A2",
    "FREQ     Freq",
    "EFFECT   b",
    "PVAL     p",
    "",
    paste("PROCESS /path/meta/THSBC/GWAS/summary/sdpheno_", THSBCID, "_TSBC.mlma", sep = ""),
    
    
    "AVERAGEFREQ ON",
    "MINMAXFREQ ON",
    "SCHEME   STDERR",
    "STDERR   se",
    "MARKER   SNP",
    "ALLELE   A1 A2",
    "FREQ     Freq",
    "EFFECT   b",
    "PVAL     p",
    "",
    paste("PROCESS /path/meta/webirth/mlma/", WeID, ".mlma", sep = ""),
    paste("OUTFILE /path/meta/meta/result/",index," .tbl",sep = ""),
    "ANALYZE HETEROGENEITY",
    sep = "\n"
  )

  txt_path <- paste0("/path/meta/meta_txt/",index,".txt")
  writeLines(txt_content, txt_path)
}



#  =====================metal.sh =====================
#!/bin/bash
#SBATCH --mail-user=  #############################
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --constraint=haswell|skylake #sandybridge|
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --job-name=meta
#SBATCH --array=1-4

metal_path="/path/tools/generic-metal/metal"
txt_folder="/path/meta/meta_txt/"

all_txt_files=("$txt_folder"/*.txt)

total_files=${#all_txt_files[@]}
files_per_task=$((total_files / SLURM_ARRAY_TASK_COUNT))


start_index=$(( (SLURM_ARRAY_TASK_ID - 1) * files_per_task ))
end_index=$(( SLURM_ARRAY_TASK_ID * files_per_task - 1))

echo "Total files: $total_files"
echo "Files per task: $files_per_task"

echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Start index: $start_index"
echo "End index: $end_index"


for ((i=start_index; i<=end_index; i++)); do
  txt_file="${all_txt_files[$i]}"
  if [ -f "$txt_file" ]; then
    metal_command="$metal_path $txt_file"
    $metal_command
    echo "Processed: $txt_file"
  else
    echo "File not found: $txt_file"
  fi
done

echo "All is well :)"




# =====================Extraction of significant SNP (1e-5 SNP) =====================
input_dir="/path/meta/meta/result/tbl"
output_dir="/path/meta/meta/result/P1e5_result"

for tbl_file in "$input_dir"/*.tbl; do
    base_name=$(basename -- "$tbl_file")
    file_name="${base_name%.tbl}"

    awk -F'\t' 'NR==1 || ($10 < 1e-5 && ($11 == "+++" || $11 == "---")) {print}' "$tbl_file" > "$output_dir/$file_name.filtered.txt"
done





