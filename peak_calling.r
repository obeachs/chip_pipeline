library(ChIPpeakAnno)
library(dplyr)
library(tidyr)
library(reticulate)
conda_list()
conda_list()[[1]][5] %>% 
  use_condaenv(required = TRUE)

macs3_path <- '/Users/frankwellmer/miniconda3/envs/chipseq/bin/macs3'
out_directory

'''Set up where you want the output files to be made, might be 
easier to set this wherever you have the sample info sheet'''
setwd('/Users/frankwellmer/joe/chip_pipeline/')

'''Just setting up the little functions that will be used in the script'''
get_unique_comparisons <- function(data) {
  unique_values <- unique(data)
  num_values <- length(unique_values)
  
  comparisons <- list()
  
  for (i in 1:(num_values - 1)) {
    for (j in (i + 1):num_values) {
      comparisons[[length(comparisons) + 1]] <- c(unique_values[i], unique_values[j])
    }
  }
  
  return(comparisons)
}

'''Loading in the sample sheet - ensure the IPs and Controls are loaded as IP and Input, respectively'''
sample_sheet <- read.csv('./sample_sheet.csv') %>% dplyr::mutate(comparison=paste0(Genotype,'_',Condition)
)

'''All comparisons within the sample sheet'''
comps <- get_unique_comparisons(sample_sheet$comparison)

'''macs3 - peak calling and differential peak analysis'''
for(i in unique(sample_sheet$Genotype)){
    subset_sample <-  sample_sheet  %>% dplyr::filter(Genotype==i)
    for(j in unique(subset_sample$Condition)){
        data <- subset_sample %>% dplyr::filter(Condition == j)
        ips <-  data %>% dplyr::filter(Control=='IP')
        inputs <- data %>% dplyr::filter(Control=='Input')
        ip_info <-  paste(ips$bampath, collapse = " ")
        input_info <- paste(inputs$bampath, collapse = " ")
        outname <- paste0(ips$Genotype[1],'_',ips$Condition[1],'_macs3')
        macs3_command <- paste0(macs3_path,' callpeak -f BAMPE -t ', ip_info,' -c ',input_info,
        ' --gsize 119485143 --bdg -n ', outname)
        system(macs3_command)
    }
}

macs_peaks <- list.files('./',pattern='bdg')
macs3 bdgdiff -t1 treatment1.bedGraph -c1 control1.bedGraph -t2 
treatment2.bedGraph -c2 control2.bedGraph --depth1 1.0 --depth2 1.0 -o output.bedGraph --minlen 500 --maxgap 1000 --cutoff 1.0



for(i in comps){
  comp1_IP <- unique(macs_peaks[grepl(i[1], macs_peaks)])
  comp1_IP <- comp1_IP[grepl('treat', comp1_IP)]
  comp1_input <- unique(macs_peaks[grepl(i[1], macs_peaks)])
  comp1_input <- comp1_input[grepl('control', comp1_input)]
  
  comp2_IP <- unique(macs_peaks[grepl(i[2], macs_peaks)])
  comp2_IP <- comp2_IP[grepl('treat', comp2_IP)]
  comp2_input <- unique(macs_peaks[grepl(i[2], macs_peaks)])
  comp2_input <- comp2_input[grepl('control', comp2_input)]
  
  
  outname <- paste0(i[1],'_',i[2],'_macs_diffpeaks.bedgraph')
  
  sicer_command <- paste0(macs3_path, ' bdgdiff -t1 ', comp1_IP, ' -c1 ', comp1_input,' -t2 ',comp2_IP,' -c2 ',comp2_input,' -s tair10 -o ', outname ,
                          ' --depth1 1.0 --depth2 1.0 -o output.bedGraph --minlen 500 --maxgap 1000 --cutoff 1.0')
  print(sicer_command)
}


'''SICER2 - epic2 does not want to run on MAC anymore - mostly use these results as
a back-up or comparison for the macs3 results, this section will merge the bam files of each
of the replicates because it can only take sample per input and control'''

test_list <- list()



for(i in unique(sample_sheet$Genotype)){
    subset_sample <-  sample_sheet  %>% dplyr::filter(Genotype==i)
    for(j in unique(subset_sample$Condition)){
        data <- subset_sample %>% dplyr::filter(Condition == j)
        ips <-  data %>% dplyr::filter(Control=='IP')
        inputs <- data %>% dplyr::filter(Control=='Input')
        ip_info <-  paste(ips$bampath, collapse = " ")
        input_info <- paste(inputs$bampath, collapse = " ")
        merge_outname_ips <- paste0(ips$Genotype[1],'_',ips$Condition[1],'_', ips$Control[1],'_merged.bam')
        merge_command_ips <- paste0('samtools merge -o ', merge_outname_ips,' ', ip_info)

        merge_outname_input <-  paste0(inputs$Genotype[1],'_',inputs$Condition,'_', inputs$Control[1],'_merged.bam')
        merge_command_input <- paste0('samtools merge -o ', merge_outname_input,' ', input_info)

        print(merge_command_ips)
        test_list <- append(test_list, merge_outname_input)
        test_list <- append(test_list, merge_outname_ips)
    }
}


for(i in unique(sample_sheet$Genotype)){
    subset_sample <-  sample_sheet  %>% dplyr::filter(Genotype==i)
    for(j in unique(subset_sample$Condition)){
        data <- subset_sample %>% dplyr::filter(Condition == j)
        ips <-  data %>% dplyr::filter(Control=='IP')
        inputs <- data %>% dplyr::filter(Control=='Input')
        ip_info <-  paste(ips$bampath, collapse = " ")
        input_info <- paste(inputs$bampath, collapse = " ")
        merge_outname_ips <- paste0(ips$Genotype[1],'_',ips$Condition[1],'_merged.bam')
        merge_command_ips <- paste0('samtools merge -o ', merge_outname_ips,' ', ip_info)

        merge_outname_input <-  paste0(inputs$Genotype[1],'_',inputs$Condition,'_merged.bam')
        merge_command_input <- paste0('samtools merge -o ', merge_outname_input,' ', input_info)

        print(merge_command_ips)
        print(merge_command_input)
    }
}

for(i in comps){
    comp1_IP <- unique(test_list[grepl(i[1], test_list)])
    comp1_IP <- comp1_IP[!grepl('Input', comp1_IP)]
    comp1_input <- unique(test_list[grepl(i[1], test_list)])
    comp1_input <- comp1_input[grepl('Input', comp1_input)]

    comp2_IP <- unique(test_list[grepl(i[2], test_list)])
    comp2_IP <- comp2_IP[!grepl('Input', comp2_IP)]
    comp2_input <- unique(test_list[grepl(i[2], test_list)])
    comp2_input <- comp2_input[grepl('Input', comp2_input)]


    outname <- paste0(i[1],'_',i[2],'_SICER2')

    sicer_command <- paste0('sicer_df -t ', comp1_IP, ' ', comp2_IP,' -c ',comp1_input,' ',comp2_input,' -s tair10 -o ', outname )
    print(sicer_command)
}
