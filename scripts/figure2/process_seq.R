#command takes *coords input for FGS region

#
#if (length(commandArgs(trailingOnly = TRUE)) != 1) {
#  print("Usage: Rscript process_seq.R\n")
#  q(status = 1)
#}

#input <- as.character(commandArgs(trailingOnly = TRUE)[1])

result_df <- data.frame(
  hap = character(),
  start = numeric(),
  stop = numeric(),
  barcode = character(),
  individual = character()
)

file_names <- Sys.glob("/scratch/groups/schumer/tris/falsegravid/pangenome/haplotypes/*coords")
#file_names <- Sys.glob("/Users/trisdodge/Desktop/Stanford/falsegravid/pangenome/barcoding/*coords")

for (input in file_names){
  #input <- "/Users/trisdodge/Desktop/Stanford/falsegravid/pangenome/barcoding/barcoding.fa_2_xbir-COAC-16-VIII-22-M_v2023.1_FGS_hap_chr-02-FGS.fa.delta.m.coords"
  
  column_names <- c("S1","E1","S2","E2","LEN1","LEN2","IDY","SEQ1","SEQ2")
  head <-read.csv(input, sep=" ", head=FALSE, nrows=1)
  mummer<-read.csv(input, sep="\t", skip=4, head=FALSE, col.names = column_names)
  mummer$HAP <- basename(head$V2)
  
  mummer_filtered <- subset(mummer, abs(S2 - min(S2[SEQ1 == "dup1"])) <= 50000 | (abs(S2 - min(S2[SEQ1 == "dup1"])) <= 65000 & LEN1 > 1000))
  
  
  unprocessed <- c("dup1", "dup2", "insertion","syntenic1", "syntenic2", "wt_insertion")
  processed <- c("inversion")
  
  #unprocessed <- c("dup1", "dup2", "insertion")
  #processed <- c("syntenic1", "syntenic2", "inversion")
  
  for (seq in unique(mummer_filtered$SEQ1)) {
    mummer_filtered_sub <- subset(mummer_filtered, SEQ1==seq)
    if (seq %in% unprocessed) {
      print(paste(length(mummer_filtered_sub$S1), "copies", seq))
      new_rows <-
        data.frame(
          hap = mummer_filtered_sub$SEQ2,
          start = mummer_filtered_sub$S2,
          stop = mummer_filtered_sub$E2,
          barcode = mummer_filtered_sub$SEQ1,
          individual = mummer_filtered_sub$HAP
        )
    } else if (seq %in% processed & mummer_filtered_sub$S2[1] < mummer_filtered_sub$E2[1]){
      print(paste(length(mummer_filtered_sub$S1), "copies", seq))
      new_rows <-
        data.frame(
          hap = mummer_filtered_sub$SEQ2,
          start = mummer_filtered_sub$S2,
          stop = mummer_filtered_sub$E2,
          barcode = mummer_filtered_sub$SEQ1,
          individual = mummer_filtered_sub$HAP
        )
    } else if (seq %in% processed & mummer_filtered_sub$S2[1] > mummer_filtered_sub$E2[1]){
      print(paste("1", "copies", seq))
      if(mummer_filtered_sub$S2[1] < mummer_filtered_sub$E2[1]){
        normal_start <- min(mummer_filtered_sub$S2)
        normal_end <- max(mummer_filtered_sub$E2)
      } else if (mummer_filtered_sub$S2[1] > mummer_filtered_sub$E2[1]){
        normal_start <- max(mummer_filtered_sub$S2)
        normal_end <- min(mummer_filtered_sub$E2)
      }else {
        normal_start <- "NA"
        normal_end <- "NA"
      }
      new_rows <-
        data.frame(
          hap = unique(mummer_filtered_sub$SEQ2),
          start = normal_start,
          stop = normal_end,
          barcode = mummer_filtered_sub$SEQ1,
          individual = mummer_filtered_sub$HAP
        )
    } else {
      print("error")}
    result_df <- rbind(result_df, new_rows)
  }
  
}

write.table(
  result_df,
  file = paste0("./", "barcoding.fa_all_coords", "_clean"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,  # Do not write column names
  quote = FALSE)
