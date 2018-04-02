library(tibble)
library(readr)

setwd("/Users/staskowiak/Combio_project")

#load data to tibbles
data_freqs <- read_delim("data_freq_pairs_COMPLETE.txt", delim = ";", na = "NULL")
AA_content <- read_delim("AA_content.txt", delim = ";")

#make an adjustments to get decimal point system of the percentages
AA_content[,2] <- AA_content[,2]/100
AA_content[,3] <- AA_content[,3]/100

#prepare a new dataframe
AA_names <- as.vector(names(data_freqs))
AA_names <- AA_names[4:length(nazwy)]

mean_val <- NULL
std_val <- NULL

z <- length(AA_names)

#a getting-combined-frequencies-of-AA loop
for (i in 1:z){
  
  AA <- AA_names[i]
  first <- substr(AA, 1, 1)   #split AA pairs for single letters
  second <- substr(AA, 2, 2) 
  
  std_1 <- AA_content$std[AA_content$amino_acid == first]
  std_2 <- AA_content$std[AA_content$amino_acid == second]
  
  first <- AA_content$mean_freq[AA_content$amino_acid == first]   #get the frequencies for a given letter
  second <-AA_content$mean_freq[AA_content$amino_acid == second]
  
  multifreq <- first*second   #save multiplied frequecies of a pair of AA's
  mean_val[i] <- multifreq
  std_val[i] <- std_2*std_1
}

options(scipen=999) #disable scientific notation for printing


#make a new data frame
AA_expected_freqs <- tibble(AA_name = AA_names,
                            freq = round(mean_val, 6),
                             std = round(std_val, 10))




z_scores <- data_freqs
end_col <- ncol(z_scores)
end_row <- nrow(z_scores)

for (i in 1:end_row){
  
  for (j in 4:end_col){
    
    print(paste(i, ":", j))
    
    AA_pair <- names(z_scores[j])
    observed <- z_scores[i,j]
    
    expected <- AA_expected_freqs$freq[AA_expected_freqs$AA_name == AA_pair]
    
    std <- AA_expected_freqs$std[AA_expected_freqs$AA_name == AA_pair]

    
    #getting a z-score (observed minus expeted) divided by standard dev
    z_scores[i,j] <- (observed - expected)/std
      
    
  }
  
  any(z_scores[4, 401] > 0 )
  
}




any(data_freqs[2, 401] > 0)










