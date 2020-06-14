library(reticulate)
library(tidyverse)

args <- commandArgs(trailingOnly = FALSE)
print(args[6])

dataFunc <- function(x){
  myfile <- read_delim(x, delim = ",", col_names = F) %>%
    arrange(X1) %>%
    select(-X2) %>%
    separate(X1, into = c("ID","blah","Count"), sep = "\\s") %>%
    dplyr::mutate(Count = replace_na(Count, 1)) %>% # na's should be replaced by 1, this was
                                                    # realised later after running the script.
    separate(Count, into = c("blah1","blah2"), sep = "(?<=\\[)", fill = "left") %>%
    separate(blah2, into=c("ReadCount","blah3"),sep = "(?=\\])") %>%
    select(-c("blah","blah1","blah3"))
  return(myfile)
}

outfile <- str_extract(args[6], pattern = "\\d+-\\d+_.*_.*_R1_001")
mydata <- dataFunc(args[6])
numrows <- mydata %>% count() %>% pull()

source_python('CLCumi_defs.py')

mainfunc(mydata, outfile, numrows)
