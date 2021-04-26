my_data %>% 
  separate(Group, into = c("foo", "SampleID"),sep = "[/]") %>% 
  separate(SampleID, into = c("SampleID","dnainput"),sep = "[-_]") %>% 
  select(-foo) -> my_data_1

# 10,000 feet overview
my_data_1 %>%
  group_by(SampleID,dnainput) %>%
  summarise(Umitot=n(), Counttot=sum(Count),Umidistinct=n_distinct(UMI),.groups="keep") -> my_data_2

my_data_1 %>% group_by(SampleID,dnainput,UMI) %>% summarise(UMIdistinct=n(),.groups="keep") -> bar

bar %>% ungroup() %>% group_by(SampleID,dnainput) %>% summarise(UMIdistinct=n()) -> foobar

my_data_2 %>% ggplot(aes(x=SampleID,y=Umidistinct)) + geom_bar(aes(fill=dnainput),stat = "identity",position = "dodge") 
  

my_data_2 %>% ggplot(aes(x=SampleID,y=Umidistinct)) + geom_bar(stat = "identity") +
  facet_wrap(~dnainput)

my_data_2 %>% ggplot(aes(x=SampleID,y=Counttot)) + geom_bar(aes(fill=dnainput),stat = "identity",position = "dodge")

# 10,000 feet overview with atleast 20 supporting reads
my_data_1 %>%
  filter(Count >= 10) %>%
  group_by(SampleID,dnainput) %>%
  summarise(Umitot=n(), Counttot=sum(Count),Umidistinct=n_distinct(UMI),.groups="keep") -> my_data_3


my_data_3 %>% filter(dnainput==100) %>% ggplot(aes(x=SampleID,y=Umidistinct)) + geom_bar(aes(fill=dnainput),stat = "identity",position = "dodge")
  ylim(0,200000)

my_data_3 %>% ggplot(aes(x=SampleID,y=Counttot)) + geom_bar(aes(fill=dnainput),stat = "identity",position = "dodge")


# read multicov data (samples aligned to ref and counts)
my_data_multicov <- read_delim("/home/data/LabData/MPS_files/BAM_Gallery/UMIs/BAMS_phase2/multicov.tsv", 
                               col_names = T, delim = "\t")
my_data_multicov_1 <- my_data_multicov %>% select(-(Start:Locus)) %>% separate(Chrom, into = c("Chrom", "foo"),sep = 3) %>%
  select(-foo) %>%
  group_by(Chrom) %>% summarise_all(list(sum)) %>% t() %>%
  as_tibble(rownames = "Sample") %>% # convert rownames to a column
  slice(-c(1)) %>% # remove first row 
  separate(Sample, into = c("SampleID","dnainput"), sep = "[-_]") %>% 
  mutate(mulitcov=as.integer(V1)) %>% select(-V1)


# read summary stats
my_data_summ <- read_delim("data_hisp_2021/summary_stats.txt", 
             col_names = c("Sample",
                           "CSMatch",
                           "CS_PrimerMatch",
                           "CS_Primer_AnchorMatch",
                           "NoCSMatch"),delim = ",") %>%
    separate(Sample, into = c("SampleID","dnainput"),sep = "[-_]") %>%
    separate(CSMatch, into = c("foo","CSMatch"),sep = "=") %>% select(-foo) %>%
    mutate(CSMatch=as.integer(CSMatch)) %>%
    separate(CS_PrimerMatch,into = c("foo","CS_PrimerMatch"),sep = "=") %>%
    separate(CS_Primer_AnchorMatch, into = c("bar","CS_Primer_AnchorMatch"),sep = "=") %>%
    separate(NoCSMatch,into = c("foobar","NoCSMatch"),sep = "=") %>%
    select(-c(foo,bar,foobar)) %>%
    mutate(CS_PrimerMatch=as.integer(CS_PrimerMatch ),CS_Primer_AnchorMatch=as.integer(CS_Primer_AnchorMatch),
           NoCSMatch=as.integer(NoCSMatch))


mydata_2gether <- my_data_2 %>% left_join(my_data_multicov_1, by=c("SampleID"="SampleID","dnainput"="dnainput")) %>%
  left_join(my_data_summ,by = c("SampleID"="SampleID","dnainput"="dnainput"))

mydata_2gether %>% select(-c(Umitot,Counttot)) -> mydata_2gether_long

mydata_2gether_long <- mydata_2gether %>% pivot_longer(cols = Umidistinct:NoCSMatch, names_to="Type", values_to="Count")

mydata_2gether_long %>% ggplot(aes(x=SampleID,y=Count)) + geom_bar(aes(fill=Type),stat = "identity",position = "dodge") +
  facet_wrap(~dnainput, scales = "free_y")