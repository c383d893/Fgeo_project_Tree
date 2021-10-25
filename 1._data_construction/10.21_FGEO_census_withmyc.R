###################################################
###################################################
################## BRING IN DAT ###################
############# MYC SPECIES AND GENUS ###############
###################################################
###################################################

mycstat<-read.table("data/Funroot_species_08092021.csv",header=TRUE,sep=",")     # Species myc list

mycstat.g<-read.table("data/Funroot_genus_08092021.csv",header=TRUE,sep=",")     # Genus myc list

###################################################
###################################################
################## BRING IN DAT ###################
#################### PREP DAT #####################
###################################################
###################################################

# Read data and extract census 1 only and alive only.
dat<- readRDS("data/allcensusdata_10.2021.rds") %>%
  group_by(site) %>%
  slice_max(census) %>% # get latest census
  ungroup() %>%
  filter(status=="alive") %>%
  filter(!site %in% c("Traunstein","Amacayacu")) # remove these site for dbh and gx-gy issues
 
# Add mycorrhizal category sp, then genus:
dat.myc<-dat %>%                           
  left_join(mycstat,by="latin")                                                  # Merge serc.dat with mycstat (species)

dat.myc.g<-dat %>% 
  left_join(mycstat.g,by="genus")%>%                                             # Merge sp list (alive) with by genus
  select(c("latin","myc")) %>%                                                   # Subset just latin and myc
  rename(myc.g=myc) %>%                                                          # Rename cols
  distinct(latin,.keep_all=TRUE)                                                 # Keep only on representative from each species

dat2<-dat.myc %>% left_join(dat.myc.g, by="latin")%>%                            # Merge sp myc list with genus myc list
  mutate(myc=as.character(myc),myc.g=as.character(myc.g)) %>%                    # Convert to characters
  mutate(consensus_myc = ifelse(is.na(myc), myc.g, myc)) %>%                     # if myc is na, assign myc.g, otherwise assigned, myc
  select(-c(myc,myc.g)) %>%
  rename(myc= consensus_myc) %>%
  filter(myc=="AM"| myc=="EM") %>% # only keep two cat
  as.data.frame() %>%
  select(-c("quadrat","census","status","genus","species"))

# Find sites that don't have both myc types: can't include for conmyc heteromyc
dat2.myc<- dat2 %>% 
  group_by(site,myc,.drop=FALSE) %>% 
  distinct(latin) %>% 
  summarize(N = n())%>%
  mutate(freq = N / sum(N)) %>%
  filter(freq==1)

# Subset out those that can't run CM
dat2 <- dat2 %>% filter(!site %in% dat2.myc$site)

# Save as .rds object.
saveRDS(dat2, 'data/allcensusdata_withmyc.rds')

