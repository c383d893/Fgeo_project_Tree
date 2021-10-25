#### CLEANING UP MYC STATUS/ FUNGALROOT DATABASE WITH REPEATS ####

#### SPECIES VERSION #####
mycstat<-read.table("data/Funroot_species_04142021.csv",header=TRUE,sep=",") %>%    # Species myc list
  select(-myc_org) %>%                                                              # Remove col
  rename(latin=species)                                                             # Rename species to latin

# find species that have multiple designations
# number of latin-myc combinations
mycrepeat<-mycstat %>% drop_na() %>% group_by(latin,myc) %>% count()

# which latins are occur more than once (aka have different myc categories above)
mycrepeat2<-mycrepeat %>% group_by(latin) %>% count() %>% filter (n>1) 
mycnorepeat2<-mycrepeat %>% group_by(latin) %>% count() %>% filter (n==1) %>% select(-n) # no repeats

# look at those from mycrepeat 2 within mycrepeat (aka double myc categories)
mycrepeat3<-mycrepeat2 %>% select(-n) %>% left_join(mycrepeat)

# order by greatest n, and take distinct as representative 
mycrepeat4<-mycrepeat3 %>% group_by(latin) %>% arrange(desc(n)) %>% distinct(latin, .keep_all=TRUE) %>% select(-n)

#join myrepeat4 with mycnorepeat2
mycnorepeat3<- mycnorepeat2 %>% left_join(mycstat,by="latin") %>% slice(1) # add mycstat
mycstat2<-rbind(mycrepeat4, mycnorepeat3)

# sanity check: mycstat should have no counts above 1
mycstat2 %>% drop_na()%>% group_by(latin,myc) %>% count() %>% filter(n>1)

# write table: mycstat2
write_csv(mycstat2,"data/Funroot_species_08092021.csv")

#### GENUS VERSION #####

mycstat.g<-read.table("data/Funroot_genus_04142021.csv",header=TRUE,sep=",") %>%    # Genus myc list
  select(-myc_org)%>%                                                               # Remove col
  rename(genus=Genus)                                                               # Rename Genus to genus

# sanity check: mycstat should have no counts above 1
mycstat.g %>% drop_na%>% group_by(genus,myc) %>% count() %>% filter(n>1)

# write table: mycstat.g
write_csv(mycstat.g, "data/Funroot_genus_08092021.csv")

