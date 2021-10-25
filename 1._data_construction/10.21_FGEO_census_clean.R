#### Cleaning multiple census data ####
#### UPDATED Oct 22, 2021 ####

# manually select 8 relevant cols (excel)
# manually replace all " " with "_" (excel)
# manually replace all empty cells with NA (excel)
# don't filter out latin with odd species names we can still assign these to myc by genus

### Read each file
### Keep/ rename cols: dbh, gx, gy, sp, treeID, DBH, status
### Create Genus col: genus
### For now, left all undeternmined species (noted out)
### For now, kept all status tree, recoded so alive = alive, otherwise dead
### checked dbh in cm
### Checked all site areas

### Issues
#1. Sinjaraha 4: need latin names;only to check mortality form 3.
#2. SC, LP, and HSD: no quadrats
#3. San Lorenzo should have 5 censuses, only have 4
#4. Mo Singto: coming up as too large (45,70 ha)
#5. Amacayacu: status unclear

library(tidyverse)

# set global column order
col_order <- c("sp", "genus", "species","gx", "gy", "treeID","dbh","status","quadrat","census","site")

##################
## Danum Valley ##
##################
# Removed all "¬†" with "_"

DV.1<- read.table("data/multiple_census_data/Danum_Valley/Danum_Valley_Census_1_Accessed_02-17-2021_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin, quadrat=Quadrat, gx=PX, gy=PY, treeID=TreeID,dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Danum_Valley") %>%
  mutate(dbh=dbh*.10) %>%   
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("sp","","unknown mixed")) %>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

DV.2<- read.table("data/multiple_census_data/Danum_Valley/Danum_Valley_Census_2_Accessed_02-17-2021_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin,  quadrat=Quadrat, gx=PX, gy=PY, treeID=TreeID,dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=2,site="Danum_Valley") %>%
  mutate(dbh=dbh*.10) %>%    
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("sp","","unknown mixed")) %>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

##################
#### Luquillo ####
##################
# Removed all "¬†" with ""

LQ.4<- read.table("data/multiple_census_data/Luquillo/Luquillo_Census_4_Accessed_02-16-2021_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin, quadrat=Quadrat,gx=PX, gy=PY, treeID=TreeID, dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=4,site="Luquillo") %>%
  mutate(dbh=dbh*.10) %>%   
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("spp","SPECIES","X","sp.")) %>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

LQ.5<- read.table("data/multiple_census_data/Luquillo/Luquillo_Census_5_Accessed_02-16-2021_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin, quadrat=Quadrat, gx=PX, gy=PY, treeID=TreeID,dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=5,site="Luquillo") %>%
  mutate(dbh=dbh*.10) %>%   
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("spp","SPECIES","X","sp.")) %>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

LQ.6<- read.table("data/multiple_census_data/Luquillo/Luquillo_Census_6_Accessed_02-16-2021_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin, quadrat=Quadrat, gx=PX, gy=PY, treeID=TreeID,dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=6,site="Luquillo") %>%
  mutate(dbh=dbh) %>%    # No need to convert DBH, in cm
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("spp","SPECIES","X","sp.")) %>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

##################
###### SERC ######
##################

SERC.1<- read.table("data/multiple_census_data/SERC/SERC_Census_1_Accessed_02-16-2021_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin, quadrat=Quadrat, gx=PX, gy=PY, treeID=TreeID, dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="SERC") %>%
  mutate(dbh=dbh) %>%   
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("sp.","unknown")) %>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

SERC.2<- read.table("data/multiple_census_data/SERC/SERC_Census_2_Accessed_02-16-2021_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin, quadrat=Quadrat, gx=PX, gy=PY, treeID=TreeID, dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=2,site="SERC") %>%
  mutate(dbh=dbh) %>%    
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("sp.","unknown")) %>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

##################
###### SCBI ######
##################

SCBI.1<- read.table("data/multiple_census_data/SCBI/SCBI_Census_1_Accessed_02-16-2021_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin, quadrat=Quadrat, gx=PX, gy=PY, treeID=TreeID, dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="SCBI") %>%
  mutate(dbh=dbh*.10) %>%   
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("sp","unk"))%>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

SCBI.2<- read.table("data/multiple_census_data/SCBI/SCBI_Census_2_Accessed_02-16-2021_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin, quadrat=Quadrat, gx=PX, gy=PY, treeID=TreeID, dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=2,site="SCBI") %>%
  mutate(dbh=dbh*.10) %>%    # No need to convert DBH, in cm
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("sp","unk"))%>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

SCBI.3<- read.table("data/multiple_census_data/SCBI/SCBI_Census_3_Accessed_02-16-2021_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin, quadrat=Quadrat, gx=PX, gy=PY, treeID=TreeID, dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=3,site="SCBI") %>%
  mutate(dbh=dbh*.10) %>%    # No need to convert DBH, in cm
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("sp","unk"))%>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

##################
###### BCI #######
##################

BCI.1<- read.table("data/multiple_census_data/BCI/BCI_Census_1_accessed20210216_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin, quadrat=Quadrat, gx=PX, gy=PY, treeID=TreeID, dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="BCI") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("species","sp.4")%>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

BCI.2<- read.table("data/multiple_census_data/BCI/BCI_Census_2_accessed20210216_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin, quadrat=Quadrat, gx=PX, gy=PY, treeID=TreeID, dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=2,site="BCI") %>%
  mutate(dbh=dbh*.10) %>%   
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("species","sp.4")%>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

BCI.3<- read.table("data/multiple_census_data/BCI/BCI_Census_3_accessed20210216_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin, quadrat=Quadrat, gx=PX, gy=PY, treeID=TreeID, dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=3,site="BCI") %>%
  mutate(dbh=dbh*.10) %>%    
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("species","sp.4")%>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

BCI.4<- read.table("data/multiple_census_data/BCI/BCI_Census_4_accessed20210216_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin, quadrat=Quadrat, gx=PX, gy=PY, treeID=TreeID, dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=4,site="BCI") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("species","sp.4")%>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

BCI.5<- read.table("data/multiple_census_data/BCI/BCI_Census_5_accessed20210216_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin, quadrat=Quadrat, gx=PX, gy=PY, treeID=TreeID, dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=5,site="BCI") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("species","sp.4")%>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

BCI.6<- read.table("data/multiple_census_data/BCI/BCI_Census_6_accessed20210216_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin, quadrat=Quadrat, gx=PX, gy=PY, treeID=TreeID, dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=6,site="BCI") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("species","sp.4")%>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

BCI.7<- read.table("data/multiple_census_data/BCI/BCI_Census_7_accessed20210216_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin, quadrat=Quadrat, gx=PX, gy=PY, treeID=TreeID, dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=7,site="BCI") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("species","sp.4")%>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

BCI.8<- read.table("data/multiple_census_data/BCI/BCI_Census_8_accessed20210216_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin, quadrat=Quadrat, gx=PX, gy=PY, treeID=TreeID, dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=8,site="BCI") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("species","sp.4")%>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

##################
##### Cocoli #####
##################

load("data/multiple_census_data/Cocoli/marena3cns.stem1.rdata") # includes multiple sites
load("data/multiple_census_data/Cocoli/marena3cns.stem2.rdata") # includes multiple sites
load("data/multiple_census_data/Cocoli/marena3cns.stem3.rdata") # includes multiple sites
plot.1 <-marena3cns.stem1
plot.2 <-marena3cns.stem2
plot.3 <-marena3cns.stem3

load("data/multiple_census_data/Cocoli/bci.spptable.rdata")
Cocoli.sp <- bci.spptable %>%
  mutate(latin=paste(Genus,Species, sep="_")) %>%
  rename(spcode=sp) %>%
  select(c("latin","spcode"))

Cocoli.plot <- rbind(plot.1,plot.2,plot.3) %>%
  filter(plot=="cocoli") %>% # filter to only include Cocoli
  rename(spcode=sp)%>%
  left_join(Cocoli.sp, by="spcode")

Cocoli.1<- Cocoli.plot %>%
  filter(census==1) %>%
  rename(sp=latin, quadrat=quadrat, gx=gx, gy=gy, treeID=treeID, dbh=dbh, status=status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Cocoli") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("species","sp.2","sp.1")) %>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

unique(Cocoli.1$species)
Cocoli.2<- Cocoli.plot %>%
  filter(census==2) %>%
  rename(sp=latin, quadrat=quadrat, gx=gx, gy=gy, treeID=treeID, dbh=dbh, status=status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=2,site="Cocoli") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("species","sp.2","sp.1")) %>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame() 

Cocoli.3<- Cocoli.plot %>%
  filter(census==3) %>%
  rename(sp=latin, quadrat=quadrat, gx=gx, gy=gy, treeID=treeID, dbh=dbh, status=status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=3,site="Cocoli") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("species","sp.2","sp.1")) %>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()

##################
### San_Lorenzo ##
##################

load("data/multiple_census_data/San_Lorenzo__Sherman/marena4cns.stem1.rdata") # includes multiple sites
load("data/multiple_census_data/San_Lorenzo__Sherman/marena4cns.stem2.rdata") # includes multiple sites
load("data/multiple_census_data/San_Lorenzo__Sherman/marena4cns.stem3.rdata") # includes multiple sites
load("data/multiple_census_data/San_Lorenzo__Sherman/marena4cns.stem4.rdata") # includes multiple sites
plot.1 <-marena4cns.stem1
plot.2 <-marena4cns.stem2
plot.3 <-marena4cns.stem3
plot.4 <-marena4cns.stem4

load("data/multiple_census_data/San_Lorenzo__Sherman/bci.spptable.rdata")
SanLo.sp <- bci.spptable %>%
  mutate(latin=paste(Genus,Species, sep="_")) %>%
  rename(spcode=sp) %>%
  select(c("latin","spcode"))

SanLo.plot <- rbind(plot.1,plot.2,plot.3,plot.4) %>%
  filter(plot=="sherman") %>% # filter to only include San Lorenzo
  rename(spcode=sp)%>%
  left_join(Cocoli.sp, by="spcode")

SanLo.1<- SanLo.plot %>%
  filter(census==1) %>%
  rename(sp=latin, quadrat=quadrat, gx=gx, gy=gy, treeID=treeID, dbh=dbh, status=status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="San_Lorenzo") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("(sherman)","species","sp.1","sp.2","sp.3","sp.4","sp.5","sp.6")) %>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()  

SanLo.2<- SanLo.plot %>%
  filter(census==2) %>%
  rename(sp=latin, quadrat=quadrat, gx=gx, gy=gy, treeID=treeID, dbh=dbh, status=status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=2,site="San_Lorenzo") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("(sherman)","species","sp.1","sp.2","sp.3","sp.4","sp.5","sp.6")) %>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()  

SanLo.3<- SanLo.plot %>%
  filter(census==3) %>%
  rename(sp=latin, quadrat=quadrat, gx=gx, gy=gy, treeID=treeID, dbh=dbh, status=status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=3,site="San_Lorenzo") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("(sherman)","species","sp.1","sp.2","sp.3","sp.4","sp.5","sp.6")) %>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()

SanLo.4<- SanLo.plot %>%
  filter(census==4) %>%
  rename(sp=latin, quadrat=quadrat, gx=gx, gy=gy, treeID=treeID, dbh=dbh, status=status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=4,site="San_Lorenzo") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("(sherman)","species","sp.1","sp.2","sp.3","sp.4","sp.5","sp.6")) %>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()

##################
#### Sinharaja ###
##################
#!!!! 4 left to do- need latin names

SH.1<- read.table("data/multiple_census_data/Sinharaja/Sinharaja_Census_1_Accessed_03-09-2021_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin, quadrat=Quadrat, gx=PX, gy=PY, treeID=TreeID, dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Sinharaja") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("","Unidentified","unknown","*")%>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

SH.2<- read.table("data/multiple_census_data/Sinharaja/Sinharaja_Census_2_Accessed_03-09-2021_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin, quadrat=Quadrat, gx=PX, gy=PY, treeID=TreeID, stemID=StemID, dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=2,site="Sinharaja") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("","Unidentified","unknown","*")%>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

SH.3<- read.table("data/multiple_census_data/Sinharaja/Sinharaja_Census_3_Accessed_03-09-2021_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin, quadrat=Quadrat, gx=PX, gy=PY, treeID=TreeID, stemID=StemID, dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=3,site="Sinharaja") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("","Unidentified","unknown","*")%>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

load("data/multiple_census_data/Sinharaja/sinharaja.full4.rdata")

##################
## Michigan BW ###
##################
# alive = M, AL, B, R

MBW.sp<- read.table("data/multiple_census_data/Michigan_Big_Woods/species_CSD_07.31.21.txt", header=TRUE) %>%
  mutate(latin=paste(genus,species, sep="_")) %>%
  select(c("latin","spcode"))

MBW.1<- read.table("data/multiple_census_data/Michigan_Big_Woods/2003census_cortag_gxy_CSD_07.31.21.txt",header=TRUE) %>%
  left_join(MBW.sp, by="spcode") %>%
  select(-spcode) %>%
  rename(sp=latin, quadrat=quadrat, gx=gx, gy=gy, treeID=treeid, dbh=dbh, status=codes) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Michigan_Big_Woods") %>%
  mutate(dbh=dbh) %>%   
  mutate(status= ifelse(status=="M"|status=="AL"|status=="B"|status=="R", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("x","sp1")%>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

MBW.2<- read.table("data/multiple_census_data/Michigan_Big_Woods/2008census_cortag_gxy_CSD_07.31.21.txt",header=TRUE) %>%
  left_join(MBW.sp, by="spcode") %>%
  select(-spcode) %>%
  rename(sp=latin, quadrat=quadrat, gx=gx, gy=gy, treeID=treeid, dbh=dbh, status=codes) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=2,site="Michigan_Big_Woods") %>%
  mutate(dbh=dbh) %>%   
  mutate(status= ifelse(status=="M"|status=="AL"|status=="B"|status=="R", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("x","sp1")%>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

MBW.3<- read.table("data/multiple_census_data/Michigan_Big_Woods/2014census_cortag_gxy_CSD_07.31.21.txt",header=TRUE) %>%
  left_join(MBW.sp, by="spcode") %>%
  select(-spcode) %>%
  rename(sp=latin, quadrat=quadrat, gx=gx, gy=gy, treeID=treeid, dbh=dbh, status=codes) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=3,site="Michigan_Big_Woods") %>%
  mutate(dbh=dbh) %>% 
  mutate(status= ifelse(status=="M"|status=="AL"|status=="B"|status=="R", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("x","sp1")%>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

##################
### Wind River ###
##################

WR.sp<- read.table("data/multiple_census_data/Wind_River/WFDP_species_CSD_07.31.21.txt", header=TRUE) 

WR.1<- read.table("data/multiple_census_data/Wind_River/WFDP_Tree_Census1_CSD_07.31.21.csv",header=TRUE,sep=",") %>%
  left_join(WR.sp, by="SPECIES") %>%
  select(-SPECIES) %>%
  mutate(GX=UTM_X-min(UTM_X)) %>%
  mutate(GY=UTM_Y-min(UTM_Y)) %>%
  rename(sp=Latin, quadrat=QUADRAT, gx=GX, gy=GY, treeID=TREE_TAG,dbh=DBH, status=DA) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Wind_River") %>%
  mutate(dbh=dbh) %>%   
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   
  
WR.2<- read.table("data/multiple_census_data/Wind_River/WFDP_Tree_Census2_CSD_07.31.21.csv",header=TRUE,sep=",") %>%
    left_join(WR.sp, by="SPECIES") %>%
    select(-SPECIES) %>%
    mutate(GX=UTM_X-min(UTM_X)) %>%
    mutate(GY=UTM_Y-min(UTM_Y)) %>%
    rename(sp=Latin, quadrat=QUADRAT, gx=GX, gy=GY, treeID=TREE_TAG,dbh=DBH, status=DA) %>%
    select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
    separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
    mutate(sp=paste(genus,species,sep="_")   )%>%
    mutate(census=2,site="Wind_River") %>%
    mutate(dbh=dbh) %>%  
    mutate(status= ifelse(status=="A", "alive","dead")) %>%
    filter(dbh>0) %>%
    #filter(!species=="x")%>%
    #filter(!species=="sp1")%>%
    drop_na()%>%
    relocate(all_of(col_order))%>%
  as.data.frame()   

##################
#### Yosemite ####
##################
# plot burned September 1, 2013

Yos.sp<- read.table("data/multiple_census_data/Yosemite/YFDP_species_CSD_07.31.21.txt", header=TRUE)

Yos.1<- read.table("data/multiple_census_data/Yosemite/YFDP_Tree_Census1_CSD_07.31.21.txt",header=TRUE) %>%
  left_join(Yos.sp, by="SPECIES") %>%
  select(-SPECIES) %>%
  rename(sp=Latin, quadrat=QUADRAT, gx=PLOT_X, gy=PLOT_Y, treeID=TREE_TAG,dbh=DBH, status=DA) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Yosemite") %>%
  mutate(dbh=dbh) %>%   
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   
  
Yos.2<- read.table("data/multiple_census_data/Yosemite/YFDP_Tree_Census2_CSD_07.31.21.txt",header=TRUE) %>%
  left_join(Yos.sp, by="SPECIES") %>%
  select(-SPECIES) %>%
  rename(sp=Latin, quadrat=QUADRAT, gx=PLOT_X, gy=PLOT_Y, treeID=TREE_TAG,dbh=DBH, status=DA) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=2,site="Yosemite") %>%
  mutate(dbh=dbh) %>%    
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   
  
Yos.3<- read.table("data/multiple_census_data/Yosemite/YFDP_Tree_Census2_CSD_07.31.21.txt",header=TRUE) %>%
  left_join(Yos.sp, by="SPECIES") %>%
  select(-SPECIES) %>%
  rename(sp=Latin, quadrat=QUADRAT, gx=PLOT_X, gy=PLOT_Y, treeID=TREE_TAG, dbh=DBH, status=DA) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=3,site="Yosemite") %>%
  mutate(dbh=dbh) %>%    
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()   
  
##################
###### Utah ######
##################
# tag taken as treeID

Utah.sp<- read.table("data/multiple_census_data/Utah/UFDP_species_CSD_07.31.21.txt", header=TRUE) 

Utah.1<- read.table("data/multiple_census_data/Utah/UFDP_Tree_Census1_CSD_07.31.21.csv",header=TRUE,sep=",") %>%
  left_join(Utah.sp, by="SPECIES") %>%
  select(-SPECIES) %>%
  mutate(GX=UTM_X-min(UTM_X)) %>%
  mutate(GY=UTM_Y-min(UTM_Y)) %>%
  rename(sp=Latin, quadrat=QUADRAT, gx=GX, gy=GY, treeID=TREE_TAG,dbh=DBH, status=DA) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Utah") %>%
  mutate(dbh=as.numeric(as.character(dbh)),gx=as.numeric(as.character(gx)),gy=as.numeric(as.character(gy))) %>%    #  convert to numeric
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame() 

Utah.2<- read.table("data/multiple_census_data/Utah/UFDP_Tree_Census2_CSD_07.31.21.csv",header=TRUE,sep=",") %>%
  left_join(Utah.sp, by="SPECIES") %>%
  select(-SPECIES) %>%
  mutate(GX=UTM_X-min(UTM_X)) %>%
  mutate(GY=UTM_Y-min(UTM_Y)) %>%
  rename(sp=Latin, quadrat=QUADRAT, gx=GX, gy=GY, treeID=TREE_TAG,dbh=DBH, status=DA) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=2,site="Utah") %>%
  mutate(dbh=as.numeric(as.character(dbh)),gx=as.numeric(as.character(gx)),gy=as.numeric(as.character(gy))) %>%    #  convert to numeric
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  drop_na()%>%
  relocate(all_of(col_order))%>%
  as.data.frame()  

##################
#### Wabikon #####
##################
# tag taken as treeID

Wab.sp<- read.table("data/multiple_census_data/Wabikon/WabikonTreeSpeciesList_v20180801_CSD_07.31.21.txt", header=TRUE) %>%
  rename(sp=SpCode)

Wab.123<- read.table("data/multiple_census_data/Wabikon/Wabikon123_CTFS.format_20210201_CSD_07.31.21.txt",header=TRUE) %>%
  left_join(Wab.sp, by="sp") %>%
  select(-sp) %>%
  rename(sp=Latin, quadrat=quadrat, gx=gx, gy=gy, treeID=tag, dbh=dbh, status=DFstatus) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat","CensusID")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(site="Wabikon") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  drop_na() 
  
  
Wab.1<- Wab.123 %>%
        filter(CensusID=="WAB1")%>%
        select(-CensusID) %>%
        mutate(census = "1")%>%
        relocate(all_of(col_order)) %>%
       as.data.frame()   

Wab.2<- Wab.123 %>%
  filter(CensusID=="WAB2")%>%
  select(-CensusID) %>%
  mutate(census ="2")%>%
  relocate(all_of(col_order)) %>%
  as.data.frame()   

Wab.3<- Wab.123 %>%
  filter(CensusID=="WAB3")%>%
  select(-CensusID) %>%
  mutate(census ="3")%>%
  relocate(all_of(col_order)) %>%
  as.data.frame()   

## HSD ##
# tag taken as treeID

HSD.1<- read.table("data/multiple_census_data/Heishiding/HSD_CSD_7.31.21.txt",header=TRUE) %>%
  rename(sp=latin,  gx=gx, gy=gy, treeID=tag, dbh=dbh1, status=status1) %>%
  select(c("sp","gx","gy","treeID","dbh","status")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Heishiding") %>%
  mutate(dbh=dbh) %>%   
  mutate(status= ifelse(status=="Alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  drop_na() %>%
  add_column(quadrat = NA) %>%# add empty quadrat col
  relocate(all_of(col_order))%>%
  as.data.frame()   

HSD.2<- read.table("data/multiple_census_data/Heishiding/HSD_CSD_7.31.21.txt",header=TRUE) %>%
  rename(sp=latin,  gx=gx, gy=gy, treeID=tag, dbh=dbh2, status=status2) %>%
  select(c("sp","gx","gy","treeID","dbh","status")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=2,site="Heishiding") %>%
  mutate(dbh=dbh) %>%   
  mutate(status= ifelse(status=="Alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  drop_na() %>%
  add_column(quadrat = NA) %>%# add empty quadrat col
  relocate(all_of(col_order))%>%
  as.data.frame()  

##################
## Scotty Creek ##
##################

SC.1<- read.table("data/multiple_census_data/Scotty_Creek/scottycreek_data_3sep2020_CSD_7.31.21.txt",header=TRUE) %>%
  rename(sp=species,  gx=gx, gy=gy, treeID=tag, dbh=DBH1, status=status1) %>%
  select(c("sp","gx","gy","treeID","dbh","status")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Scotty_Creek") %>%
  mutate(dbh=dbh) %>%   
  mutate(status= ifelse(status=="1", "alive","dead")) %>%
  filter(dbh>0) %>%
  drop_na() %>%
  add_column(quadrat = NA) %>%# add empty quadrat col
  relocate(all_of(col_order))%>%
  as.data.frame()   

SC.2<- read.table("data/multiple_census_data/Scotty_Creek/scottycreek_data_3sep2020_CSD_7.31.21.txt",header=TRUE) %>%
  rename(sp=species,  gx=gx, gy=gy, treeID=tag, dbh=DBH2, status=status2) %>%
  select(c("sp","gx","gy","treeID","dbh","status")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=2,site="Scotty_Creek") %>%
  mutate(dbh=dbh) %>%   
  mutate(status= ifelse(status=="1", "alive","dead")) %>%
  filter(dbh>0) %>%
  drop_na() %>%
  add_column(quadrat = NA) %>% # add empty quadrat col
  relocate(all_of(col_order))%>%
  as.data.frame()   

##################
##### Zofin ######
##################
# alive = AB, AW

Zof.4<- read.table("data/multiple_census_data/Zofin/Zofin_Census_4_Accessed_03-03-2021_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin,  quadrat=Quadrat,gx=PX, gy=PY, treeID=TreeID, dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=4,site="Zofin") %>%
  mutate(dbh=dbh*.10) %>%   
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

Zof.5<- read.table("data/multiple_census_data/Zofin/Zofin_Census_5_Accessed_03-03-2021_CSD_07.31.21.txt",header=TRUE) %>%
  rename(sp=Latin,  quadrat=Quadrat,gx=PX, gy=PY, treeID=TreeID, dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=5,site="Zofin") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

##################
##### Korup ######
##################

load("data/multiple_census_data/Korup/korup.full1.rdata")
Kor.plot1<-korup.full1

load("data/multiple_census_data/Korup/korup.full2.rdata")
Kor.plot2<-korup.full2

Kor.sp<-read.table("data/multiple_census_data/Korup/KFDP_spptable_JAL_Updated_CDedit_10.06.2021.csv", header=TRUE,sep=",")

Kor.1<-  Kor.plot1%>%
  rename(Sp.code=sp)%>%
  left_join(Kor.sp, by="Sp.code") %>%
  rename(sp=Latin,  quadrat=quadrat,gx=gx, gy=gy, treeID=treeID, dbh=dbh, status=status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep=" ", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Korup") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("sp._1","sp._2","sp._3","sp._4","sp._5","sp._6","sp._7","sp._8","sp._9","sp._10",
  #                       "sp.","sp._nov_1.","","<NA>","cf","cf.","*","sp.2","cf._goldieana","sp._nov.")) %>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  as.data.frame()  

Kor.2<-  Kor.plot2%>%
  rename(Sp.code=sp)%>%
  left_join(Kor.sp, by="Sp.code") %>%
  rename(sp=Latin,  quadrat=quadrat,gx=gx, gy=gy, treeID=treeID, dbh=dbh, status=status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep=" ", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=2,site="Korup") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("sp._1","sp._2","sp._3","sp._4","sp._5","sp._6","sp._7","sp._8","sp._9","sp._10",
  #                       "sp.","sp._nov_1.","","<NA>","cf","cf.","*","sp.2","cf._goldieana","sp._nov.")) %>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  as.data.frame()   

######################
##### La Planada #####
######################
# alive = A or B

LaPl.sp<- read.table("data/multiple_census_data/La_Planada/codes_spp_CSD_8.2.21.txt", header=TRUE)

LaPl.1<- read.table("data/multiple_census_data/La_Planada/Censo_1_La Planada_CSD_8.2.21.txt",header=TRUE) %>%
  left_join(LaPl.sp, by="sp") %>%
  select(-sp) %>%
  rename(sp=Latin, gx=gx, gy=gy, treeID=tag, dbh=dbh, status=status) %>%
  select(c("sp","gx","gy","treeID","dbh","status")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="La_Planada") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="A"| status=="B", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species=="NA")%>%
  drop_na() %>%
  add_column(quadrat = NA) %>% # add empty quadrat col
  relocate(all_of(col_order))%>%
  as.data.frame()    
  
LaPl.2<- read.table("data/multiple_census_data/La_Planada/Censo_2_La Planada_CSD_8.2.21.txt",header=TRUE) %>%
  left_join(LaPl.sp, by="sp") %>%
  select(-sp) %>%
  rename(sp=Latin, gx=gx, gy=gy, treeID=tag, dbh=dbh, status=status) %>%
  select(c("sp","gx","gy","treeID","dbh","status")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=2,site="La_Planada") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="A"| status=="B", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species=="NA")%>%
  drop_na() %>%
  add_column(quadrat = NA) %>% # add empty quadrat col
  relocate(all_of(col_order))%>%
  as.data.frame()   

#####################
##### Mo Singto #####
#####################

# TBD
MoSing.sp<-read.table("data/single_census_data/Mo_Singto/MoSingto_SpeciesTable_20160719_CDedit.csv", header=TRUE, sep=",") 

MoSing.plot<-read.table("data/single_census_data/Mo_Singto/MoSingto_Census_1_2_3_20160719_CDedit.csv", header=TRUE,sep=",")

MoSing.1<- MoSing.plot %>%        # Read plot dat  
  left_join(MoSing.sp, by="speciesid") %>%
  mutate(Latin=paste(genus,species,sep="_")   ) %>%
  rename(sp=Latin, gx=X, gy=Y, treeID=tag, dbh=dbh_C1, status=status_1, quadrat=quadrat) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Mo_Singto") %>%
  mutate(dbh=dbh) %>%   
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("cf. integrifolia","sp"))%>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  as.data.frame() 

MoSing.2<- MoSing.plot %>%        # Read plot dat  
  left_join(MoSing.sp, by="speciesid") %>%
  mutate(Latin=paste(genus,species,sep="_")   ) %>%
  rename(sp=Latin, gx=X, gy=Y, treeID=tag, dbh=dbh_C2, status=status_2, quadrat=quadrat) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=2,site="Mo_Singto") %>%
  mutate(dbh=dbh) %>%   
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("cf. integrifolia","sp"))%>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  as.data.frame()

MoSing.3<- MoSing.plot %>%        # Read plot dat  
  left_join(MoSing.sp, by="speciesid") %>%
  mutate(Latin=paste(genus,species,sep="_")   ) %>%
  rename(sp=Latin, gx=X, gy=Y, treeID=tag, dbh=dbh_C3, status=status_3, quadrat=quadrat) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=3,site="Mo_Singto") %>%
  mutate(dbh=dbh) %>%   
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("cf. integrifolia","sp"))%>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  as.data.frame()
#### Cleaning single census data ####
#### Aug 2, 2021 ####

##################
### Amacayacu ####
##################

Amac.1<- read.table("data/single_census_data/Amacayacu/Amacayacu1_CDedit_10.05.21.csv", header=TRUE, sep=",") %>%
  select(-sp) %>%
  rename(sp=Species, gx=gx, gy=gy, treeID=tag, dbh=dbh.cm, status=Ha, quadrat=quadrat) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Amacayacu") %>%
  mutate(dbh=dbh) %>%   
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("sp,1","sp,2","sp,3","sp,4","sp,5","sp,6","sp,7","sp,8","sp,9",
  #                       "sp,10","sp,11","sp,12","sp,13","sp,14","sp,15","sp,16","sp,17",
  #                       "sp,17","sp,18","sp,19","sp,nov,02","sp,nov,","sp,"))#%>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  mutate(treeID= gsub("_", "", treeID)) %>%
  as.data.frame()  

##################
###### TRC #######
##################

TRC.1<-read.table("data/single_census_data/Tyson_Research_Center/trc_PlotDataReport01-13-2021_CSD.txt", header=TRUE) %>%        # Read plot dat  
  rename(sp=Latin, gx=PX, gy=PY, treeID=TreeID, dbh=DBH, status=Status, quadrat=Quadrat) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Tyson_Research_Center") %>%
  mutate(dbh=dbh) %>%   
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species=="NA")%>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  as.data.frame()  

######################
### Harvard Forest ###
######################

load("data/single_census_data/Harvard_Forest/harvardforest.full1.rdata")
hf.plot<-harvardforest.full1
hf.tax<-read.table("data/single_census_data/Harvard_Forest/hf253-02-species-codes_nospace.csv", header=TRUE,sep=",")  # Read tax dat
HF.1<-hf.plot %>%        # Read plot dat  
  left_join(hf.tax, by= "sp") %>%                                                  # Merge with tax
  select(-sp) %>%
  rename(sp=latin, gx=gx, gy=gy, treeID=treeID, dbh=dbh, status=status, quadrat=quadrat) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Tyson_Research_Center") %>%
  mutate(dbh=dbh) %>%   
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("spp.","unk))%>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  as.data.frame()  

######################
####### LDW ##########
######################

load("data/single_census_data/Little_Dickie_Woods/ldw.stem1.rdata")
ldw.plot<-ldw.stem1
load("data/single_census_data/Little_Dickie_Woods/ldw.spptable.rdata")
ldw.tax<- ldw.spptable %>% mutate(sp = as.character(sp))
LDW.1<-ldw.plot %>%        # Read plot dat  
  left_join(ldw.tax, by= "sp") %>%   # Merge with tax
  select(-sp) %>%
  rename(sp=Latin, gx=gx, gy=gy, treeID=treeID, dbh=dbh, status=status, quadrat=quadrat) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  mutate_if(is.character, str_replace_all, pattern = " ", replacement = "_")  %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Little_Dickie_Woods") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species=="NA")%>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  as.data.frame()  

######################
##### Traunstein #####
######################

# Read plot, tax, join and simplify:
ts.tax<-read.table("data/single_census_data/Traunstein/TS_TaxonomyReport07-09-2021_1152371225_CSD.txt", header=TRUE)            # Read tax dat
TS.1.5cm<- # read in low cm values
TS.1<-read.table("data/single_census_data/Traunstein/TS_PlotDataReport07-09-2021_918548146_CSD.txt", header=TRUE) %>%        # Read plot dat  
  left_join(ts.tax, by= "Mnemonic") %>%                                                  # Merge with tax
  rename(sp = Latin, dbh = DBH, gx = PX, gy = PY, treeID = TreeID, quadrat = Quadrat, status=Status) %>% # Rename cols
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Traunstein") %>%
  mutate(dbh=dbh*.10) %>%  
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  filter(dbh>0) %>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  as.data.frame()  

######################
#### Indian Cave #####
######################

IC_sp<-read.table("data/single_census_data/Indian_Cave/20210318_ICSPWoodySpecies_CSD_8.2.21.txt", header=TRUE)

IC.1<- read.table("data/single_census_data/Indian_Cave/20210305_ICSP_main-stem-data_CSD_8.2.21.txt",header=TRUE) %>%
  left_join(IC_sp, by="Species_Code") %>%
  select(-Species_Code) %>%
  rename(sp=Latin_binomial, gx=x, gy=y,quadrat=Quadrat, treeID=PtID, dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Indian_Cave") %>%
  mutate(dbh=dbh) %>%    
  filter(status=="LI") %>%
  filter(dbh>0) %>%
  #filter(!species=="NA")%>%
  mutate(status= ifelse(status=="alive", "alive","dead")) %>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  as.data.frame() 

######################
###### Niobrara ######
######################

Nio_sp<-read.table("data/single_census_data/Niobrara/20200325_NVPWoodySpecies_CSD_8.2.21.txt", header=TRUE)

Nio.1<- read.table("data/single_census_data/Niobrara/20210223-NVP_main-stem-data_CSD_8.2.21.txt",header=TRUE) %>%
  left_join(Nio_sp, by="Species_Code") %>%
  select(-Species_Code) %>%
  rename(sp=Latin_binomial, gx=x, gy=y,quadrat=Quadrat, treeID=PtID, dbh=DBH, status=Status) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Niobrara") %>%
  mutate(dbh=dbh) %>%    
  mutate(status= ifelse(status=="LI", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species=="NA")%>%
  drop_na()%>%    
  relocate(all_of(col_order))%>%
  as.data.frame()  

#####################
###### Fushan #######
#####################

Fushan_sp<-read.table("data/single_census_data/Fushan/fushan_sp_list_CDedit.csv", header=TRUE, sep=",") %>%
  select(c("spcode","Latin"))

load("data/single_census_data/Fushan/fushan.full2.rdata")
Fushan.1<- fushan.full2 %>%        # Read plot dat 
  rename(spcode=sp) %>%
  left_join(Fushan_sp, by="spcode") %>%
  rename(sp=Latin, gx=gx, gy=gy, treeID=treeID, dbh=dbh, status=status, quadrat=quadrat) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Fushan") %>%
  mutate(dbh=dbh) %>%   
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species=="NA")%>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  as.data.frame()  

#####################
## Huai Kha Khaeng ##
#####################

HKK_sp<-read.table("data/single_census_data/HKK/HKK_splist_CDedit.csv", header=TRUE, sep=",") %>%
  select(c("spcode","latin"))

load("data/single_census_data/HKK/hkk.full4.rdata")
HKK.1<- hkk.full4 %>%        # Read plot dat  
  rename(spcode=sp) %>%
  left_join(HKK_sp, by="spcode") %>%
  rename(sp=latin, gx=gx, gy=gy, treeID=treeID, dbh=dbh, status=status, quadrat=quadrat) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep="_", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Huai_Kha_Khaeng") %>%
  mutate(dbh=dbh*.10) %>%   
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("sp.","cf.","sp.1))%>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  as.data.frame() 

#####################
##### Khao Chong ####
#####################

load("data/single_census_data/Khao_Chong/khaochong.spptable.rdata")
KhaoC.sp<- khaochong.spptable %>%
  rename(spcode=sp) %>%
  select(c("spcode","Latin"))

load("data/single_census_data/Khao_Chong/khaochong.full3.rdata")
KhoaC.plot<-khaochong.full3

KhoaC.1<- KhoaC.plot %>%        # Read plot dat  
  rename(spcode=sp) %>%
  left_join(KhaoC.sp, by="spcode") %>%
  rename(sp=Latin, gx=gx, gy=gy, treeID=treeID, dbh=dbh, status=status, quadrat=quadrat) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep=" ", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Khoa_Chong") %>%
  mutate(dbh=dbh*.10) %>%   
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("Indet","sp.3","sp.2","cf.","sp."))%>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  as.data.frame() 

#####################
#### Laupahoehoe ####
#####################

load("data/single_census_data/Laupahoehoe/Laupahoehoe.spptable.rdata")
Laup.sp<- Laupahoehoe.spptable %>%
  rename(spcode=sp)%>%
  select(c("spcode","Latin"))

load("data/single_census_data/Laupahoehoe/Laupahoehoe.full5.rdata")
Laup.plot<-Laupahoehoe.full5

Laup.1<- Laup.plot %>%        # Read plot dat  
  rename(spcode=sp) %>%
  left_join(Laup.sp, by="spcode") %>%
  rename(sp=Latin, gx=gx, gy=gy, treeID=treeID, dbh=dbh, status=status, quadrat=quadrat) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep=" ", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Laupahoehoe") %>%
  mutate(dbh=dbh) %>%   
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species=="sp.")%>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  as.data.frame() 

#####################
#### Lienhuachih ####
#####################

LHC.sp<- read.table("data/single_census_data/LHC/lhc_sp_list_CDedit.csv", header=TRUE,sep=",")

load("data/single_census_data/LHC/lhc_census.rdata")
LHC.plot<-census1

LHC.1<- LHC.plot %>%        # Read plot dat  
  rename(spcode=sp) %>%
  left_join(LHC.sp, by="spcode") %>%
  rename(sp=Latin, gx=gx, gy=gy, treeID=treeID, dbh=dbh, status=status, quadrat=quadrat) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep=" ", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Lienhuachih") %>%
  mutate(dbh=dbh) %>%   
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species=="sp.")%>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  as.data.frame() 

#####################
##### Palamanui #####
#####################

load("data/single_census_data/Palamanui/Palamanui.spptable.rdata")
Palam.sp<-Palamanui.spptable %>%
  rename(spcode=sp)

load("data/single_census_data/Palamanui/Palamanui.full5.rdata")
Palam.plot<-Palamanui.full5

Palam.1<- Palam.plot %>%        # Read plot dat  
  rename(spcode=sp) %>%
  left_join(Palam.sp, by="spcode") %>%
  rename(sp=Latin, gx=gx, gy=gy, treeID=treeID, dbh=dbh, status=status, quadrat=quadrat) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep=" ", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Palamanui") %>%
  mutate(dbh=dbh) %>%   
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species=="sp.")%>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  as.data.frame() 

#####################
###### Palan ########
#####################

load("data/single_census_data/Palanan/palanan.spptable.rdata")
Palan.sp<-palanan.spptable %>%
  rename(spcode=sp)

load("data/single_census_data/Palanan/palanan.full4.rdata")
Palan.plot<-palanan.full4

Palan.1<- Palan.plot %>%        # Read plot dat  
  rename(spcode=sp) %>%
  left_join(Palan.sp, by="spcode") %>%
  rename(sp=Latin, gx=gx, gy=gy, treeID=treeID, dbh=dbh, status=status, quadrat=quadrat) %>%
  rename(genus=Genus, species=Species) %>%
  select(c("genus","species","gx","gy","treeID","dbh","status","quadrat")) %>%
  mutate(sp=paste(genus,species,sep="_"))%>%
  mutate(census=1,site="Palanan") %>%
  mutate(dbh=dbh) %>%   
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("cf.","sp.","cf.ruficaulis","hairy"))%>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  as.data.frame() 

#####################
####### Rabi ########
#####################

load("data/single_census_data/Rabi/rabisp.rdata")
Rabi.sp<-rabisp %>%
  rename(spcode=sp)

load("data/single_census_data/Rabi/rabifull.rdata")
Rabi.plot<-rabifull

Rabi.1<- Rabi.plot %>%        # Read plot dat  
  rename(spcode=sp) %>%
  left_join(Rabi.sp, by="spcode") %>%
  rename(sp=Latin, gx=gx, gy=gy, treeID=treeID, dbh=dbh, status=status, quadrat=quadrat) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep=" ", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Rabi") %>%
  mutate(dbh=dbh*.10) %>%   
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  #filter(!species %in% c("cf.","sp.","sp.1","sp.2","sp.3","sp.4","sp.5","sp.6","sp.7","sp.8","NA","Unknown"))%>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  as.data.frame() 

#####################
### Speulderbos #####
#####################

load("data/single_census_data/Speulderbos/speulderbos.spptable.rdata")
Speul.sp<-speulderbos.spptable %>%
  rename(spcode=sp)

load("data/single_census_data/Speulderbos/speulderbos.full1.rdata")
Speul.plot<-speulderbos.full1

Speul.1<- Speul.plot %>%        # Read plot dat  
  rename(spcode=sp) %>%
  left_join(Speul.sp, by="spcode") %>%
  rename(sp=Latin, gx=gx, gy=gy, treeID=treeID, dbh=dbh, status=status, quadrat=quadrat) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep=" ", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Speulderbos") %>%
  mutate(dbh=dbh) %>%   
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  filter(!species=="x")%>%
  filter(!species=="sp.")%>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  as.data.frame() 

#####################
###### Wanang #######
#####################

load("data/single_census_data/Wanang/wanang.spptable.rdata")
Wan.sp<-wanang.spptable %>%
  rename(spcode=sp)

load("data/single_census_data/Wanang/wanang.full1.rdata")
Wan.plot<-wanang.full1

Wan.1<- Wan.plot %>%        # Read plot dat  
  rename(spcode=sp) %>%
  left_join(Wan.sp, by="spcode")  %>%
  rename(sp=Latin, gx=gx, gy=gy, treeID=treeID, dbh=dbh, status=status, quadrat=quadrat) %>%
  select(c("sp","gx","gy","treeID","dbh","status","quadrat")) %>%
  separate(sp, c("genus", "species"), sep=" ", remove=TRUE) %>%
  mutate(sp=paste(genus,species,sep="_")   )%>%
  mutate(census=1,site="Wanang") %>%
  mutate(dbh=dbh) %>%   
  mutate(status= ifelse(status=="A", "alive","dead")) %>%
  filter(dbh>0) %>%
  filter(!species=="")%>%
  drop_na() %>%
  relocate(all_of(col_order))%>%
  as.data.frame() 


#################################
#### Joining all census data ####
#################################

allcensusdata<-rbind(DV.1,DV.2,LQ.4,LQ.5, LQ.6,SERC.1,SERC.2,SCBI.1,SCBI.2, SCBI.3, 
      BCI.1,BCI.2,BCI.3,BCI.4,BCI.5,BCI.6,BCI.7,BCI.8,Cocoli.1, Cocoli.2, Cocoli.3, 
      SanLo.1, SanLo.2, SanLo.3, SanLo.4,SH.1,SH.2,SH.3,MBW.1,MBW.2,MBW.3,WR.1,WR.2, 
      Yos.1,Yos.2,Yos.3,Utah.1,Utah.2, Wab.1,Wab.2,Wab.3,HSD.1,HSD.2,SC.1,SC.2, 
      Zof.4,Zof.5,Kor.1, Kor.2,LaPl.1,LaPl.2, MoSing.1, MoSing.2, MoSing.3, 
      Amac.1, TRC.1, HF.1,LDW.1,TS.1,IC.1, Nio.1, Fushan.1, HKK.1, KhoaC.1, 
      Laup.1, LHC.1, Palam.1, Palan.1, Rabi.1, Speul.1, Wan.1) %>% 
      rename(latin=sp) %>%  # rename sp as latin
      mutate(gx=as.numeric(as.character(gx)),gy=as.numeric(as.character(gy))) # make all gx and gy numeric

# write table

saveRDS(allcensusdata, "data/allcensusdata_10.2021.rds")

####CHECK DBH
range(DV.1$dbh) # *.10
range(DV.2$dbh) # *.10
range(LQ.4$dbh) # *.10
range(LQ.5$dbh) # *.10
range(LQ.6$dbh) # *.10
range(SERC.1$dbh) # ok
range(SERC.2$dbh) # ok
range(SCBI.1$dbh) # *.10
range(SCBI.2$dbh) # *.10
range(SCBI.3$dbh) # *.10
range(BCI.1$dbh) # *.10
range(BCI.2$dbh) # *.10
range(BCI.3$dbh) # *.10
range(BCI.4$dbh) # *.10
range(BCI.5$dbh) # *.10
range(BCI.6$dbh) # *.10
range(BCI.7$dbh) # *.10
range(BCI.8$dbh) # *.10
range(Cocoli.1$dbh) # *.10
range(Cocoli.2$dbh) # *.10
range(Cocoli.3$dbh) # *.10
range(SanLo.1$dbh) # *.10
range(SanLo.2$dbh) # *.10
range(SanLo.3$dbh) # *.10
range(SanLo.4$dbh) # *.10
range(SH.1$dbh) # *.10
range(SH.2$dbh) # *.10
range(SH.3$dbh) # *.10
range(MBW.1$dbh) # ok
range(MBW.2$dbh) # ok
range(MBW.3$dbh) # ok
range(WR.1$dbh) # ok
range(WR.2$dbh) # ok
range(Yos.1$dbh) # ok
range(Yos.2$dbh) # ok
range(Yos.3$dbh) # ok
range(Utah.1$dbh)# ok
range(Utah.2$dbh)# ok
range(Wab.1$dbh) # *.10
range(Wab.2$dbh) # *.10
range(Wab.3$dbh) # *.10
range(HSD.1$dbh)# ok
range(HSD.2$dbh)# ok
range(SC.1$dbh)# ok
range(SC.2$dbh)# ok
range(Zof.4$dbh)# *.10
range(Zof.5$dbh)# *.10
range(Kor.1$dbh)# *.10
range(Kor.2$dbh)# *.10
range(LaPl.1$dbh) # *.10
range(LaPl.2$dbh) # *.10
range(MoSing.1$dbh) #ok
range(MoSing.2$dbh) #ok
range(MoSing.3$dbh) #ok
range(Amac.1$dbh) # ok
range(TRC.1$dbh) # ok
range(HF.1$dbh) # ok
range(LDW.1$dbh) # *.10
range(TS.1$dbh) # *.10
range(IC.1$dbh) # ok
range(Nio.1$dbh) #ok
range(Fushan.1$dbh) # ok
range(HKK.1$dbh) # *.10
range(KhoaC.1$dbh) # *.10
range(Laup.1$dbh) # ok 
range(LHC.1$dbh) # ok 
range(Palam.1$dbh) # ok 
range(Palan.1$dbh) # ok 
range(Rabi.1$dbh) # *.10
range(Speul.1$dbh) # ok 
range(Wan.1$dbh) # ok 

### CHECK X-Y RANGE/AREA
Lx = max(lp.dat$gx) ; Ly = max(lp.dat$gy); Area = (Lx*Ly)/10000	  # East-west width (Lx) and north-south height (Ly) of plot along with number of total ha of plot

### CHECK unidentified species to remove