## Relative abundance plots code for generating the results related with the 18S dataset generated and presented in 
## the manuscript entitled "Unveiling a rich and functionally diverse soil funga in the endangered South Brazilian 
## Highland Grasslands for promoting conservation".

#autor: Kelmer Martins-Cunha
#contact: kelmermartinscunha@gmail.com


require(ggplot2)
require(dplyr)
require(tidyr)
require(glue)
require(httr)
require(jsonlite)
require(stringr)


set.seed(2134567, kind = 'Mersenne-Twister')

#### Phylum + order ####
ddOTUtable <- read.csv('./data/18s/ddOTU_table.csv', header=TRUE)
rownames(ddOTUtable) <- ddOTUtable$X
ddOTUtable$X <- NULL
rownames(ddOTUtable) <- gsub('.blast', '', as.character(rownames(ddOTUtable)))
ddOTUtable <- as.data.frame(ddOTUtable)
ddOTUtable$sample <- rownames(ddOTUtable)
sh_list <- as.list(colnames(ddOTUtable[1:(length(ddOTUtable)-1)]))
ddOTUtable <- ddOTUtable %>% 
  pivot_longer(cols=-sample)


sh_tax <- read.csv('./data/18s/ddOTU_taxonomy.csv')
correct <- sh_tax %>% 
  filter(grepl('les',order))
incorrect <- sh_tax %>% 
  filter(!grepl('les',order)) %>% 
  filter(!order%in%c('Incertae Sedis','uncultured', 'uncultured fungus', 
                     'metagenome', 'unidentified', 'fungal', ' ', '',
                     'uncultured Chytridiomycota', 'uncultured Basidiomycota',
                     'uncultured eukaryote', 'uncultured Glomus', 'uncultured Glomeromycotina',
                     'eukaryotic picoplankton environmental sample', 'uncultured Ascomycota',
                     'uncultured marine eukaryote','uncultured glomeraceous AM fungus',
                     'uncultured Eimeriidae', 'uncultured fungal contaminant',
                     'uncultured Plakinidae sp.', 'uncultured Tremellaceae', 'uncultured Holtermanniella',
                     'uncultured Sagenomella')) %>% 
  na.exclude()
incorrect$order <- as.data.frame(sub(" .*", "", incorrect$order))
#write.csv(incorrect, 'incorrect_sh_tax.csv')
# Passed on species matching GBIF tool (https://www.gbif.org/tools/species-lookup)
incorrect_updated <- read.csv('./data/18s/incorrect_sh_tax_updated.csv', sep=';')
incorrect_updated <- select(incorrect_updated, -X.1)

sh_tax <- rbind(correct, incorrect_updated)


## Checking phyla relative abundance numbers
phyla_18 <- sh_tax %>% 
  filter(sh %in% sh_list) %>% 
  filter(!phylum=='') %>% 
  select(phylum, sh) %>% 
  distinct() %>% 
  filter(!phylum%in%c('Incertae Sedis'))

phyla_18 <- left_join(ddOTUtable, phyla_18, by=c('name'='sh'))
phyla_18 <- phyla_18 %>%
  na.exclude() %>%
  group_by(sample) %>% 
  mutate(percentage = value/sum(value)*100) %>% 
  group_by(sample, phylum) %>% 
  summarise(percentage = sum(percentage)) %>%
  group_by(phylum) %>% 
  summarise(perc=mean(percentage)) ## Mean relative abundance in samples

## Phylum plot
phyla_18 <- sh_tax %>% 
  filter(sh %in% sh_list) %>% 
  filter(!phylum=='') %>% 
  select(phylum, sh) %>% 
  distinct() %>% 
  filter(!phylum%in%c('Incertae Sedis'))

phyla_18 <- left_join(ddOTUtable, phyla_18, by=c('name'='sh'))
phyla_18 <- phyla_18 %>%
  na.exclude() %>%
  mutate(phylum = case_when(
    phylum == "Entomophthoromycotina" ~ "Entomophthoromycota",
    phylum == "Glomeromycotina" ~ "Glomeromycota",
    phylum == "Kickxellomycotina" ~ "Kickxellomycota",
    phylum == "Mortierellomycotina" ~ "Mortierellomycota",
    phylum == "Mucoromycotina" ~ "Mucoromycota",
    phylum == "Zoopagomycotina" ~ "Zoopagomycota",
    TRUE ~ phylum
  )) %>% 
  group_by(sample) %>% 
  mutate(percentage = value/sum(value)*100) %>% 
  group_by(sample, phylum) %>% 
  summarise(percentage = sum(percentage)) %>% 
  filter(percentage > 0)

phyla_18$sample <- factor(
  phyla_18$sample,
  levels = c('MP1', 'MP2', 'MP3', 'MP4', 'MP5', 'MP6', 'MP7', 'MP8', 'MP9',
             'MP10', 'MP12', 'MP13', 'MP14', 'MP15', 'MP16', 'MP17', 'MP18',
             'MP19', 'MP20', 'MP21', 'MP23', 'MP24', 'MP25', 'MP26', 'MP27',
             'MP28', 'MP29', 'MP30', 'MP31', 'MP32', 'MP33', 'MP34', 'MP35',
             'MP36') 
)

b1 <- c('MP13', 'MP14', 'MP15', 'MP16', 'MP17', 'MP18', 'MP19', 'MP20', 'MP21',
        'MP23', 'MP24')
b2 <- c('MP25', 'MP26', 'MP27', 'MP28', 'MP29', 'MP30', 'MP31', 'MP32', 'MP33',
        'MP34', 'MP35', 'MP36')
b5 <- c('MP1', 'MP2', 'MP3', 'MP4', 'MP5', 'MP6', 'MP7', 'MP8', 'MP9', 'MP10',
        'MP12')

phyla_18 <- phyla_18 %>% 
  mutate(group = case_when(
    sample %in% b1 ~ 'B1',
    sample %in% b2 ~ 'B2',
    sample %in% b5 ~ 'B5',
    TRUE ~ NA_character_
  ))

phyla_18 %>%
  ggplot(aes(x=sample, y=reorder(phylum, percentage), size=percentage, color=phylum)) +
  geom_point(alpha=0.80) +
  scale_size(range=c(1,5), name='Relative Abundance (%)', breaks=c(10,30,50)) +
  scale_color_manual(values = c(
                               '#9F1D1Dff',
                               "#D07118ff",
                               
                               '#8f50ab',
                               "#DE85D1",
                               '#98823c',

                               '#99c161',
                               "#A9AC48",

                               '#9a5ea1')) +
  theme_bw() +
  theme(
    legend.position = ('bottom'),
    legend.text = element_text(size=10),
    legend.title = element_text(size=10),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=10, angle=90, vjust=0.5, hjust=0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=12),
    strip.text = element_text(size=13),
    axis.ticks.length=unit(.25,'cm')
  ) +
  guides(color='none') +
  facet_grid(.~group, scales='free_x', space='free_x', switch='y')



## Proceeding to order plot
order <- sh_tax %>% 
  filter(sh %in% sh_list) %>% 
  filter(!order=='') %>% 
  select(order, sh) %>% 
  distinct()
order <- left_join(ddOTUtable, order, by=c('name'='sh'))
order <- order %>%
  na.exclude() %>%
  filter(!grepl('_',order)) %>% 
  group_by(sample) %>% 
  mutate(percentage = value/sum(value)*100) %>% 
  group_by(sample, order) %>% 
  summarise(percentage = sum(percentage)) %>%
  mutate(order = ifelse(percentage <= 2, "Others (<2%)", order)) %>% 
  group_by(sample, order) %>% 
  summarise(percentage=sum(percentage)) %>%
  group_by(order) %>%
  mutate(sample_count = n()) %>%
  filter(sample_count >= 3) %>%
  select(-sample_count) %>%
  ungroup()

allOrders <- unique(order$order)
allOrders <- as.data.frame(allOrders)
colnames(allOrders) <- 'scientificName'
write.csv(allOrders, 'data/orders.csv')
## Passed the list on the GBIF species matching tool (https://www.gbif.org/tools/species-lookup)
allOrders <- read.csv('./data/18s/orders_GBIF.csv')
allOrders[7,10] <- 'Zzzz'
allOrders[24,10] <- 'Ascomycota'
orderedOrders <- allOrders$verbatimScientificName[order(allOrders$phylum, allOrders$verbatimScientificName)]


order$order <- factor(
  order$order,
  levels = orderedOrders
)
order$sample <- factor(
  order$sample,
  levels = c('MP1', 'MP2', 'MP3', 'MP4', 'MP5', 'MP6', 'MP7', 'MP8', 'MP9',
             'MP10', 'MP12', 'MP13', 'MP14', 'MP15', 'MP16', 'MP17', 'MP18',
             'MP19', 'MP20', 'MP21', 'MP23', 'MP24', 'MP25', 'MP26', 'MP27',
             'MP28', 'MP29', 'MP30', 'MP31', 'MP32', 'MP33', 'MP34', 'MP35',
             'MP36') 
)

b1 <- c('MP13', 'MP14', 'MP15', 'MP16', 'MP17', 'MP18', 'MP19', 'MP20', 'MP21',
        'MP23', 'MP24')
b2 <- c('MP25', 'MP26', 'MP27', 'MP28', 'MP29', 'MP30', 'MP31', 'MP32', 'MP33',
        'MP34', 'MP35', 'MP36')
b5 <- c('MP1', 'MP2', 'MP3', 'MP4', 'MP5', 'MP6', 'MP7', 'MP8', 'MP9', 'MP10',
        'MP12')

order <- order %>% 
  mutate(group = case_when(
    sample %in% b1 ~ 'B1',
    sample %in% b2 ~ 'B2',
    sample %in% b5 ~ 'B5',
    TRUE ~ NA_character_
  ))


ggplot(order, aes(fill=order, y=percentage, x=forcats::fct_rev(sample))) +
  geom_bar(position='stack', stat='identity', width=.95, show.legend = TRUE) +
  scale_fill_manual(values = c('#F38080ff',
                               '#F07C7Cff',
                               '#EC7878ff',
                               '#E97474ff',
                               '#E57070ff',
                               '#E26B6Bff',
                               '#DE6767ff',
                               '#DB6363ff',
                               '#D75F5Fff',
                               '#D45B5Bff',
                               '#D05757ff',
                               '#CD5353ff',
                               '#C94F4Fff',
                               '#C64A4Aff',
                               '#C24646ff',
                               '#BF4242ff',
                               '#BB3E3Eff',
                               '#B83A3Aff',
                               '#B43636ff',
                               '#B13232ff',
                               '#AD2E2Eff',
                               '#AA2929ff',
                               '#A62525ff',
                               '#A32121ff',
                               '#9F1D1Dff',
                               
                               "#F2B780ff",
                               "#F0B279ff",
                               "#EDAD71ff",
                               "#EBA86Aff",
                               "#E8A362ff",
                               "#E69E5Bff",
                               "#E39953ff",
                               "#E1944Cff",
                               "#DF8F45ff",
                               "#DC8A3Dff",
                               "#DA8536ff",
                               "#D7802Eff",
                               "#D57B27ff",
                               "#D2761Fff",
                               "#D07118ff",
                               
                               "#94DFE7",
                               
                               "#A9E485",
                               '#8AC864',
                               "#6AAB42",
                               
                               "#DE85D1",
                               
                               "#E2E681",
                               "#C6C965",
                               "#A9AC48",
                               
                               "#9aabba")) +
  #  theme_bw() +
  theme(
    legend.position = ('none'),
    legend.title = element_blank(),
    legend.text = element_text(size=10),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=20, face='bold', angle=90, vjust=0.5, hjust=0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=40),
    strip.text = element_text(size=13),
    panel.background = element_blank(),
    axis.ticks.length=unit(.25,'cm')
  ) +
  scale_y_continuous(limits = c(0, 100.05), expand = c(0, 0.5)) +
  ylab('Relative abundance (%)') +
  facet_grid(.~group, scales='free_x', space='free_x', switch='y')


ggsave(
  '18S_StackedBarplot_orderPhylum_abund.pdf', 
  dpi = 500, width=10000, height=5400,
  units='px', plot=last_plot())

## Bubble plot

order %>%
  ggplot(aes(x=sample, y=order, size=percentage, color=order)) +
  geom_point(alpha=0.80) +
  scale_size(range=c(1,10), name='Relative Abundance (%)', breaks=c(10,20,30,40,50)) +
  scale_color_manual(values=c(
                              "#D55D5D",
                              "#D05757",
                              "#CB5151",
                              "#C74C4C",
                              "#C24646",
                              "#BD4040",
                              "#B83A3A",
                              "#B33434",
                              "#AE2E2E",
                              "#A92929",
                              "#A42323",
                              "#9F1D1D",
                              "#8b1415",
                              "#780a0d",
                              "#650101",
                              
                              "#E39953",
                              "#DF8F45",
                              "#D07118",
                              "#a75b14",
                              "#8f4e11",
                              "#83470f",
                              "#77410e",
                              
                              "#DE85D1",
                              
                              "#E2E681",
                              "#C6C965",
                              "#A9AC48",
                              
                              "#9aabba"
  )) +
  theme_bw() +
  theme(
    legend.position = ('bottom'),
    legend.text = element_text(size=10),
    legend.title = element_text(size=10),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=10, angle=90, vjust=0.5, hjust=0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=12),
    strip.text = element_text(size=13),
    axis.ticks.length=unit(.25,'cm')
  ) +
  guides(color='none') +
  scale_y_discrete(limits=rev) +
  facet_grid(.~group, scales='free_x', space='free_x', switch='y')

ggsave(
  '18S_Bubble_orderPhylum_abund.pdf', 
  dpi = 500, width=6000, height=3500,
  units='px', plot=last_plot())

## Relative percentage of phyla
phyla_18 <- left_join(order, allOrders[,c('verbatimScientificName', 'phylum')], by=c('order'='verbatimScientificName')) %>% 
  group_by(sample, phylum) %>%
  summarise(perc=sum(percentage)) %>% 
  group_by(phylum) %>% 
  summarise(perc=mean(perc))

#### Genera ####
ddOTUtable <- read.csv('./data/18s/ddOTU_table.csv', header=TRUE)
rownames(ddOTUtable) <- ddOTUtable$X
ddOTUtable$X <- NULL
rownames(ddOTUtable) <- gsub('.blast', '', as.character(rownames(ddOTUtable)))
ddOTUtable <- as.data.frame(ddOTUtable)
ddOTUtable$sample <- rownames(ddOTUtable)
sh_list <- as.list(colnames(ddOTUtable[1:(length(ddOTUtable)-1)]))
ddOTUtable <- ddOTUtable %>% 
  pivot_longer(cols=-sample)


sh_tax <- read.csv('./data/18s/ddOTU_taxonomy.csv')

genera <- sh_tax %>% 
  filter(sh %in% sh_list) %>% 
  filter(!genus=='') %>% 
  select(genus, sh) %>% 
  distinct()
genera <- left_join(ddOTUtable, genera, by=c('name'='sh'))
genera <- genera %>%
  na.exclude() %>%
  filter(!grepl('les',genus)) %>% 
  filter(!grepl('eta',genus)) %>%
  filter(!grepl('etes',genus)) %>%
  filter(!grepl('dae',genus)) %>%
  filter(!genus%in%c('insect','uncultured', 'cf.', 'euascomycete',
                     'metagenome', 'unidentified', 'fungal', ' ', '',
                     'eukaryotic', 'JCM8259', 'Leotiomycetes', 'Leotiomyceta',
                     'Glomeromycotina', 'Mucoromycotina')) %>%  
  group_by(sample) %>% 
  mutate(percentage = value/sum(value)*100) %>%
  group_by(sample, genus) %>% 
  summarise(percentage = sum(percentage)) %>%
  mutate(genus = ifelse(percentage < 2, "Others (<2%)", genus)) %>% 
  group_by(sample, genus) %>% 
  summarise(percentage=sum(percentage)) %>%
  group_by(genus) %>%
  mutate(sample_count = n()) %>%
  filter(sample_count >= 3) %>%
  select(-sample_count) %>%
  ungroup()

allGen <- unique(genera$genus)
allGen <- as.data.frame(allGen)
colnames(allGen) <- 'scientificName'
write.csv(allGen, './data/18s/genera.csv')
## Passed the list on the GBIF species matching tool
allGen <- read.csv('./data/18s/genera_GBIF.csv')
allGen[3,10] <- 'Zzzz'
orderedGen <- allGen[order(allGen$phylum),] %>% select(verbatimScientificName) %>% unlist()


genera$genus <- factor(
  genera$genus,
  levels = orderedGen
)
genera$sample <- factor(
  genera$sample,
  levels = c('MP1', 'MP2', 'MP3', 'MP4', 'MP5', 'MP6', 'MP7', 'MP8', 'MP9',
             'MP10', 'MP12', 'MP13', 'MP14', 'MP15', 'MP16', 'MP17', 'MP18',
             'MP19', 'MP20', 'MP21', 'MP23', 'MP24', 'MP25', 'MP26', 'MP27',
             'MP28', 'MP29', 'MP30', 'MP31', 'MP32', 'MP33', 'MP34', 'MP35',
             'MP36') 
)

b1 <- c('MP13', 'MP14', 'MP15', 'MP16', 'MP17', 'MP18', 'MP19', 'MP20', 'MP21',
        'MP23', 'MP24')
b2 <- c('MP25', 'MP26', 'MP27', 'MP28', 'MP29', 'MP30', 'MP31', 'MP32', 'MP33',
        'MP34', 'MP35', 'MP36')
b5 <- c('MP1', 'MP2', 'MP3', 'MP4', 'MP5', 'MP6', 'MP7', 'MP8', 'MP9', 'MP10',
        'MP12')

genera <- genera %>% 
  mutate(group = case_when(
    sample %in% b1 ~ 'B1',
    sample %in% b2 ~ 'B2',
    sample %in% b5 ~ 'B5',
    TRUE ~ NA_character_
  ))

#### stacked barplot ####
ggplot(genera, aes(fill=genus, y=percentage, x=forcats::fct_rev(sample))) +
  geom_bar(position='stack', stat='identity', width=.95, show.legend = TRUE) +
  scale_fill_manual(values = c('#f38181',
                               '#ee7b7a',
                               '#e97574',
                               '#e46f6e',
                               '#df6967',
                               '#da6361',
                               '#d55d5a',
                               '#cf5754',
                               '#ca514e',
                               '#c54b48',
                               '#c04541',
                               '#ba3f3b',
                               '#b53935',
                               '#b0322f',
                               '#aa2c29',
                               '#a52523',
                               '#9f1d1d',
                               
                               '#f2b882',
                               '#efb27a',
                               '#edad72',
                               '#eaa76a',
                               '#e7a262',
                               '#e59c5b',
                               '#e29753',
                               '#df914b',
                               '#dc8c43',
                               '#d9863b',
                               '#d68133',
                               '#d37b2b',
                               '#d07622',
                               '#cd7018',
                               "#6AAB42",
                               "#C6C965",
                               '#A9AC48',                   
                               "#9aabba")) +
  #  theme_bw() +
  theme(
    legend.position = ('right'),
    legend.title = element_blank(),
    legend.text = element_text(size=20),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=20, face='bold', angle=90, vjust=0.5, hjust=0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=40),
    strip.text = element_text(size=13),
    panel.background = element_blank(),
    axis.ticks.length=unit(.25,'cm')
  ) +
  scale_y_continuous(limits = c(0, 100.05), expand = c(0, 0.5)) +
  ylab('Relative abundance (%)') +
  facet_grid(.~group, scales='free_x', space='free_x', switch='y')

ggsave(
  '18S_StackedBarplot_orderGenera_abund_legend.pdf', 
  dpi = 500, width=10000, height=5400,
  units='px', plot=last_plot()) ## Not used in the manuscript.

## Bubble plot
genera %>%
  ggplot(aes(x=sample, y=genus, size=percentage, color=genus)) +
  geom_point(alpha=0.80) +
  scale_size(range=c(1,8), name='Relative Abundance (%)', breaks=c(10,20,40,60)) +
  scale_color_manual(values=c(
                              "#C74C4C",
                              "#C24646",
                              "#BD4040",
                              "#B83A3A",
                              "#B33434",
                              "#AE2E2E",
                              "#A92929",
                              "#A42323",
                              "#9F1D1D",
                              "#8b1415",
                              "#780a0d",
                              "#650101",
                              
                              "#DF8F45",
                              "#D07118",
                              "#a75b14",
                              "#8f4e11",
                              "#83470f",
                              "#77410e",

                              "#A9AC48",
                              
                              "#9aabba"
  )) +
  theme_bw() +
  theme(
    legend.position = ('bottom'),
    legend.text = element_text(size=10),
    legend.title = element_text(size=10),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=10, angle=90, vjust=0.5, hjust=0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=12),
    strip.text = element_text(size=13),
    axis.ticks.length=unit(.25,'cm')
  ) +
  guides(color='none') +
  scale_y_discrete(limits=rev) +
  facet_grid(.~group, scales='free_x', space='free_x', switch='y')

ggsave(
  '18s_Bubble_genera.pdf', 
  dpi = 500, width=4000, height=2500,
  units='px', plot=last_plot())

## Relative percentage of phyla
phyla <- left_join(genera, allGen[,c('verbatimScientificName', 'phylum')], by=c('genus'='verbatimScientificName')) %>% 
  group_by(sample, phylum) %>%
  summarise(perc=sum(percentage)) %>% 
  group_by(phylum) %>% 
  summarise(perc=mean(perc))


#### Functional relative abundance plots ####
## Genera
ddOTUtable <- read.csv('./data/18s/ddOTU_table.csv', header=TRUE)
rownames(ddOTUtable) <- ddOTUtable$X
ddOTUtable$X <- NULL
rownames(ddOTUtable) <- gsub('.blast', '', as.character(rownames(ddOTUtable)))
ddOTUtable <- as.data.frame(ddOTUtable)
ddOTUtable$sample <- rownames(ddOTUtable)
sh_list <- as.list(colnames(ddOTUtable[1:(length(ddOTUtable)-1)]))
ddOTUtable <- ddOTUtable %>% 
  pivot_longer(cols=-sample)


sh_tax <- read.csv('./data/18s/ddOTU_taxonomy.csv')

genera <- sh_tax %>% 
  filter(sh %in% sh_list) %>% 
  filter(!genus=='') %>% 
  select(genus, sh) %>% 
  distinct()
genera <- left_join(ddOTUtable, genera, by=c('name'='sh'))
genera <- genera %>%
  na.exclude() %>%
  filter(!grepl('les',genus)) %>% 
  filter(!grepl('eta',genus)) %>%
  filter(!grepl('etes',genus)) %>%
  filter(!grepl('dae',genus)) %>%
  filter(!genus%in%c('insect','uncultured', 'cf.', 'euascomycete',
                     'metagenome', 'unidentified', 'fungal', ' ', '',
                     'eukaryotic', 'JCM8259', 'Leotiomycetes', 'Leotiomyceta'))

genList <- unique(genera$genus)
funguild <- data.frame(gen = character(),
                       trophicMode = character(),
                       guild = character(),
                       confidenceRanking = character(),
                       stringsAsFactors = FALSE)

for(genus in genList) {
  res <- GET(paste('https://mycoportal.org/funguild/services/api/db_return.php?qDB=funguild_db&qField=taxon&qText=', genus, sep=''))
  data <- fromJSON(rawToChar(res$content))
  
  if(length(data) != 0) {
    funguild <- rbind(funguild, data.frame(gen = genus,
                                           trophicMode = data$trophicMode,
                                           guild = data$guild,
                                           confidenceRanking = data$confidenceRanking,
                                           stringsAsFactors = FALSE))
  } else {
    funguild <- rbind(funguild, data.frame(gen = genus,
                                           trophicMode = 'Unknown',
                                           guild = 'Unknown',
                                           confidenceRanking = 'Unknown',
                                           stringsAsFactors = FALSE))
  }
}

row.names(funguild) <- NULL
funguild <- funguild[!duplicated(funguild$gen), ]

names(funguild) <- c('genus', 'trophicMode', 'guild', 'confidenceRanking')

genFunc <- genera %>% left_join(funguild, by=c('genus')) %>% 
  dplyr::select(genus, value, sample, trophicMode, guild, confidenceRanking)
colnames(genFunc)[2] <- 'abundance'

genTM <- genFunc %>%
  filter(confidenceRanking!='Possible') %>% 
  separate_rows(trophicMode, sep='-') %>% 
  group_by(sample, trophicMode) %>% 
  summarise(count=sum(abundance)) %>% 
  group_by(sample) %>% 
  mutate(percentage = (count/sum(count))*100) %>% 
  ungroup()

genG <- genFunc %>%
  filter(confidenceRanking!='Possible') %>% 
  separate_rows(guild, sep='-') %>%
  mutate(guild = gsub("\\|", "", guild)) %>% 
  group_by(sample, guild) %>% 
  summarise(count=sum(abundance)) %>% 
  group_by(sample) %>% 
  mutate(percentage = (count/sum(count))*100) %>%
  filter(percentage > 0) %>% 
  ungroup()

# Checking if samples sums to 1
genTM %>% 
  filter(sample=='MP1') %>% 
  select(percentage) %>% 
  sum()

genG %>% 
  filter(sample=='MP1') %>% 
  select(percentage) %>% 
  sum()

# Adding block column
b1 <- c('MP13', 'MP14', 'MP15', 'MP16', 'MP17', 'MP18', 'MP19', 'MP20', 'MP21',
        'MP23', 'MP24')
b2 <- c('MP25', 'MP26', 'MP27', 'MP28', 'MP29', 'MP30', 'MP31', 'MP32', 'MP33',
        'MP34', 'MP35', 'MP36')
b5 <- c('MP1', 'MP2', 'MP3', 'MP4', 'MP5', 'MP6', 'MP7', 'MP8', 'MP9', 'MP10',
        'MP12')

genTM <- genTM %>% 
  mutate(block = case_when(
    sample %in% b1 ~ 'B1',
    sample %in% b2 ~ 'B2',
    sample %in% b5 ~ 'B5',
    TRUE ~ NA_character_
  )) %>% mutate(across(where(is.character), str_trim))

genG <- genG %>% 
  mutate(block = case_when(
    sample %in% b1 ~ 'B1',
    sample %in% b2 ~ 'B2',
    sample %in% b5 ~ 'B5',
    TRUE ~ NA_character_
  )) %>% mutate(across(where(is.character), str_trim))

genTM$percentage <- genTM$percentage/100
genG$percentage <- genG$percentage/100


genG %>%
  filter(percentage > 0.02) %>%
  filter(!guild == 'Unknown') %>% 
  ggplot(aes(x=sample, y=reorder(guild,percentage), size=percentage, color=guild)) +
  geom_point(alpha=0.80) +
  scale_size(range=c(1,10), name='Relative Abundance (%)') +
  scale_color_manual(values=c('Undefined Saprotroph'='#565e21',
                              'Wood Saprotroph'='#6c7629',
                              'Plant Saprotroph'='#839032',
                              'Plant Pathogen'='#742a2f',
                              'Lichenized'='#742a6d',
                              'Endophyte'='#84327c',
                              'Animal Pathogen'='#84373c',
                              'Algal Parasite'='#954549',
                              'Animal Parasite'='#a65257',
                              'Ectomycorrhizal'='#933a8b',
                              'Fungal Parasite'='#b86065',
                              'Animal Symbiotroph'='#a4429a',
                              'Orchid Mycorrhizal'='#b44aaa',
                              'Ericoid Mycorrhizal'='#c552ba',
                              'Undefined Symbiotroph'='#d65bca',
                              'Epiphyte'='#e763da',
                              'Pollen Saprotroph'='#9baa3b',
                              'Bryophyte Parasite'='#c96e74',
                              'Algal Saprotroph'='#b3c444',
                              'Dung Saprotroph'='#cce04d',
                              'Arbuscular Mycorrhizal'='#e763da',
                              'Lichen Parasite'='#c96e60')) +
  theme_bw() +
  theme(
    legend.position = ('bottom'),
    legend.text = element_text(size=10),
    legend.title = element_text(size=10),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=10, angle=90, vjust=0.5, hjust=0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=12),
    strip.text = element_text(size=13),
    axis.ticks.length=unit(.25,'cm')
  ) +
  guides(color='none') +
  facet_grid(.~block, scales='free_x', space='free_x', switch='y')

ggsave('18S_TrophicMode.pdf',
       dpi = 500, width=10000, height=5400,
       units='px', plot=last_plot())

genTM %>%
  group_by(sample, trophicMode) %>% 
  mutate(perc=sum(percentage)) %>%
  filter(perc > 0) %>% 
  distinct(sample, trophicMode, block, perc) %>% 
  ggplot(aes(x=sample, y=reorder(trophicMode,perc), size=perc, color=trophicMode)) +
  geom_point(alpha=0.80) +
  scale_size(range=c(2,10), name='Relative Abundance (%)', breaks=c(0.2, 0.3, 0.5)) +
  scale_color_manual(values=c('Saprotroph'='#565e21',
                              'Pathotroph'='#742a2f',
                              'Symbiotroph'='#742a6d')) +
  theme_bw() +
  theme(
    legend.position = ('bottom'),
    legend.text = element_text(size=10),
    legend.title = element_text(size=10),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=10, angle=90, vjust=0.5, hjust=0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=12),
    strip.text = element_text(size=13),
    axis.ticks.length=unit(.25,'cm')
  ) +
  guides(color='none') +
  facet_grid(.~block, scales='free_x', space='free_x', switch='y')
