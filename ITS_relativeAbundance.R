## Relative abundance plots code for generating the results related with the ITS dataset generated and presented in 
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
ddOTUtable <- read.csv('./data/its/ddOTU_table.csv', header=TRUE)
rownames(ddOTUtable) <- ddOTUtable$X
ddOTUtable$X <- NULL
rownames(ddOTUtable) <- gsub('.blast', '', as.character(rownames(ddOTUtable)))
ddOTUtable <- as.data.frame(ddOTUtable)
ddOTUtable$sample <- rownames(ddOTUtable)
sh_list <- as.list(colnames(ddOTUtable[1:(length(ddOTUtable)-1)]))
ddOTUtable <- ddOTUtable %>% 
  pivot_longer(cols=starts_with('SH'))


sh_tax <- read.csv('./data/its/ddOTU_taxonomy.csv')

## Checking phyla relative abundance numbers
phyla_its <- sh_tax %>% 
  filter(sh %in% sh_list) %>% 
  filter(!phylum=='') %>% 
  select(phylum, sh) %>% 
  distinct()

phyla_its$phylum <- gsub('p__', '', as.character(phyla_its$phylum))
phyla_its <- left_join(ddOTUtable, phyla_its, by=c('name'='sh'))
phyla_its <- phyla_its %>%
  na.exclude() %>%
  group_by(sample) %>% 
  mutate(percentage = value/sum(value)*100) %>% 
  group_by(sample, phylum) %>% 
  summarise(percentage = sum(percentage)) %>%
  group_by(phylum) %>% 
  summarise(perc=mean(percentage)) ## Mean relative abundance in samples

## Phylum plot
phyla_its <- sh_tax %>% 
  filter(sh %in% sh_list) %>% 
  filter(!phylum=='') %>% 
  select(phylum, sh) %>% 
  distinct()

phyla_its$phylum <- gsub('p__', '', as.character(phyla_its$phylum))
phyla_its <- left_join(ddOTUtable, phyla_its, by=c('name'='sh'))
phyla_its <- phyla_its %>%
  na.exclude() %>%
  group_by(sample) %>% 
  mutate(percentage = value/sum(value)*100) %>% 
  group_by(sample, phylum) %>% 
  summarise(percentage = sum(percentage)) %>% 
  filter(percentage > 0)


phyla_its$sample <- factor(
  phyla_its$sample,
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

phyla_its <- phyla_its %>% 
  mutate(group = case_when(
    sample %in% b1 ~ 'B1',
    sample %in% b2 ~ 'B2',
    sample %in% b5 ~ 'B5',
    TRUE ~ NA_character_
  ))

ordering <- c('Ascomycota', 'Basidiomycota', 'Chytridiomycota', 'Glomeromycota',
              'Mucoromycota', 'Mortierellomycota', 'Zoopagomycota', 'Rozellomycota',
              'Neocallimastigomycota', 'Blastocladiomycota', 'Olpidiomycota',
              'Aphelidiomycota', 'Entomophthoromycota', 'Calcarisporiellomycota',
              'Entorrhizomycota', 'Kickxellomycota')
phyla_its$phylum <- factor(
  phyla_its$phylum,
  levels = ordering
)

phyla_its %>%
  ggplot(aes(x=sample, y=reorder(phylum,percentage), size=percentage, color=phylum)) +
  geom_point(alpha=0.80) +
  scale_size(range=c(1,5), name='Relative Abundance (%)', breaks=c(10,30,50)) +
  scale_color_manual(values = c(
    '#9F1D1Dff',
    "#D07118ff",
    "#6AAB42",
    "#DE85D1",
    "#A9AC48",
    '#99c161',
    '#9a5ea1',
    '#4e3348',
    '#84caa8',
    '#9d98c2',
    '#c6763f',
    '#ca568b',
    '#8f50ab',
    '#ff83a2',
    '#870015',
    '#98823c'
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
  facet_grid(.~group, scales='free_x', space='free_x', switch='y')



## Proceeding to order plot 
order <- sh_tax %>% 
  filter(sh %in% sh_list) %>% 
  filter(!order=='') %>% 
  select(order, sh) %>% 
  distinct()
order$order <- gsub('o__', '', as.character(order$order))
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
write.csv(allOrders, './data/its/orders.csv')
## Passed the list on the GBIF species matching tool (https://www.gbif.org/tools/species-lookup)
allOrders <- read.csv('./data/its/orders_GBIF.csv')
allOrders[3,10] <- 'Glomeromycota'
allOrders[30,10] <- 'Zygomycota s.lat.'
allOrders[8,10] <- 'Zzzz'


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
  scale_fill_manual(values = c("#F38080",
                               "#EE7A7A",
                               "#E97474",
                               "#E46F6F",
                               "#DF6969",
                               "#DA6363",
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
                               
                               "#F2B780",
                               "#EDAD71",
                               "#E8A362",
                               "#E39953",
                               "#DF8F45",
                               "#DA8536",
                               "#D57B27",
                               "#D07118",
                               
                               "#94DFE7",
                               
                               "#A9E485",
                               "#9CD978",
                               "#90CD6A",
                               "#83C25D",
                               "#77B64F",
                               "#6AAB42",
                               
                               "#DE85D1",
                               
                               "#E2E681",
                               "#C6C965",
                               "#A9AC48",
                               
                               "#2a5b6f",
                               
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
  'ITS_StackedBarplot_orderPhylum_abund.pdf', 
  dpi = 500, width=10000, height=5400,
  units='px', plot=last_plot())



## Bubble plot
order %>%
  ggplot(aes(x=sample, y=order, size=percentage, color=order)) +
  geom_point(alpha=0.80) +
  scale_size(range=c(1,10), name='Relative Abundance (%)', breaks=c(5,10,15,20)) +
  scale_color_manual(values = c(
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
                               
                               "#E39953",
                               "#DF8F45",
                               "#D07118",
                               
                               "#9CD978",
                               "#90CD6A",
                               "#83C25D",
                               "#77B64F",
                               "#6AAB42",
                               
                               "#DE99D9",
                               "#DE85D1",
                               
                               "#9aabba")) +
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
  'ITS_Bubble_orderPhylum_abund.pdf', 
  dpi = 500, width=6000, height=3500,
  units='px', plot=last_plot())

## Relative percentage of phyla
phyla <- left_join(order, allOrders[,c('verbatimScientificName', 'phylum')], by=c('order'='verbatimScientificName')) %>% 
  group_by(sample, phylum) %>%
  summarise(perc=sum(percentage)) %>% 
  group_by(phylum) %>% 
  summarise(perc=mean(perc))

#### Genera ####
ddOTUtable <- read.csv('./data/its/ddOTU_table.csv', header=TRUE)
rownames(ddOTUtable) <- ddOTUtable$X
ddOTUtable$X <- NULL
rownames(ddOTUtable) <- gsub('.blast', '', as.character(rownames(ddOTUtable)))
ddOTUtable <- as.data.frame(ddOTUtable)
ddOTUtable$sample <- rownames(ddOTUtable)
sh_list <- as.list(colnames(ddOTUtable[1:(length(ddOTUtable)-1)]))
ddOTUtable <- ddOTUtable %>% 
  pivot_longer(cols=starts_with('SH'))


sh_tax <- read.csv('./data/its/ddOTU_taxonomy.csv')

genera <- sh_tax %>% 
  filter(sh %in% sh_list) %>% 
  filter(!genus=='') %>% 
  select(genus, sh) %>% 
  distinct()
genera$genus <- gsub('g__', '', as.character(genera$genus))
genera <- left_join(ddOTUtable, genera, by=c('name'='sh'))
genera <- genera %>%
  na.exclude() %>%
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
write.csv(allGen, './data/its/genera.csv')
## Passed the list on the GBIF species matching tool (https://www.gbif.org/tools/species-lookup)
allGen <- read.csv('./data/its/genera_GBIF.csv')
allGen[6,10] <- 'Zzzz'
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

#### Stacked barplot ####
ggplot(genera, aes(fill=genus, y=percentage, x=forcats::fct_rev(sample))) +
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
                               '#9F1D1Dff',
                               
                               "#F2B780ff",
                               "#e0944f",
                               "#D07118ff",
                               
                               "#A9E485",
                               "#6AAB42",
                               
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
  './StackedBarplot_orderGenera_abund_legend.pdf', 
  dpi = 500, width=10000, height=5400,
  units='px', plot=last_plot()) # Not used in the manuscript
#### Bubble plot ####
genera %>%
  ggplot(aes(x=sample, y=genus, size=percentage, color=genus)) +
  geom_point(alpha=0.80) +
  scale_size(range=c(1,8), name='Relative Abundance (%)', breaks=c(10,20,40,60)) +
  scale_color_manual(values = c(
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
                               '#9F1D1Dff',

                               "#e0944f",
                               "#D07118ff",
                               
                               "#A9E485",
                               "#6AAB42",
                               
                               "#9aabba")) +
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
  './ITS_Bubble_genera.pdf', 
  dpi = 500, width=4000, height=2500,
  units='px', plot=last_plot())

#### Functional relative abundance plots ####
ddOTUtable <- read.csv('./data/its/ddOTU_table.csv', header=TRUE)
rownames(ddOTUtable) <- ddOTUtable$X
ddOTUtable$X <- NULL
rownames(ddOTUtable) <- gsub('.blast', '', as.character(rownames(ddOTUtable)))
ddOTUtable <- as.data.frame(ddOTUtable)
ddOTUtable$sample <- rownames(ddOTUtable)
sh_list <- as.list(colnames(ddOTUtable[1:(length(ddOTUtable)-1)]))
ddOTUtable <- ddOTUtable %>% 
  pivot_longer(cols=starts_with('SH'))


sh_tax <- read.csv('./data/its/ddOTU_taxonomy.csv')

genera <- sh_tax %>% 
  filter(sh %in% sh_list) %>% 
  filter(!genus=='') %>% 
  select(genus, sh) %>% 
  distinct()
genera$genus <- gsub('g__', '', as.character(genera$genus))
genera <- left_join(ddOTUtable, genera, by=c('name'='sh'))
genera <- genera %>% na.exclude()

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

# Checking if samples sums to 1
genTM %>% 
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

genTM$percentage <- genTM$percentage/100
ggplot(genTM, aes(fill=trophicMode, y=percentage, x=forcats::fct_rev(sample))) +
  geom_bar(position='fill', stat='identity') +
  scale_fill_manual(values = c("#b45948", "#755531", "#99c15f", "#9aabba")) +
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
  ylab('Relative read abundance') +
  facet_grid(.~block, scales='free_x', space='free_x', switch='y')


ggsave('./ITStrophicMode.pdf',
       dpi = 500, width=10000, height=5400,
       units='px', plot=last_plot())




