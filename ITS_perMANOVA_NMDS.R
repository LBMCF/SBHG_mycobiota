## perMANOVA and NMDS code for ITS dataset analysis generated and presented in  the manuscript entitled 
## "Unveiling the rich and functionally diverse Southern Brazilian Highland Grasslands soil funga for promoting conservation".

#autor: Kelmer Martins-Cunha
#contact: kelmermartinscunha@gmail.com


require(ggplot2)
require(dplyr)
require(tidyr)
require(vegan)
require(stringr)
require(glue)
require(car)
require(httr)
require(jsonlite)
require(ggpubr)

set.seed(2134567, kind = 'Mersenne-Twister')

#### Loading and formatting data ####
## 1. ddOTU table
ddOTUtable <- read.csv('./data/its/ddOTU_table.csv', header=TRUE)
rownames(ddOTUtable) <- ddOTUtable$X
ddOTUtable$X <- NULL
rownames(ddOTUtable) <- gsub('.blast', '', as.character(rownames(ddOTUtable)))
ddOTUtable <- as.data.frame(ddOTUtable)
ddOTUtable$sample <- rownames(ddOTUtable)

b1 <- c('MP13','MP14','MP15','MP16','MP17','MP18','MP19','MP20','MP21','MP23','MP24')
b2 <- c('MP25', 'MP26','MP27','MP28','MP29','MP30','MP31','MP32','MP33','MP34','MP35','MP36')
b5 <- c('MP1','MP2','MP3','MP4','MP5','MP6','MP7','MP8','MP9','MP10','MP12')
ddOTUtable <- ddOTUtable %>% 
  mutate(block = case_when(
    rownames(ddOTUtable) %in% b1 ~ 'B1',
    rownames(ddOTUtable) %in% b2 ~ 'B2',
    rownames(ddOTUtable) %in% b5 ~ 'B5',
    TRUE ~ NA_character_
  ))
colsums <- colSums(ddOTUtable[,1:(length(ddOTUtable)-2)][-10,])
seleccol <- colsums > 5
ddOTUtable_filtered <- ddOTUtable[,seleccol]

## 2. Edaphic variables (as received by the lab, need some tinkering)
soils <- readxl::read_excel('./data/soils-data.xlsx')
soils[c('block', 'treatment', 'transect', 'plot')] <- str_split_fixed(soils$`Amostras Peld Urubici`,
                                                                      '-', 4)

duplicatedSamples <- soils %>% 
  filter(`Amostras Peld Urubici` %in% soils$`Amostras Peld Urubici`[duplicated(soils$`Amostras Peld Urubici`)]) %>% 
  filter(str_detect(plot, 'q')) %>% 
  group_by(plot) 

duplicatedSamples <- duplicatedSamples %>% 
  mutate_at(colnames(duplicatedSamples[2:15]), mean) %>% 
  unique(.)

soils <- soils %>% 
  group_by(`Amostras Peld Urubici`) %>%
  filter(n()==1) %>% 
  bind_rows(duplicatedSamples) %>% 
  filter(str_detect(plot, 'q')) %>% 
  arrange(block, transect, plot)

soils$X <- c('MP13', 'MP14', 'MP15', 'MP16', 'MP17', 'MP18', 'MP19', 'MP20', 'MP22', 'MP23', 'MP24',
             'MP21', 'MP25', 'MP26', 'MP27', 'MP28', 'MP29', 'MP30', 'MP31', 'MP36', 'MP34', 'MP35',
             'MP33', 'MP1', 'MP2', 'MP3', 'MP4', 'MP5', 'MP6', 'MP7', 'MP8', 'MP10', 'MP11', 'MP12',
             'MP9')

cn <- readxl::read_excel('./data/soils-data2.xls')
cn <- cn %>% 
  select('COT %','MOS %', 'NT %', '...6') %>% 
  filter(!is.na(...6))
cn[c('block', 'treatment', 'transect', 'plot')] <- str_split_fixed(cn$`...6`,
                                                                      '-', 4)
cn$sample <- c('MP1', 'MP2', 'MP3', 'MP4', 'MP5', 'MP6', 'MP7', 'MP8', 'MP9', 'MP10', 'MP11', 'MP12',
               'MP13', 'MP14', 'MP15', 'MP16', 'MP17', 'MP18', 'MP19', 'MP20', 'MP21', 'MP22', 'MP23',
               'MP24', 'MP25', 'MP26', 'MP27', 'MP28', 'MP29', 'MP30', 'MP31', 'MP33', 'MP34', 'MP35',
               'MP36')
cn <- cn %>% 
  select('COT %','MOS %', 'NT %', 'sample') %>% 
  left_join(soils, by=c('sample'='X'))
cn <- cn %>% 
  select('sample','COT %','MOS %','NT %',"% Argila","pH-SMP","P (mg/kg)","K (mg/kg)",
         "Al (cmolc/kg)","Ca (cmolc/kg)","Mg (cmolc/kg)","H + Al (cmolc/kg)","CTC pH7.0  (cmolc/kg)")

metadata_soils <- ddOTUtable %>% 
  dplyr::select(sample, block) %>% 
  left_join(cn, by='sample')
metadata_soils <- metadata_soils[-10,] ##Removing MP32
rownames(metadata_soils) <- metadata_soils$sample
metadata_soils$block <- NULL
metadata_soils$sample <- NULL
str(metadata_soils) ## Just checking column types — all should be numeric
metadata_soils <- metadata_soils %>% 
                    rename('C' = 'COT %',
                           'OM' = 'MOS %',
                           'N' = 'NT %',
                           'Clay' = '% Argila',
                           'pH' = 'pH-SMP',
                           'P' = 'P (mg/kg)',
                           'K' = 'K (mg/kg)',
                           'Al' = 'Al (cmolc/kg)',
                           'Ca' = 'Ca (cmolc/kg)',
                           'Mg' = 'Mg (cmolc/kg)',
                           'H+Al' = 'H + Al (cmolc/kg)',
                           'CEC' = 'CTC pH7.0  (cmolc/kg)')

## PCA to extract first axis for permanova
pca <- rda(metadata_soils, scale=TRUE)
### Extracting scores and loadings to plot
scores <- as.data.frame(scores(pca, display = "sites")[, 1:2])
scores <- scores %>% 
  mutate(block = case_when(
    rownames(scores) %in% b1 ~ 'B1',
    rownames(scores) %in% b2 ~ 'B2',
    rownames(scores) %in% b5 ~ 'B5'
  ))
scores$block <- case_when(
  scores$block == 'B1' ~ 'S1',
  scores$block == 'B2' ~ 'S2',
  scores$block == 'B5' ~ 'S3',
  TRUE ~ scores$block
)

perc_exp <- 100*pca$CA$eig/sum(pca$CA$eig)
perc_exp[1:2] ## Checking explained variation

round_pe <- format(round(perc_exp[1:2], digits=1), nsmall=1)
labs <- c(glue("PC1 ({round_pe[1]}%)"),
          glue("PC2 ({round_pe[2]}%)"))

loadings <- as.data.frame(scores(pca, display = "species"))
loadings$variable <- rownames(loadings)

# Plot biplot
ggplot() +
  # loadings
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "blue") +
  geom_text(data = loadings, aes(x = PC1, y = PC2, label = variable),
            vjust=-0.6, hjust=-0.1, color = "blue", size=5) +
  # scores
  geom_point(data = scores, aes(x = PC1, y = PC2, color=block), size = 5) +
  scale_color_manual(values=c('#549a3e',
                              '#8f8b39',
                              '#d0d09a')) +
  stat_ellipse(data = scores, aes(x = PC1, y = PC2, group = block, color = block), level = 0.95) +
  labs(x=labs[1], y=labs[2], color=NULL) +
  theme_bw() +
  theme(
    legend.text = element_text(size=12),
    axis.title.y = element_text(size=25),
    axis.text.y = element_text(size=13),
    axis.title.x = element_text(size=25),
    axis.text.x = element_text(size=12)
  )
ggsave('./17.01.25/PCA_soil.pdf',
       width=6000, height=5400,
       units='px', plot = last_plot(), dpi = 500)

scores$sample <- rownames(scores)

# Barplot
loadings %>% 
  ggplot(aes(x = reorder(variable,-PC1), y = PC1)) +
  geom_bar(stat='identity', fill='steelblue', color='black') +
  labs(
    x = 'Variable',
    y = 'Contribution to PC1'
  ) +
  theme_bw() +
  theme(
    legend.text = element_text(size=12),
    axis.title.y = element_text(size=25),
    axis.text.y = element_text(size=13),
    axis.title.x = element_text(size=25),
    axis.text.x = element_text(size=12)
  ) +
  coord_flip()
ggsave('./17.01.25/PCA_loadigs_soil.pdf',
       width=3000, height=5400,
       units='px', plot = last_plot(), dpi = 500)


## 3. Vegetational and abiotic variables
veg_data <- read.csv('./data/veg_data.csv', sep=';')
veg_data$mp_samples <- c('MP13','MP14','MP15','MP16','MP17','MP18','MP19','MP20','MP21','MP22','MP23','MP24',
                         'MP25','MP26','MP27','MP28','MP29','MP30','MP31','MP32','MP33','MP34','MP35','MP36',
                         'MP1','MP2','MP3','MP4','MP5','MP6','MP7','MP8','MP9','MP10', 'MP11', 'MP12')
veg_data <- veg_data[-20,] ## Removing MP32.
veg_data <- veg_data %>% 
  select(c('grass_volume','richness','simpson_dom','mp_samples', 'heat_load', 'fire_2022', 'block'))

## Creating final metadata dataframe for permanova
metadata <- veg_data %>% 
  dplyr::select("grass_volume","richness","simpson_dom","mp_samples","heat_load") %>% 
  left_join(scores, by=c('mp_samples'='sample')) %>% 
  filter(!is.na(PC1))


metadata <- ddOTUtable %>% 
  dplyr::select(sample, block) %>% 
  left_join(metadata, by=c('sample'='mp_samples'))
metadata <- metadata[-10,] ## Removing MP32.
metadata <- metadata %>% 
  select(!block.y)
colnames(metadata)[2] <- 'block'
rownames(metadata) <- metadata$sample
metadata$sample <- NULL

str(metadata) ## checking columns classes
metadata$block <- as.factor(metadata$block)

#### Rarefaction ####

# List to store rarefied lists (rarefaction to the sample level)
samples_rarefied_data <- list()
# List to store Hellinger transformed distances
samples_rarefied_data_dist_hell <- list()

for(i in 1:1000){ #Rarefying the data
  
  samples_rarefied_data[[i]] <- samples_rarefied <- rrarefy(ddOTUtable_filtered[1:(length(ddOTUtable_filtered)-2)], 
                                                            sample=min(rowSums(ddOTUtable_filtered[1:(length(ddOTUtable_filtered)-2)])))
  samples_rarefied_data_dist_hell[[i]] <- vegdist(decostand(samples_rarefied_data[[i]][-10,], method='hellinger'), method='bray', na.rm=TRUE) #removing MP32 as we do not have edaphic data for this sample.
  
}

samples_rarefied_data_dist_hell_mean <- plyr::aaply(plyr::laply(samples_rarefied_data_dist_hell, 
                                                                as.matrix), c(2,3), mean)


#### FunGuild ####
## Skip if already executed in another script.

# Creating lists with detected genera within each block
sh_list <- as.list(colnames(ddOTUtable[1:(length(ddOTUtable)-2)]))
sh_tax <- read.csv('./data/its/ddOTU_taxonomy.csv')

genera <- sh_tax %>% 
  filter(sh %in% sh_list) %>% 
  filter(!genus=='') %>% 
  select(genus, sh) %>%
  filter(!str_detect(genus, '-')) %>% 
  distinct()
genera$genus <- gsub('g__', '', as.character(genera$genus))

## Matching genera against FunGuild

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

genFunc <- funguild %>% 
  filter(confidenceRanking != 'Possible') %>% 
  separate_rows(trophicMode)

sh_genera <- inner_join(genera, genFunc, by=c('genus' = 'gen'))
#### Creating set up for constrain permutations only within blocks ####
perm.ctrl <- how(
  within = Within(type="free"),
  plots = Plots(strata=metadata$block, type="none"),
  nperm=9999)

#### NMDS and perMANOVA - rarefied ####
trophic_list <- c('All', 'Saprotroph', 'Symbiotroph', 'Pathotroph')
for(mode in trophic_list){
  
  if (mode=='All'){
    final_matrix <- samples_rarefied_data_dist_hell_mean
  } else {
    sh_filtered <- sh_genera %>% 
      filter(trophicMode==mode)
    rarefied_dist <- list()
    for (i in 1:1000){
      matrix <- as.data.frame(samples_rarefied_data[[i]][-10,])
      matrix <- matrix %>% 
        select(matches(sh_filtered$sh))
      matrix <- as.matrix(matrix)
      matrix <- decostand(matrix, method='hellinger')
      
      rarefied_dist[[i]] <- vegdist(matrix, method='bray')
      
    }
    final_matrix <- rarefied_dist[[1]]
    for (i in 2:1000){
      final_matrix <- final_matrix + rarefied_dist[[i]]
    }
    final_matrix <- final_matrix/1000
  }

  nmds <- metaMDS(final_matrix, k=2, trymax=500)
  data.scores <- as.data.frame(scores(nmds))
  data.scores$block <- ddOTUtable[-10,]$block

  p <- ggplot(data.scores, aes(x=NMDS1, y=NMDS2, color=block)) +
    geom_point(size=8) +
    stat_ellipse(level=0.95) +
    scale_color_manual(values=c('#549a3e',
                                '#8f8b39',
                                '#d0d09a')) +
    theme_bw() +
    theme(
      legend.text = element_text(size=12),
      axis.title.y = element_text(size=25),
      axis.text.y = element_text(size=13),
      axis.title.x = element_text(size=25),
      axis.text.x = element_text(size=12)
    ) +
    guides(colour=guide_legend(title='Block'),
           fill=guide_legend(title='Block'))

  ggsave(glue('./ITS_R_NMDS_{mode}_S={nmds$stress}.pdf'),
         width=6000, height=5400,
         units='px', plot = p, dpi = 500)
  
  perMANOVA <- adonis2(final_matrix~block+PC1+heat_load+grass_volume+simpson_dom, 
                       data=metadata, permutations=perm.ctrl, method='bray')
   write.csv(perMANOVA, glue('./ITS_R_perMANOVA_{mode}.csv'))
  
   bd <- betadisper(as.dist(final_matrix), metadata$block)
   bdboxplot <- data.frame(unlist(metadata$block), unlist(bd$distances))
   names(bdboxplot) <- c('block', 'distance')
   b <- bdboxplot %>%
     ggplot(aes(x=block, y=distance)) +
     geom_boxplot(outlier.shape = NA) +
     theme_bw() +
     geom_boxplot(outlier.shape = NA, aes(colour=block)) +
     scale_colour_manual(values=c('#549a3e',
                                  '#8f8b39',
                                  '#d0d09a')) +
     geom_jitter(size=3.4, alpha=0.3, aes(colour=block)) +
     theme_bw() +
     theme(
       axis.title.x = element_blank(),
       axis.text.x = element_text(size=20, vjust=0.5, hjust=0.5),
       axis.text.y = element_text(size=20),
       axis.title.y = element_text(size=25)) +
     ylab('Distance to centroid')
  
   compare <- list(c('B1', 'B2'), c('B2', 'B5'), c('B1', 'B5'))
   b <- b + stat_compare_means(test='kruskal.test', comparisons=compare) +
     stat_compare_means(label.y=0.4, label.x=0)

   ggsave(glue('./ITS_R_BDISPER_{mode}.pdf'),
          width=6000, height=5400,
          units='px', plot = b, dpi = 500)
  
}

#### NMDS and perMANOVA - non-rarefied ####
trophic_list <- c('All', 'Saprotroph', 'Symbiotroph', 'Pathotroph')
for(mode in trophic_list){
  if (mode=='All'){
    sh_filtered <- as.data.frame(rownames(t(ddOTUtable_filtered[1:(length(ddOTUtable_filtered)-2)])))
    colnames(sh_filtered) <- 'sh'
  } else {
    sh_filtered <- sh_genera %>% 
      filter(trophicMode==mode)
  }
  
  matrix <- ddOTUtable %>% 
    select(matches(sh_filtered$sh))
  hell <- vegdist(decostand(matrix[-10,], method='hellinger'), method='bray', na.rm=TRUE)
  
 nmds <- metaMDS(hell, k=2, trymax=500)
 data.scores <- as.data.frame(scores(nmds))
 data.scores$block <- ddOTUtable[-10,]$block

 p <- ggplot(data.scores, aes(x=NMDS1, y=NMDS2, color=block)) +
   geom_point(size=8) +
   stat_ellipse(level=0.95) +
   scale_color_manual(values=c('#549a3e',
                               '#8f8b39',
                               '#d0d09a')) +
   theme_bw() +
   theme(
     legend.text = element_text(size=12),
     axis.title.y = element_text(size=25),
     axis.text.y = element_text(size=13),
     axis.title.x = element_text(size=25),
     axis.text.x = element_text(size=12)
   ) +
   guides(colour=guide_legend(title='Block'),
          fill=guide_legend(title='Block'))
 ggsave(glue('./ITS_NR_NMDS_{mode}_S={nmds$stress}.pdf'),
        width=6000, height=5400,
        units='px', plot = p, dpi = 500)
 
   
 perMANOVA <- adonis2(hell~block+PC1+heat_load+grass_volume+simpson_dom,
                      data=metadata, permutations=perm.ctrl, method='bray')
 write.csv(perMANOVA, glue('./ITS_NR_perMANOVA_{mode}.csv'))

 bd <- betadisper(as.dist(hell), metadata$block)
 bdboxplot <- data.frame(unlist(metadata$block), unlist(bd$distances))
 names(bdboxplot) <- c('block', 'distance')
 b <- bdboxplot %>%
   ggplot(aes(x=block, y=distance)) +
   geom_boxplot(outlier.shape = NA) +
   theme_bw() +
   geom_boxplot(outlier.shape = NA, aes(colour=block)) +
   scale_colour_manual(values=c('#549a3e',
                                '#8f8b39',
                                '#d0d09a')) +
   geom_jitter(size=3.4, alpha=0.3, aes(colour=block)) +
   theme_bw() +
   theme(
     axis.title.x = element_blank(),
     axis.text.x = element_text(size=20, vjust=0.5, hjust=0.5),
     axis.text.y = element_text(size=20),
     axis.title.y = element_text(size=25)) +
   ylab('Distance to centroid')

 compare <- list(c('B1', 'B2'), c('B2', 'B5'), c('B1', 'B5'))
 b <- b + stat_compare_means(test='kruskal.test', comparisons=compare) +
   stat_compare_means(label.y=0.4, label.x=0)

 ggsave(glue('./ITS_NR_BDISPER_{mode}.pdf'),
        width=6000, height=5400,
        units='px', plot = b, dpi = 500)
  
  
}
