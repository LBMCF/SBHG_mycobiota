## Alpha diversity metrics code for 18S dataset analysis generated and presented in  the manuscript entitled 
## "Unveiling the rich and functionally diverse Southern Brazilian Highland Grasslands soil funga for promoting conservation".

#autor: Kelmer Martins-Cunha
#contact: kelmermartinscunha@gmail.com


require(iNEXT)
require(dplyr)
require(tidyverse)
require(glue)

#### Loading and formatting data ####
ddOTUtable <- read.csv('./data/18s/ddOTU_table.csv', header=TRUE)
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

blocks <- ddOTUtable %>%
  group_by(block) %>%
  summarise(across(where(is.numeric), sum)) %>%
  as.data.frame()
rownames(blocks) <- blocks$block
blocks$block <- NULL





#### FunGuild ####
# Creating lists with detected genera within each block
sh_list <- as.list(colnames(ddOTUtable[1:(length(ddOTUtable)-2)]))
sh_tax <- read.csv('./data/18s/ddOTU_taxonomy.csv')

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


#### Running iNEXT ####
b1 <- c('MP13','MP14','MP15','MP16','MP17','MP18','MP19','MP20','MP21','MP23','MP24')
b2 <- c('MP25', 'MP26','MP27','MP28','MP29','MP30','MP31','MP32','MP33','MP34','MP35','MP36')
b5 <- c('MP1','MP2','MP3','MP4','MP5','MP6','MP7','MP8','MP9','MP10','MP12')
trophic_list <- c('All', 'Saprotroph', 'Symbiotroph', 'Pathotroph')

for(mode in trophic_list){
  if (mode=='All'){
    sh_filtered <- as.data.frame(rownames(t(ddOTUtable[1:(length(ddOTUtable)-2)])))
    colnames(sh_filtered) <- 'sh'
  } else {
    sh_filtered <- sh_genera %>% 
      filter(trophicMode==mode)
  }
  
  matrix <- ddOTUtable %>% 
    select(matches(sh_filtered$sh))
  
  matrix <- matrix %>% 
    mutate(block = case_when(
      rownames(matrix) %in% b1 ~ 'B1',
      rownames(matrix) %in% b2 ~ 'B2',
      rownames(matrix) %in% b5 ~ 'B5',
      TRUE ~ NA_character_
    ))
  
  blocks <- matrix %>%
    group_by(block) %>%
    summarise(across(where(is.numeric), sum)) %>%
    as.data.frame()
  rownames(blocks) <- blocks$block
  blocks$block <- NULL
  
  re_uncons <- iNEXT(t(blocks), q=c(0,1,2), datatype='abundance', se=TRUE)
  asypmtotic_all <- re_uncons$AsyEst
  
  re_cons <- iNEXT(t(blocks), q=c(0,1,2), datatype='abundance', endpoint=(min(rowSums(blocks)))*2)
  estimates_cons <- re_cons$iNextEst$coverage_based
  estimates_cons$se <- (estimates_cons$qD.UCL-estimates_cons$qD.LCL)/(2*1.96)
  
  ggiNEXT(re_cons, type=1, se=TRUE, facet.var = 'Order.q') +
    scale_color_manual(values=c('#549a3e',
                                '#d0d09a',
                                '#8f8b39')) +
    scale_fill_manual(values=c('#549a3e',
                               '#d0d09a',
                               '#8f8b39')) +
    facet_wrap(~Order.q, scales='free', strip.position='left',
               labeller=as_labeller(c('0'='Effective ddOTU richness', '1'='Effective Shannon diversity',
                                      '2'='Effective Simpson diversity'))) +
    theme_bw() +
    theme(
      legend.text = element_text(size=12),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=20),
      axis.title.x = element_text(size=20),
      axis.text.x = element_text(size=20),
      strip.background = element_blank(),
      strip.placement = 'outside',
      strip.text.y.left = element_text(size=20)
    ) +
    guides(shape='none') +
    labs(x='Sequences')
  
  ggsave(glue('./18S_q_orders_{mode}.pdf'), 
         width=10000, height=3000,
         units='px', plot = last_plot(), dpi = 500)
  
}
