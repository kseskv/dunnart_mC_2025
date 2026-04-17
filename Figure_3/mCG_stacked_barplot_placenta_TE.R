###--- load libraries 
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(wesanderson)

###--- read in data
setwd("/Users/allang/Desktop/dunnart_reprogramming/placenta_te_mCG/")

# squirrel monkey placenta
squirrel_monkey_placenta <- fread("BS_Seeker_squirrelmonkey_placenta_percentmeth.bed", header=F)
input_squirrel_monkey_placenta <- squirrel_monkey_placenta %>%
  select(V4) %>% 
  rename(squirrel_monkey_placenta = V4) %>%
  gather()
input_squirrel_monkey_placenta$value <- input_squirrel_monkey_placenta$value*100
input_squirrel_monkey_placenta_grouped <- input_squirrel_monkey_placenta %>% 
  group_by(key, value=cut(input_squirrel_monkey_placenta$value, breaks=c(0, 0.01, 20, 80, 100), include.lowest = TRUE)) %>% 
  summarise(n=(n()/nrow(input_squirrel_monkey_placenta)*100))
vec <- c("x0", "x0-20", "x20-80", "x80-100")  
input_squirrel_monkey_placenta_grouped$value <- vec
squirrel_monkey_placenta_mean_mCG <- mean(input_squirrel_monkey_placenta$value)

# rhesus trophoblast
rhesus_trophoblast <- fread("BS_Seeker_rhesus_trophoblast_percentmeth.bed", header=F)
input_rhesus_trophoblast <- rhesus_trophoblast %>%
  select(V4) %>% 
  rename(rhesus_trophoblast = V4) %>%
  gather()
input_rhesus_trophoblast$value <- input_rhesus_trophoblast$value*100
input_rhesus_trophoblast_grouped <- input_rhesus_trophoblast %>% 
  group_by(key, value=cut(input_rhesus_trophoblast$value, breaks=c(0, 0.01, 20, 80, 100), include.lowest = TRUE)) %>% 
  summarise(n=(n()/nrow(input_rhesus_trophoblast)*100))
vec <- c("x0", "x0-20", "x20-80", "x80-100")  
input_rhesus_trophoblast_grouped$value <- vec
rhesus_trophoblast_mean_mCG <- mean(input_rhesus_trophoblast$value)

# dog placenta
dog_placenta <- fread("BS_Seeker_dog_plac1_percentmeth.bed", header=F)
input_dog_placenta <- dog_placenta %>%
  select(V4) %>% 
  rename(dog_placenta = V4) %>%
  gather()
input_dog_placenta$value <- input_dog_placenta$value*100
input_dog_placenta_grouped <- input_dog_placenta %>% 
  group_by(key, value=cut(input_dog_placenta$value, breaks=c(0, 0.01, 20, 80, 100), include.lowest = TRUE)) %>% 
  summarise(n=(n()/nrow(input_dog_placenta)*100))
vec <- c("x0", "x0-20", "x20-80", "x80-100")  
input_dog_placenta_grouped$value <- vec
dog_placenta_mean_mCG <- mean(input_dog_placenta$value)

# cow placenta
cow_placenta <- fread("BS_Seeker_cow_placenta_percentmeth.bed", header=F)
input_cow_placenta <- cow_placenta %>%
  select(V4) %>% 
  rename(cow_placenta = V4) %>%
  gather()
input_cow_placenta$value <- input_cow_placenta$value*100
input_cow_placenta_grouped <- input_cow_placenta %>% 
  group_by(key, value=cut(input_cow_placenta$value, breaks=c(0, 0.01, 20, 80, 100), include.lowest = TRUE)) %>% 
  summarise(n=(n()/nrow(input_cow_placenta)*100))
vec <- c("x0", "x0-20", "x20-80", "x80-100")  
input_cow_placenta_grouped$value <- vec
cow_placenta_mean_mCG <- mean(input_cow_placenta$value)

# horse placenta
horse_placenta <- fread("BS_Seeker_horse_placenta_percentmeth.bed", header=F)
input_horse_placenta <- horse_placenta %>%
  select(V4) %>% 
  rename(horse_placenta = V4) %>%
  gather()
input_horse_placenta$value <- input_horse_placenta$value*100
input_horse_placenta_grouped <- input_horse_placenta %>% 
  group_by(key, value=cut(input_horse_placenta$value, breaks=c(0, 0.01, 20, 80, 100), include.lowest = TRUE)) %>% 
  summarise(n=(n()/nrow(input_horse_placenta)*100))
vec <- c("x0", "x0-20", "x20-80", "x80-100")  
input_horse_placenta_grouped$value <- vec
horse_placenta_mean_mCG <- mean(input_horse_placenta$value)

# mouse E3.5 TE
mouse_te <- fread("GSM2229984_GSM2229985_TE_mCG_merged_avg.dec.bed", header=F)
input_mouse_te <- mouse_te %>%
  select(V4) %>% 
  rename(mouse_te = V4) %>%
  gather()
input_mouse_te$value <- input_mouse_te$value*100
input_mouse_te_grouped <- input_mouse_te %>% 
  group_by(key, value=cut(input_mouse_te$value, breaks=c(0, 0.01, 20, 80, 100), include.lowest = TRUE)) %>% 
  summarise(n=(n()/nrow(input_mouse_te)*100))
vec <- c("x0", "x0-20", "x20-80", "x80-100")  
input_mouse_te_grouped$value <- vec
mouse_te_mean_mCG <- mean(input_mouse_te$value)

# opossum E7.5M TE
opossum_te <- fread("GSM6252473_LEE307A178_merged_fastq.gz_trimmed_bismark_bt2.deduplicated.bedGraph", header=F)
input_opossum_te <- opossum_te %>%
  select(V4) %>% 
  rename(opossum_te = V4) %>%
  gather()
input_opossum_te_grouped <- input_opossum_te %>% 
  group_by(key, value=cut(input_opossum_te$value, breaks=c(0, 0.01, 20, 80, 100), include.lowest = TRUE)) %>% 
  summarise(n=(n()/nrow(input_opossum_te)*100))
vec <- c("x0", "x0-20", "x20-80", "x80-100")  
input_opossum_te_grouped$value <- vec
opossum_te_mean_mCG <- mean(input_opossum_te$value)

# dunnart E9 TE (pseudobulk scBSseq)
dunnart_te <- fread("scBSseq_E9_pseudobulk_TE.clean.bismark.cov", header=F)
input_dunnart_te <- dunnart_te %>%
  select(V4) %>% 
  rename(dunnart_te = V4) %>%
  gather()
input_dunnart_te$value <- input_dunnart_te$value*100
input_dunnart_te_grouped <- input_dunnart_te %>% 
  group_by(key, value=cut(input_dunnart_te$value, breaks=c(0, 0.01, 20, 80, 100), include.lowest = TRUE)) %>% 
  summarise(n=(n()/nrow(input_dunnart_te)*100))
vec <- c("x0", "x0-20", "x20-80", "x80-100")  
input_dunnart_te_grouped$value <- vec
dunnart_te_mean_mCG <- mean(input_dunnart_te$value)

# human placenta 
human_placenta <- fread("GSM6455952_human.CpG.0rm.txt", header=F)
input_human_placenta <- human_placenta %>%
  select(V5) %>% 
  rename(human_placenta = V5) %>%
  gather()
input_human_placenta$value <- input_human_placenta$value*100
input_human_placenta_grouped <- input_human_placenta %>% 
  group_by(key, value=cut(input_human_placenta$value, breaks=c(0, 0.01, 20, 80, 100), include.lowest = TRUE)) %>% 
  summarise(n=(n()/nrow(input_human_placenta)*100))
vec <- c("x0", "x0-20", "x20-80", "x80-100")  
input_human_placenta_grouped$value <- vec
human_placenta_mean_mCG <- mean(input_human_placenta$value)

# rat placenta 
rat_placenta <- fread("GSM6455949_B45.CpG.0rm.txt", header=F)
input_rat_placenta <- rat_placenta %>%
  select(V5) %>% 
  rename(rat_placenta = V5) %>%
  gather()
input_rat_placenta$value <- input_rat_placenta$value*100
input_rat_placenta_grouped <- input_rat_placenta %>% 
  group_by(key, value=cut(input_rat_placenta$value, breaks=c(0, 0.01, 20, 80, 100), include.lowest = TRUE)) %>% 
  summarise(n=(n()/nrow(input_rat_placenta)*100))
vec <- c("x0", "x0-20", "x20-80", "x80-100")  
input_rat_placenta_grouped$value <- vec
rat_placenta_mean_mCG <- mean(input_rat_placenta$value)

# opossum EEM
opossum_eem <- fread("BS_Seeker_opossum_EEM_percentmeth.bed", header=F)
input_opossum_eem <- opossum_eem %>%
  select(V4) %>% 
  rename(opossum_eem = V4) %>%
  gather()
input_opossum_eem$value <- input_opossum_eem$value*100
input_opossum_eem_grouped <- input_opossum_eem %>% 
  group_by(key, value=cut(input_opossum_eem$value, breaks=c(0, 0.01, 20, 80, 100), include.lowest = TRUE)) %>% 
  summarise(n=(n()/nrow(input_opossum_eem)*100))
vec <- c("x0", "x0-20", "x20-80", "x80-100")  
input_opossum_eem_grouped$value <- vec
opossum_eem_mean_mCG <- mean(input_opossum_eem$value)

# mouse placenta
mouse_placenta <- fread("BS_Seeker_mouse_placenta_percentmeth.bed", header=F)
input_mouse_placenta <- mouse_placenta %>%
  select(V4) %>% 
  rename(mouse_placenta = V4) %>%
  gather()
input_mouse_placenta$value <- input_mouse_placenta$value*100
input_mouse_placenta_grouped <- input_mouse_placenta %>% 
  group_by(key, value=cut(input_mouse_placenta$value, breaks=c(0, 0.01, 20, 80, 100), include.lowest = TRUE)) %>% 
  summarise(n=(n()/nrow(input_mouse_placenta)*100))
vec <- c("x0", "x0-20", "x20-80", "x80-100")  
input_mouse_placenta_grouped$value <- vec
mouse_placenta_mean_mCG <- mean(input_mouse_placenta$value)

###--- merge datasets
mCG <- rbind(input_squirrel_monkey_placenta_grouped, input_rhesus_trophoblast_grouped, 
             input_dog_placenta_grouped, input_cow_placenta_grouped, 
             input_horse_placenta_grouped, input_mouse_te_grouped, 
             input_opossum_te_grouped, input_dunnart_te_grouped,
             input_human_placenta_grouped, input_rat_placenta_grouped,
             input_opossum_eem_grouped, input_mouse_placenta_grouped)

fill.order <- factor(mCG$key, levels = c(
  'squirrel_monkey_placenta',  # 0.6553342
  'human_placenta',            # 0.6266247
  'horse_placenta',            # 0.5987895
  'dog_placenta',              # 0.5674493
  'opossum_eem',               # 0.5586966
  'mouse_placenta',            # 0.5355082
  'rat_placenta',              # 0.4682574
  'opossum_te',                # 0.4654497
  'rhesus_trophoblast',        # 0.4616784
  'cow_placenta',              # 0.3057065
  'dunnart_te',                # 0.2324198
  'mouse_te'                   # 0.1917322
))

mCG_means <- data.frame(
  sample = c(
    'squirrel_monkey_placenta',  # 0.6553342
    'human_placenta',            # 0.6266247
    'horse_placenta',            # 0.5987895
    'dog_placenta',              # 0.5674493
    'opossum_eem',               # 0.5586966
    'mouse_placenta',            # 0.5355082
    'rat_placenta',              # 0.4682574
    'opossum_te',                # 0.4654497
    'rhesus_trophoblast',        # 0.4616784
    'cow_placenta',              # 0.3057065
    'dunnart_te',                # 0.2324198
    'mouse_te'                   # 0.1917322
  ),
  mean_mCG = c(
    0.6553342,
    0.6266247,
    0.5987895,
    0.5674493,
    0.5586966,
    0.5355082,
    0.4682574,
    0.4654497,
    0.4616784,
    0.3057065,
    0.2324198,
    0.1917322
  )
)


###--- set colour palette
palette <- wes_palette("Zissou1", 4, type = "continuous")

###--- plot table as stacked barplot, save as PDF 
pdf("evo_placenta_TE_mCG_stacked_barplot.pdf", width=6, height=3)
ggplot(mCG, aes(x = fill.order, y = n, fill = factor(value))) + 
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = palette) + 
  xlab("") +
  ylab("CpG sites (%)") +
  coord_flip() +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  # overlay means as dots (adjust x/y accordingly due to coord_flip)
  geom_point(data = mCG_means, 
             aes(x = sample, y = mean_mCG), 
             inherit.aes = FALSE,
             color = "black", size = 3)
dev.off()
