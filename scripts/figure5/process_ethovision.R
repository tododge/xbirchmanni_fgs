# load packages
detach("package:tidyverse", unload=TRUE)
detach("package:Rmisc", unload=TRUE)
library(Rmisc)
library(tidyverse)
library(dplyr)
library(readxl)
library(gridExtra)

# figuring out how to process these files in a reproducible way!

# EthoVision column names
cols <- c("trial_time", "recording_time", "x", "y", "area", "area_change", "elongation", "distance", "velocity", "zone1_25", "zone2_25", "zone2_5", "zone1_5", "result1")
##### NOTE: some of the other trials will have multiple "zone" names, so will have more columns
##### NOTE: will need to add an extra step later to combine them into a single column!


# make tibble to add data to
summary_tibble <- tibble(total_time=numeric(), trial_time=numeric(), 
                         z1_total=numeric(), z2_total=numeric(), 
                         z1_latency=numeric(), z2_latency=numeric(), 
                         z1_trial=numeric(), z2_trial=numeric(),
                         distance=numeric(), speed=numeric(), 
                         file_name=character(), trial_lane=character())

# list out all data file names to process
# these are the raw data output files from EthoVision
file_names <- Sys.glob("/Users/trisdodge/Desktop/Stanford/falsegravid/behavior/female_pref_data2/Raw*.xlsx")

# frequency of tracking data collected in EthoVision
scaling_factor <- 0.089 #in ethovision


# loop to process each file and append data
for (filename in file_names) {
  #filename <- "/Users/trisdodge/Desktop/Stanford/falsegravid/behavior/female_pref_data/Raw data-Xbir_FGS_female_preference_tracking-Trial     1.xlsx"
  
  video_name <- read_excel(filename, range = "B19", col_names = FALSE)[[1]]
  video_name <- gsub("F:\\\\xbir_FGS_female_preference_raw_videos\\\\Media Files\\\\trimmed_videos\\\\", "", video_name)
  arena <- read_excel(filename, range = "B10", col_names = FALSE)[[1]]
  trial_lane <- substr(arena, 1, 1)
  
  trial_data <- read_excel(filename, skip=36, col_names=cols, na="NA") %>%
##### NOTE: can add a step here to combine columns into `zone1` and `zone2` only
    dplyr::select(trial_time, x, y, zone1_25, zone2_25, zone1_5, zone2_5) %>%
    dplyr::filter(trial_time <= 600)

  #calculate distance
  trial_data$distance <- NA
  
  trial_data <- trial_data %>% 
  mutate(x = as.double(x)) %>% 
  mutate(y = as.double(y)) %>% 
  mutate(zone1_25 = as.double(zone1_25)) %>% 
  mutate(zone2_25 = as.double(zone2_25)) %>% 
  mutate(zone1_5 = as.double(zone1_5)) %>% 
  mutate(zone2_5 = as.double(zone2_5))
  
  for (i in 1:nrow(trial_data)) {
    trial_data$distance[i] <- ifelse(i==1, 0, sqrt(abs(trial_data$x[i] - trial_data$x[i-1])^2 + abs(trial_data$y[i] - trial_data$y[i-1])^2))
  }
  
  #calculate speed
  trial_data <- trial_data %>%
    mutate(speed = distance/scaling_factor)  

    
#calculate summary stats for the file
  
  #total time spent in each zone
  #validated with EthoVision summary output
  z1_total_time <- sum(trial_data$zone1_25 + trial_data$zone1_5, na.rm=TRUE)*scaling_factor
  z2_total_time <- sum(trial_data$zone2_25 + trial_data$zone2_5, na.rm=TRUE)*scaling_factor

  #latency to enter each zone
  #prints NA if fish never entered the zone
  #validated with EthoVision summary output
  z1_latency <- ifelse(dim(filter(trial_data, zone1_25==1))[1] == 0, 
                       NA, 
                       (min(dplyr::select(filter(trial_data, zone1_25==1), trial_time))))
  z2_latency <- ifelse(dim(dplyr::filter(trial_data, zone2_25==1))[1] == 0, 
                       NA, 
                       (min(dplyr::select(filter(trial_data, zone2_25==1), trial_time))))
  max_latency <- max(z1_latency, z2_latency)
  
  #cut down trial df to include times AFTER the latest latency time only
  
  ####SHOULD WE CHANGE THIS TO ONLY INCLUDE TIME AFTER LEAVING THIS SECOND LATENCY ZONE????####
  
  trim_trial_data <- trial_data %>%
    filter(trial_time >= max_latency)
  #calculate total time in zones AFTER the latest latency time
  z1_post_latency <- sum(trim_trial_data$zone1_25 + trim_trial_data$zone1_5, na.rm=TRUE)*scaling_factor
  z2_post_latency <- sum(trim_trial_data$zone2_25 + trim_trial_data$zone2_5, na.rm=TRUE)*scaling_factor

  #total time of the trial (max=300, but could be slightly less depending on video editing)
  total_length <- max(trial_data$trial_time)
  #calculate length of trial after the latest latency time
  trial_total = total_length - max_latency
  
  #calculate total distance traveled (cm) AFTER the latest latency time
  total_distance <- sum(trim_trial_data$distance, na.rm=TRUE)
  #calculate average speed (cm/s) AFTER the latest latency time
  avg_speed <- mean(trim_trial_data$speed, na.rm=TRUE)
  
  #append the  summary stats  from the individual file
  summary_tibble <- add_row(summary_tibble, total_time=total_length, trial_time=trial_total, z1_total=z1_total_time, z2_total=z2_total_time, 
                            z1_latency=z1_latency, z2_latency=z2_latency, 
                            z1_trial=z1_post_latency, z2_trial=z2_post_latency,
                            distance=total_distance, speed=avg_speed, 
                            file_name=video_name, trial_lane=trial_lane)
}



# load trial keyfile
design <- read_excel("/Users/trisdodge/Desktop/Stanford/falsegravid/behavior/female_pref_data_collection.xlsx")
animations2zones <- read_excel("/Users/trisdodge/Desktop/Stanford/falsegravid/behavior/female_pref_animations_2_zones.xlsx")
keyfile <- read_excel("/Users/trisdodge/Desktop/Stanford/falsegravid/behavior/female_pref_keyfile.xlsx")

# combine the summary tibble and the keyfile
#summary_tibble2 <- summary_tibble %>% drop_na()
trial_summary <- left_join(summary_tibble, keyfile, by=c("file_name","trial_lane"))
trial_summary2 <- left_join(trial_summary, animations2zones, by="animation_segment_lane")
trial_summary3 <- left_join(trial_summary2, design, by=c('trial_num'='trial', 'trial_lane.x'='trial_lane...3'))

####TESTING####

FGSvsWT <- subset(trial_summary3, comparision != "BvJ")

#FGSvsWT <- subset(trial_summary3, comparision != "BvJ" & origin == "WC")

#FGSvsWT <- subset(trial_summary3, comparision == "JvJ" & origin == "WC")
#FGSvsWT <- subset(trial_summary3, comparision == "BvB" & origin == "WC")

FGSvsWT$FGS_zone <- ifelse(FGSvsWT$FGS_zone_1=="FGS", FGSvsWT$z1_trial, FGSvsWT$z2_trial)
FGSvsWT$WT_zone <- ifelse(FGSvsWT$FGS_zone_1=="non-FGS", FGSvsWT$z1_trial, FGSvsWT$z2_trial)

FGSvsWT$FGS_zone_prop <- FGSvsWT$FGS_zone / FGSvsWT$trial_time
#FGSvsWT$FGS_zone_sec <- FGSvsWT$FGS_zone

FGSvsWT$WT_zone_prop <- FGSvsWT$WT_zone / FGSvsWT$trial_time
#FGSvsWT$WT_zone_sec <- FGSvsWT$WT_zone

FGSvsWT$prop_diff <- FGSvsWT$FGS_zone_prop - FGSvsWT$WT_zone_prop
#FGSvsWT$sec_diff <- FGSvsWT$WT_zone_sec - FGSvsWT$FGS_zone_sec
FGSvsWT$trial_length <- FGSvsWT$trial_time

FGSvsWT$diff <- (FGSvsWT$FGS_zone - FGSvsWT$WT_zone) / (FGSvsWT$WT_zone + FGSvsWT$FGS_zone)



FGSvsWT_mean <- FGSvsWT %>% group_by(trial_lane...4) %>%
  dplyr::summarise(across(FGS_zone_prop:trial_length, ~ mean(.x, na.rm = TRUE)))

FGSvsWT_mean2 <- FGSvsWT %>% group_by(trial_lane...4, comparision, origin) %>%
  dplyr::summarise(across(FGS_zone_prop:trial_length, ~ mean(.x, na.rm = TRUE)))

names(FGSvsWT)
FGSvsWT_metadata <- FGSvsWT[,c("trial_lane...4","trial_lane.x","trial_num","animation.x","exclude", "recorder.x","date","trial_day","start_time","origin","pigmentation")]
FGSvsWT_metadata <- unique(FGSvsWT_metadata)

data <- merge(FGSvsWT_mean2,FGSvsWT_metadata)
data$origin <- gsub("LB","lab-reared females",data$origin)
data$origin <- gsub("WC","wild-caught females",data$origin)
data$comparision <- gsub("JvJ","unornamented\nFGS vs non-FGS",data$comparision)
data$comparision <- gsub("BvB","ornamented\nFGS vs non-FGS",data$comparision)


data$comparision <- fct_relevel(data$comparision, "unornamented\nFGS vs non-FGS", "ornamented\nFGS vs non-FGS")
data$origin <- fct_relevel(data$origin, "wild-caught females", "lab-reared females")
data_long <- data[,c("trial_lane...4","origin","comparision","FGS_zone_prop","WT_zone_prop", "prop_diff")] %>% 
  pivot_longer(
    cols = FGS_zone_prop:WT_zone_prop, 
    names_to = "zone",
    values_to = "prop"
  )
data_long$zone <- gsub("FGS_zone_prop", "FGS", data_long$zone)
data_long$zone <- gsub("WT_zone_prop", "non-FGS", data_long$zone)

names(data_long)
data_long_summary <- summarySE(na.omit(data_long[,c("trial_lane...4","origin","comparision","zone","prop")]), measurevar="prop", groupvars=c("origin","comparision","zone"))


fig6e <- 
  ggplot() +
    geom_point(data = data_long, aes(x = comparision, y = prop, group = zone, color = zone), position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.1), alpha=0.2) +
    geom_errorbar(data = data_long_summary, aes(x = comparision, ymin = prop - 2 * se, ymax = prop + 2 * se, group = zone, color=zone), width = 0.2, position = position_dodge(width = 0.5)) +
    geom_point(data = data_long_summary, aes(x = comparision, y = prop, group = zone, color = zone), position = position_dodge(width = 0.5), size = 3) + 
    facet_wrap(~origin, nrow = 1, strip.position = "bottom") +
    xlab("") +
    theme_classic(base_size = 7) +
    scale_color_manual(values = c("#55185D", "darkgrey")) +
    theme(legend.position="none",
      panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          strip.background = element_blank(), 
          strip.placement = "outside",
        strip.text.x = element_text(size=7)) +
    ylab("Female association (trial prop.)") #+
  #stat_compare_means(data=data_long, aes(x=comparision, y=prop, group = zone),label = "p.format", paired=TRUE)


fig6e_origional <- 
  ggplot() +
  geom_jitter(data=subset(data), aes(x=comparision,y=prop_diff, group=comparision, fill=prop_diff), pch=21, width=0.1, color="grey") +
    geom_errorbar(data=data_summary, aes(x=comparision, ymin=prop_diff-2*se, ymax=prop_diff+2*se, group=comparision), width=0.1) +
    geom_point(data=data_summary, aes(x=comparision, y=prop_diff, group=comparision),  color="black", size=3) + 
  facet_wrap(~origin, nrow = 1, strip.position = "bottom") +
  xlab("") +
    #labs(fill="FGS pref.") +
  scale_y_continuous(limits=c(-0.75,0.75),breaks = seq(-0.75, 0.75, 0.25)) +
  theme_classic(base_size = 8) +
  #scale_color_manual(values = c("#ECB602","#8538a8")) +
  theme(#legend.position="none",
    panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_),
        strip.background = element_blank(), 
    strip.placement = "outside",
    legend.key.height=unit(7, 'mm'),
    legend.key.width=unit(2, 'mm'),
    legend.title = element_blank(),
    legend.text=element_blank()
    ) +
  ylab("female FGS pref. index") +
  geom_hline(yintercept = 0, lty="dashed", alpha = 0.5) +
    scale_fill_gradient2(
      low = "darkgrey", 
      mid = "white", 
      high = "#55185D", 
      midpoint = 0
    )

ggplot(data=data, aes(x=origin,y=prop_diff)) + 
  geom_violin(aes(col=comparision)) + 
  geom_point(aes(col=comparision), position = position_jitterdodge()) + 
  stat_summary(fun = "mean",geom = "crossbar", width = 0.5, colour = "red")



hist(data$prop_diff) #after we run linear model test residuals
shapiro.test(data$prop_diff) #after we run linear model test residuals
ggpubr::ggqqplot(data$prop_diff)

####FINAL ANALYSIS####

#add ID in mixed model?
lm1 <- lm(prop_diff~comparision+origin+animation.x+trial_num+trial_day+start_time+trial_lane.x, data=data)
lm2 <- lm(prop_diff~comparision+origin+animation.x+trial_num+start_time, data=data)
lm3 <- lm(prop_diff~comparision*origin, data=data)

ggplot(data=data, aes(x=origin,y=prop_diff)) + 
  geom_violin(aes(col=comparision)) + 
  geom_point(aes(col=comparision), position = position_jitterdodge()) + 
  stat_summary(fun = "mean",geom = "crossbar", width = 0.5, colour = "red")



summary(lm1)
summary(lm2)
summary(lm3)

data2 <- subset(data,origin=="WC" & trial_lane...4 !="12_B")

t.test(prop_diff~comparision, data=data2, paired=TRUE)



####
FGSvsWT_mean$comparision <- "all males"
#data4molly <- rbind(FGSvsWT_mean,FGSvsWT_mean2)
data4molly$comparision <- gsub("JvJ","non-showy males",data4molly$comparision)
data4molly$comparision <- gsub("BvB","showy males",data4molly$comparision)

#data4molly_renamed <- data4molly %>% 
#  dplyr::rename(
#    individual = trial_lane...4,
#    prop_time_FGS = FGS_zone_,
#    prop_time_WT = WT_zone)

#data4molly_renamed$difference <- data4molly_renamed$prop_time_WT -data4molly_renamed$prop_time_FGS
  
write_csv(data, file.path("~/Desktop/Stanford/falsegravid/behavior/", "female_preference_proportions_power_curve.csv"))

wilcox.test(subset(data4molly, comparision=="all males")$FGS_zone_prop, subset(data4molly, comparision=="all males")$WT_zone_prop, paired = TRUE)
#t.test(subset(data4molly, comparision=="all males")$FGS_zone_sec, subset(data4molly, comparision=="all males")$WT_zone_sec, paired = TRUE)


FGSvsWT_mean_long <- FGSvsWT_mean2 %>% gather(phenotype, zone_time, FGS_zone_prop:WT_zone_prop, factor_key=TRUE)

FGSvsWT_mean_long$phenotype <- recode_factor(FGSvsWT_mean_long$phenotype, "FGS_zone_prop" = "FGS animation", 
                                             "WT_zone_prop" = "non-FGS animation")


t.test(FGSvsWT_mean$FGS_zone, FGSvsWT_mean$WT_zone, paired = TRUE)
#t.test(subset(FGSvsWT_mean_long, phenotype == "ratio")$zone_time, mu = 1)

FGSvsWT_mean <- FGSvsWT %>% group_by(trial_lane...4,) %>%
  dplyr::summarise(across(diff, ~ mean(.x, na.rm = TRUE)))

sort(FGSvsWT_mean$diff)
t.test(FGSvsWT_mean$diff, mu = 0, alternative = "two.sided")

hist(FGSvsWT_mean$diff, breaks=10)

####power test####

p_values <- NULL
FGSvsWT_mean_diff <- na.omit(FGSvsWT_mean$diff)
for (i in 1:10000) {
  t_test_result <- t.test(
    rnorm(length(FGSvsWT_mean_diff),
          mean = mean(FGSvsWT_mean_diff),
          sd = sd(FGSvsWT_mean_diff)), 
    mu = 0, alternative = "two.sided")
  p_values[i] <- t_test_result$p.value
}

#print(p_values)
sum(p_values < 0.05)/length(p_values)



#install.packages("Rmisc")
library(Rmisc)
FGSvsWT_mean_long_summary <- summarySE(FGSvsWT_mean_long, measurevar="zone_time", groupvars=c("phenotype"))

ggplot() +
  geom_point(data=FGSvsWT_mean_long_summary, aes(x=phenotype, y=zone_time, col=phenotype), size=6) + 
  geom_errorbar(data=FGSvsWT_mean_long_summary, aes(x=phenotype, ymin=zone_time-2*se, ymax=zone_time+2*se, col=phenotype), width=.1) + 
  geom_jitter(data=subset(FGSvsWT_mean_long), aes(x=phenotype,y=zone_time, col=phenotype), width = 0.1, alpha=0.5) +
  ylab("prop trial spent in zone") +
  xlab("zone") +
  scale_color_manual(values = c("#8538a8","darkgrey")) +
  theme_bw()

ggplot() +
  geom_point(data=FGSvsWT_mean_long_summary, aes(x=phenotype, y=zone_time, col=phenotype), size=3) + 
  geom_errorbar(data=FGSvsWT_mean_long_summary, aes(x=phenotype, ymin=zone_time-2*se, ymax=zone_time+2*se, col=phenotype), width=.1) + 
  geom_jitter(data=subset(FGSvsWT_mean_long), aes(x=phenotype,y=zone_time, col=phenotype), width = 0.1, alpha=0.5) +
  ylab("female preference (association time)") +
  xlab("") +
  scale_y_continuous(limits=c(0,0.6),breaks = c(0,0.2,0.4,0.6)) +
  theme(legend.position = "NA") +
  scale_color_manual(values = c("#8538a8","#ECB602")) +
  theme_bw() +
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_))

  
ggplot(data=subset(FGSvsWT_mean_long)) +
  geom_jitter(aes(x=phenotype,y=zone_time, col=phenotype), width = 0.1, alpha=0.5) +
  geom_boxplot() +
  ylab("prop trial spent in zone") +
  xlab("zone") +
  scale_color_manual(values = c("#8538a8","#ECB602")) +
  theme_bw()





FGSvsWT_mean_long2 <- FGSvsWT_mean2 %>% gather(phenotype, zone_time, FGS_zone:WT_zone, factor_key=TRUE)


ggplot(FGSvsWT_mean_long2, aes(col=phenotype, y=zone_time, x=comparision)) +
  geom_boxplot()+ geom_point(position = position_jitterdodge())

library(ggpubr)

ggboxplot(FGSvsWT_mean_long2, x = "comparision", y = "zone_time",
          color = "phenotype", palette = c("#8538a8","#ECB602"),
          add = "jitter") +
  stat_compare_means(aes(group = phenotype), paired= TRUE, label = "p.format")


model <- lm(zone_time ~ phenotype + comparision, data=FGSvsWT_mean_long2)
anova(model)


BvJ <- subset(trial_summary3, comparision == "BvJ")
BvJ$B_zone <- ifelse(BvJ$lifestage_zone_1=="bourgeois", BvJ$z1_trial, BvJ$z2_trial)
BvJ$J_zone <- ifelse(BvJ$lifestage_zone_1=="juvenile", BvJ$z1_trial, BvJ$z2_trial)

BvJ_long <- BvJ %>% dplyr::select(trial_num,B_zone:J_zone) %>% dplyr::gather(phenotype, zone_time, B_zone:J_zone, factor_key=TRUE)

ggplot(FGSvsWT_mean_long2, aes(col=phenotype, y=zone_time, x=comparision)) +
  geom_boxplot()+ geom_point(position = position_jitterdodge())

t.test(BvJ$B_zone, BvJ$J_zone, paired=TRUE)


ggboxplot(BvJ_long, x = "phenotype", y = "zone_time",
          color = "phenotype", palette = c("#8538a8","#ECB602"),
          add = "jitter") +
  stat_compare_means(aes(group = phenotype), paired=TRUE, label = "p.format")


# beginning to combine the two trials for each individual

# first need to remove individuals that were unresponsive in both trials
# and any individuals that are marked to "exclude" for other reasons (e.g., issues during trials)
# then test for side bias and remove any fish that spend >80% of their time on one side

# mark unresponsive individuals
# total trial time adds up to 0 when NAs are removed (unresponsive fish have NAs, so total==0)
filter_df <- trial_summary %>%
  group_by(fish_ID) %>%
  summarize(responsiveness = ifelse(sum(trial_time, na.rm=TRUE)==0, "unresponsive", "responsive"))
trial_summary <- left_join(trial_summary, filter_df, by="fish_ID")

# remove individuals that are unresponsive or marked to exclude
trial_summary_keep <- trial_summary %>%
  filter(responsiveness=="responsive") %>%
  filter(exclude=="N")
  
# mark individuals that exhibit side bias and exclude
filter_df2 <- trial_summary_keep %>%
  group_by(fish_ID) %>%
  summarize(z1_bias = sum(z1_total)/sum(total_time),
            z2_bias = sum(z2_total/sum(total_time)),
            side_biased = ifelse(max(z1_bias, z2_bias)>0.80, "biased", "unbiased"))
trial_summary_keep <- left_join(trial_summary_keep, filter_df2, by="fish_ID")
trial_summary_keep <- trial_summary_keep %>%
  filter(side_biased=="unbiased")




# combining the different trials to calculate association times

# combine the preferences from the two trials for A trial lanes
#A trial lanes for trial #1
combined_trials_A1 <- trial_summary_keep %>%
  group_by(fish_ID) %>%
  filter(trial_lane=="A") %>%
  filter(trial_num==1) %>%
  summarize(cortezi1=z2_trial, birchmanni1=z1_trial, trial_time1=trial_time, 
            trial_lane=trial_lane, species=species, date=date, time=time,
            distance1=distance, speed1=speed) %>%
  replace(is.na(.),0)
#A trial lanes for trial #2
combined_trials_A2 <- trial_summary_keep %>%
  group_by(fish_ID) %>%
  filter(trial_lane=="A") %>%
  filter(trial_num==2) %>%
  summarize(cortezi2=z1_trial, birchmanni2=z2_trial, trial_time2=trial_time,
            distance2=distance, speed2=speed) %>%
  replace(is.na(.),0)
#combine the two dfs
combined_trials_A <- left_join(combined_trials_A1, combined_trials_A2, by="fish_ID") %>%
  mutate(time_cortezi=cortezi1+cortezi2,
         time_birchmanni=birchmanni1+birchmanni2,
         total_time=trial_time1+trial_time2,
         prop_cortezi=time_cortezi/total_time,
         prop_birchmanni=time_birchmanni/total_time,
         prop_either=(time_cortezi+time_birchmanni)/total_time,
         SOP=(time_cortezi-time_birchmanni)/(time_cortezi+time_birchmanni),
         avg_distance=ifelse(distance1>0 & distance2>0, (distance1+distance2)/2, distance1+distance2), #addition works because one will be 0
         avg_speed=ifelse(speed1>0 & speed2>0, (speed1+speed2)/2, speed1+speed2)) %>% 
  select(date, time, fish_ID, trial_lane, species, prop_cortezi, prop_birchmanni, 
         prop_either, total_time, SOP, avg_distance, avg_speed)
      
# combine the preferences from the two trials for B trial lanes
#B trial lanes for trial #1
combined_trials_B1 <- trial_summary_keep %>%
  group_by(fish_ID) %>%
  filter(trial_lane=="B") %>%
  filter(trial_num==1) %>%
  summarize(cortezi1=z1_trial, birchmanni1=z2_trial, trial_time1=trial_time,
            trial_lane=trial_lane, species=species, date=date, time=time, 
            distance1=distance, speed1=speed) %>%
  replace(is.na(.),0)
#B trial lanes for trial #2
combined_trials_B2 <- trial_summary_keep %>%
  group_by(fish_ID) %>%
  filter(trial_lane=="B") %>%
  filter(trial_num==2) %>%
  summarize(cortezi2=z2_trial, birchmanni2=z1_trial, trial_time2=trial_time,
            distance2=distance, speed2=speed) %>%
  replace(is.na(.),0)
#combine the two dfs
combined_trials_B <- left_join(combined_trials_B1, combined_trials_B2, by="fish_ID") %>%
  mutate(time_cortezi=cortezi1+cortezi2,
         time_birchmanni=birchmanni1+birchmanni2,
         total_time=trial_time1+trial_time2,
         prop_cortezi=time_cortezi/total_time,
         prop_birchmanni=time_birchmanni/total_time,
         prop_either=(time_cortezi+time_birchmanni)/total_time,
         SOP=(time_cortezi-time_birchmanni)/(time_cortezi+time_birchmanni),
         avg_distance=ifelse(distance1>0 & distance2>0, (distance1+distance2)/2, distance1+distance2), #addition works because one will be 0
         avg_speed=ifelse(speed1>0 & speed2>0, (speed1+speed2)/2, speed1+speed2)) %>% 
  select(date, time, fish_ID, trial_lane, species, prop_cortezi, prop_birchmanni, 
         prop_either, total_time, SOP, avg_distance, avg_speed)

# combine the results from the A and B trial lanes into a single df
all_combined_trials <- rbind(combined_trials_A, combined_trials_B)







### output relevant dfs for subsequent use
#full dataset with latency, total recorded time in zones, total preference time, etc.
#each individual is separated across 2 rows
#no individuals are excluded
write.table(trial_summary, "./processed-ethovision-data/pheromone-overall-trial-summary.txt", quote=FALSE, sep="\t", row.names=FALSE)

#filtered dataset with latency, total recorded time in zones, total preference time, etc.
#each individual is separated across 2 rows
write.table(trial_summary_keep, "./processed-ethovision-data/pheromone-filtered-trial-summary.txt", quote=FALSE, sep="\t", row.names=FALSE)

#filtered dataset that includes overall cue "preference" for each individual
#excludes individuals that don't pass conditionals, etc.
write.table(all_combined_trials, "./processed-ethovision-data/pheromone-trial-preference.txt", quote=FALSE, sep="\t", row.names=FALSE)





##### MALE AGGRESSION BEHAVIOR ANALYSIS####
library(tidyverse)
library(Rmisc)
male_agg_means <- read.csv("/Users/trisdodge/Downloads/malemale_ag_means.csv")
male_agg_means <- male_agg_means %>% mutate(target = case_when(target == "fgs" ~ "FGS", target == "wt" ~ "non-FGS", TRUE ~ as.character(target)))

male_agg_means_mean <- male_agg_means %>% group_by(target) %>%
  summarise(~ mean(.x, na.rm = TRUE))

male_agg_means_mean <- summarySE(male_agg_means, measurevar="chase_tot_d", groupvars=c("target"))

male_agg_means_diff <- subset(male_agg_means, target != "female") %>% group_by(fm_id) %>%
  dplyr::summarise(across(hide_d:ag_c, ~ mean(.x, na.rm = TRUE)))

hist(subset(male_agg_means, target != "female")$chase_tot_d, breaks=10)
ggpubr::ggqqplot(subset(male_agg_means, target != "female")$chase_tot_d)
hist(male_agg_means_diff$chase_tot_d)
ggpubr::ggqqplot(male_agg_means_diff$chase_tot_d)




male_gonc_means_mean <- summarySE(male_agg_means, measurevar="app_c", groupvars=c("target"))

ggplot() +
  geom_point(data=subset(male_gonc_means_mean, target != "female"), aes(x=target, y=app_c, col=target), size=3) + 
  geom_errorbar(data=subset(male_gonc_means_mean, target != "female"), aes(x=target, ymin=app_c-2*se, ymax=app_c+2*se, col=target), width=.1) + 
  geom_jitter(data=subset(male_agg_means, target != "female"), aes(x=target,y=app_c, col=target), width = 0.1, alpha=0.5) +

  #ylab(expression(paste(italic("X. birchmanni"), " male \n aggression (chase duration)",sep=""))) +
  
  ylab("male approaches") +
  xlab("") +
  theme(legend.position = "NA") +
  scale_color_manual(values = c("#8538a8","#ECB602")) +
  theme_bw() +
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_)) +
  stat_compare_means(data=subset(male_agg_means, target != "female"), aes(x = target,y=app_c), method="wilcox.test", paired=TRUE, label = "p.format")


ggplot() +
  geom_point(data=subset(male_agg_means_mean, target != "female"), aes(x=target, y=chase_tot_d, col=target), size=3) + 
  geom_errorbar(data=subset(male_agg_means_mean, target != "female"), aes(x=target, ymin=chase_tot_d-2*se, ymax=chase_tot_d+2*se, col=target), width=.1) + 
  geom_jitter(data=subset(male_agg_means, target != "female"), aes(x=target,y=chase_tot_d, col=target), width = 0.1, alpha=0.5) +
  ylab(expression("italic('X. birchmanni') male \n aggression (chase duration)")) +
  
  #ylab(expression(paste(italic("X. birchmanni"), " male \n aggression (chase duration)",sep=""))) +
  
  #ylab("male aggression (chase duration)") +
  xlab("") +
  theme(legend.position = "NA") +
  scale_color_manual(values = c("#8538a8","#ECB602")) +
  theme_bw() +
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_))

male_agg_means_wide <- male_agg_means %>% 
  dplyr::select(c("fm_id","target","chase_tot_d"))  %>% 
  spread(target, chase_tot_d)

wilcox.test(male_agg_means_wide$FGS, male_agg_means_wide$`non-FGS`, paired = TRUE)

mean(male_agg_means_wide$fgs/male_agg_means_wide$wt)





male_court_means_mean <- summarySE(male_agg_means, measurevar="court_d", groupvars=c("target"))


ggplot() +
  geom_point(data=subset(male_court_means_mean, target != "female"), aes(x=target, y=court_d, col=target), size=3) + 
  geom_errorbar(data=subset(male_court_means_mean, target != "female"), aes(x=target, ymin=court_d-2*se, ymax=court_d+2*se, col=target), width=.1) + 
  geom_jitter(data=subset(male_agg_means, target != "female"), aes(x=target,y=court_d, col=target), width = 0.1, alpha=0.5) +
  ylab(expression("italic('X. birchmanni') male \n attention (courtship duration)")) +
  
  #ylab(expression(paste(italic("X. birchmanni"), " male \n aggression (chase duration)",sep=""))) +
  
  #ylab("male aggression (chase duration)") +
  xlab("") +
  theme(legend.position = "NA") +
  scale_color_manual(values = c("#8538a8","#ECB602")) +
  theme_bw() +
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_))


male_court_means_wide <- male_agg_means %>% 
  select(c("fm_id","target","court_d"))  %>% 
  spread(target, court_d)

wilcox.test(male_court_means_wide$FGS, male_court_means_wide$`non-FGS`, paired = TRUE)

male_court_means_wide$fgs - male_court_means_wide$wt
