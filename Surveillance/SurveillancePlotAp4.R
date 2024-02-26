############################################################
### ### ### Surveillance analysis 
### ### ### 10/13/23
### ### ### Annabelle Stanley
############################################################

### Load in r packages
library(tidyverse)
library(readxl)
library("xlsx")

### Load in data

#Appendix4Data <- read_excel("Appendix4Data.xlsx", head = F)
data <- read.xlsx(file = "Appendix4Data.xlsx", sheetIndex = "Sheet1", header=F)

colnames(data) <- paste0(data[1,], sep = "_", data[2,])
data <- data[-c(1,2),]
datal <- data %>% pivot_longer(cols=c("90_10",     "90_5" ,     "90_2" ,     "90_1"   ,   "90_0.1",    "95_10",    
                                       "95_5" ,     "95_2" ,     "95_1"  ,    "95_0.1",    "99_10"  ,   "99_5"  ,    "99_2",     
                                       "99_1"   ,   "99_0.1"),
                               names_to='Prob2',
                               values_to='samp')



temp <- datal %>% 
  separate(col = Prob2, into = c("Prob1pos", "prevalence"), sep = "_")

colnames(temp) <- c("PopnSize",  "Prob1pos",   "prevalence", "samp"  )

temp$prevalence <- as.numeric(temp$prevalence)

## Population is going to be greater than 50, but less than 10,000
temp$PopnSize <- sapply(temp$PopnSize, as.numeric)
temp <- temp %>% filter(PopnSize > 50) %>% 
  filter(PopnSize < 10000) 

ggplot(temp) +
  geom_point(aes(x = prevalence, y = samp, color =Prob1pos ))



### ### ###

## try simulating data...

namemelater <- function(prevalence = .01, confidence = .99){
  log(-(confidence - 1))/(log(1 - prevalence))
}


range_prev <- seq(0, .05, by = .0005) #0-20

startgraph <- as.data.frame(range_prev)

startgraph2 <- startgraph %>% 
  mutate(nosamp_99 = namemelater(prevalence = range_prev)) %>% ## 99%
  mutate(nosamp_95 = namemelater(prevalence = range_prev, confidence = .95)) %>% ## 95%
  mutate(nosamp_90 = namemelater(prevalence = range_prev, confidence = .9)) %>%   ## 90%
  mutate(range_prev = range_prev*100) # make into percent

colors <- c("99%" = "blue", "95%" = "red", "90%" = "orange")

ggplot(startgraph2) +
  geom_point(aes(x = range_prev, y = nosamp_99, color = "99%")) + 
  geom_point(aes(x = range_prev, y = nosamp_95, color = "95%")) + 
  geom_point(aes(x = range_prev, y = nosamp_90, color = "90%"))+
  labs(x = "Prevalence (%)", y = "No. Samples", color = "Legend")+
  scale_color_manual(values = colors)
