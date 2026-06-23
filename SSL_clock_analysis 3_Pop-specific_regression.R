library(ggplot2)
library(tidyverse)
#library(ggpmisc)

#Read data 
df = read.csv('Data_3_CpG_scores_dr.csv')
#kh adds: only 10 sites
df <- df[, c(1:3,11,12,17,20,21,22,26,28,29,31)]
#Remove rows that had any Nan
df <- na.omit(df)
#Remove the outlier - necropsied male
df <- df[-66,]

#Data in vectors
Age=df$Age
Pop=df$Location1
Sex=df$Sex
C1=df$MDGA2_58_ratio
C2=df$HSPA2_41_ratio
C3=df$MDGA2_73_ratio
C4=df$SLC12A5_61_ratio
C5=df$SLC12A5_64_ratio
C6=df$OTP_62_ratio
C7=df$SLC12A5_54_ratio
C8=df$HSPA2_60_ratio
C9=df$MDGA2_40_ratio
#C10=df$HSPA2_34_ratio  #Removed from clock
#C11=df$HCN1_24_ratio.  #Removed from clock
#C12=df$SLC12A5_59_ratio #Removed from clock
C13=df$SLC12A5_93_ratio

#Multiple regression of all points
modelM=lm(formula=Age~C1+C2+C3+C4+C5+C6+C7+C8+C9+C13)
predictions <- predict(modelM)
modelp=lm(formula=predictions~Age)

modelp2=lm(formula=predictions~Age*Pop)
modelp3=lm(formula=predictions~Age+Pop)
AIC(modelp,modelp2, modelp3)


#Plotting for multiple linear regression
pdata<-data.frame(Age,predictions)
plot1<-ggplot(pdata,aes(y=predictions,x=Age, color=Pop))+geom_smooth(method="lm", fill = NA)+
  geom_point(shape=1, size=4)+xlab('Known age')+ylab('Predicted age')+
  #stat_poly_eq(formula = y ~ x, aes(label = after_stat(eq.label))) +
  xlim(4, 26)+ylim(4, 26) +
  scale_color_manual(values = c("Eastern DPS" = "darkred", "Western DPS" = "blue")) +
  theme(axis.title = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(legend.position = c(0.8, 0.15)) 

ggsave("pop_all.pdf", plot1)

##Subset the data to maximum age of 16, repeat analysis
df <- subset(df, Age < 16.5,)

#Data in vectors
Sex=df$Sex
Pop=df$Location1
Age=df$Age
Pop=df$Location1
Sex=df$Sex
C1=df$MDGA2_58_ratio
C2=df$HSPA2_41_ratio
C3=df$MDGA2_73_ratio
C4=df$SLC12A5_61_ratio
C5=df$SLC12A5_64_ratio
C6=df$OTP_62_ratio
C7=df$SLC12A5_54_ratio
C8=df$HSPA2_60_ratio
C9=df$MDGA2_40_ratio
#C10=df$HSPA2_34_ratio  #Removed from clock
#C11=df$HCN1_24_ratio.  #Removed from clock
#C12=df$SLC12A5_59_ratio #Removed from clock
C13=df$SLC12A5_93_ratio

#Multiple regression of all points
modelM=lm(formula=Age~C1+C2+C3+C4+C5+C6+C7+C8+C9+C13)
predictions <- predict(modelM)
modelp=lm(formula=predictions~Age)

modelp2=lm(formula=predictions~Age*Pop)
modelp3=lm(formula=predictions~Age+Pop)
AIC(modelp,modelp2, modelp3)

#Plotting for multiple linear regression
pdata<-data.frame(Age,predictions)
plot2<-ggplot(pdata,aes(y=predictions,x=Age, color=Pop))+geom_smooth(method="lm", fill = NA)+
  geom_point(shape=1, size=4)+xlab('Known age')+ylab('Predicted age')+
  #stat_poly_eq(formula = y ~ x, aes(label = after_stat(eq.label))) +
  #xlim(7.5, 16.5)+ylim(7.5, 16.5) +
  scale_color_manual(values = c("Eastern DPS" = "darkred", "Western DPS" = "blue")) +
  theme(axis.title = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10)) +
  theme(legend.position = "none") 

ggsave("pop_8_16.pdf", plot2)