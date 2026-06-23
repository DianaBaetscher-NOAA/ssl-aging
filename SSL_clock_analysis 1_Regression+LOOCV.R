library(ggplot2)
library(tidyverse)

#Read data 
df = read.csv('data/Data_3_CpG_scores_dr.csv')
#kh adds: only 10 sites
df <- df[, c(1:3,11,12,17,20,21,22,26,28,29,31)]
#Remove rows that had any Nan
df <- na.omit(df)
#Remove column 'Sex' - number 2
df <- subset(df, select = -c(2:3))
#Remove the outlier - necropsied male
df <- df[-66,]


#Data in vectors
Age=df$Age
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
summary(modelp)
coeff=summary(modelp)$coefficients
yInt=coeff[1,1]
yIntP=coeff[1,4]
R2=summary(modelp)$adj.r.squared

#Plotting for multiple linear regression
pdata<-data.frame(Age,predictions)

plot1 <- ggplot(pdata, aes(y = predictions, x = Age)) + 
  geom_smooth(method="lm")+
  geom_point() +
  theme_bw() +
  xlab('Known age') + 
  ylab('Predicted age') +
  xlim(4, 26) +
  ylim(4, 26) +
  theme(
    axis.title = element_text(size = 14),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text = element_text(size = 12)
  )
 
ggsave("outputs/Calibration.png", plot1,  width = 5, height = 5)

#Perform LOOCV to get range of predictions around regression
removeOne <- function(dat,x) {
  list=seq(1,dat)
  x=x-1
  v1=c(list[1:x])
  x=x+2
  v2=c(list[x:length(list)])
  data=c(v1,v2)
  return (data)}

scale1=yInt/max(Age)
scaling=1+(1-R2)+scale1
scaling2=1+(1-R2)

nSamples=length(Age)
predictions2 <- c(1:nSamples)
for (z in 1:nSamples)
  {
  indices=removeOne(nSamples,z)
    modelM=lm(formula=Age[indices]~C1[indices]+C2[indices]+C3[indices]+C4[indices]+C5[indices]+C6[indices]+C7[indices]+C8[indices]+C9[indices]+C13[indices])
  newdata=data.frame(C1=C1[z],C2=C2[z],C3=C3[z],C4=C4[z],C5=C5[z],C6=C6[z],C7=C7[z],C8=C8[z],C9=C9[z],C13=C13[z])
  p=predict(modelM,newdata)
  predictions2[z]=p 
  }

differences <- c(1:nSamples)
for (z in 1:nSamples)
{
  diff=Age[z]-predictions2[z]
  if (diff<0){diff=diff*-1}
  differences[z]=diff
}
differences
summary(differences)
s = sd(differences)
t_star = qt(0.975, nSamples-1)
t.test(differences, conf.level=0.95)

x_bar = mean(differences)
mean=toString(signif(x_bar,4))
n = length(differences)
s = sd(differences)
ci_upper = x_bar + t_star*s/sqrt(n)
ci_lower = x_bar - t_star*s/sqrt(n)

for (z in 1:length(differences))
  {if (differences[z]<0){differences[z]=differences[z]*-1}
}
s = sd(differences)
t_star = qt(0.975, n-1)
LSD=signif(t_star*sqrt(2*s),4)
LSD=toString(LSD)

modelp2=lm(formula=predictions2~Age)
pred <- predict(modelp2) 

fit <- data.frame(x1=min(Age), y1=min(pred), x2=max(Age), y2=max(pred))
highL <- data.frame(x1=min(Age), y1=min(pred)+ci_lower, x2=max(Age), y2=max(pred)+ci_upper)
lowL <- data.frame(x1=min(Age), y1=min(pred)-ci_lower, x2=max(Age), y2=max(pred)-ci_upper)

#Plotting
pdata2<-data.frame(Age,predictions2)

plot2 <- ggplot(pdata,aes(y = predictions2, x = Age), show.legend = FALSE) +
  geom_point() +
  theme_bw() +
  xlab('Known age of n-1 samples') +
  ylab('Predicted age of each left out sample') +
  xlim(4, 26) + 
  ylim(4, 26) +
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2, colour = 'fit'), data=fit,show.legend = FALSE)+
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2, colour = 'bound', linetype = "bound"), data = lowL, show.legend = FALSE) +
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2, colour = 'bound', linetype = "bound"), data = highL, show.legend = FALSE) +
  scale_color_manual(values = c("fit" = "black", 'bound' = "blue")) + 
  scale_linetype_manual(values = c("fit" = "solid", "bound" = "dashed")) +
  theme(
    axis.title = element_text(size = 14),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text = element_text(size = 12)
    )


ggsave("outputs/loocv.png", plot2, width = 5, height = 5)

