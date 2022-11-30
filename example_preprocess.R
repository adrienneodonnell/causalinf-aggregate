########################################################
library(dplyr)

#Step 0. Read in data. For continuous covariate L, create a copy called LS without LOCF (we only want to aggregate observed L).
dat_pre1 <- read.csv(file = 'example_dat.csv')
dat_pre2 <- dat_pre1 %>% mutate(LS = ifelse(V==1, L, NA))
head(dat_pre2)

#Step 1. Assign granular timepoints (days, t) to intervals (7 days/1 week, w).
#Identify the position within an interval of V and initiation of A.
dat_pre3 <- dat_pre2 %>% 
         mutate(Alag=lag(A),
         w = (t %/% 7),
         position=t-w*7,
         pos_V=ifelse(V==1,position,NA),
         pos_A_pre=ifelse(A==1 & (Alag==0 | is.na(Alag)==T),position,NA))
head(dat_pre3)

#Step 2. Expand the position of A to the entire interval in order to compare with each V.
dat_pre4 <- dat_pre3 %>% group_by(id,w)%>%
  dplyr::summarize(pos_A=ifelse(max(pos_A_pre,na.rm=T)!=-Inf,max(pos_A_pre,na.rm=T),NA))
head(dat_pre4)

#Step 3. Merge with complete data. Drop unnecessary variables.
dat_pre5 <- merge(dat_pre3,
                 dat_pre4,
                 by=intersect(names(dat_pre3), names(dat_pre4)))
dat_pre6 <- select(dat_pre5,-c(Alag,pos_A_pre,position))
head(dat_pre5)

#Step 4. Compare position of measured confounder with treatment within an interval.
#If treatment is initiated in interval s and there is a post-treatment covariate measurement in interval s, push that covariate measurement to interval s+1.
dat_pre7 <- dat_pre6 %>%
  mutate(w_new=ifelse(V==1 & pos_V>pos_A & is.na(pos_A)==F,w+1,w))
head(dat_pre7)

#Step 5. Aggregate information across the new interval. Take the mean of LS and the max of A, V, and Y. Take the sum of the number of visits in an interval for PNA. Rename variables.
dat_pre8 <- dat_pre7 %>%
  group_by(id, w_new) %>%
  dplyr::summarize(LS_agg=mean(LS,na.rm=T),
                   A_agg=max(A),
                   V_agg=max(V),
                   Y_agg=max(Y),
                   V_sum=sum(V))
dat_pre9<-dplyr::rename(dat_pre8,t=w_new,LS=LS_agg,V=V_agg,A=A_agg,Yn=Y_agg)
head(dat_pre9)

#Step 6. Put Y at time s+1 in the row for time s. 
  #6A. Select the rows with Y(s)=1
  events <- dat_pre9[which(dat_pre9$Yn==1),]
  #6B. Create variable t-1, called t1. If Y(s+1)=1, then t1=s.
  events$t1 <- events$t-1
  #6C. Keep on necessary variables.
  events <- events[c("id","t1","Yn")]
  #6D. If Y occurred at baseline, move it to the first interval.
  events <- events %>% mutate(t1 = ifelse(t1==-1,0,t1))
  #6E. Rename variables to be compatible with merging back.
  events <- events %>% rename(t=t1,Y=Yn)
  #6F. Merge with complete data.
  dat_pre10 <- full_join(dat_pre9,events,by=c("id","t"))
  #6G. Set Y=NA to Y=0.
  dat_pre11 <- dat_pre10 %>% mutate(Y = ifelse(is.na(Y), 0, Y))

#Step 7. Drop unecessary variables and finalize aggregated dataset.
dat_agg <- select(dat_pre11,-c(Yn))



#To compute PNA:
#Numerator counts number of person-times with 0 or 1 visit
num <- sum(table(dat_agg$V_sum)[1]+table(dat_agg$V_sum)[2])
#Denominator counts number of person-times
denom <- sum(table(dat_agg$V_sum))
PNA <- num/denom 
print(PNA)
ifelse(PNA<0.75,print("Try a narrower interval size."),print("Interval size sufficient (PNA>75%)."))
