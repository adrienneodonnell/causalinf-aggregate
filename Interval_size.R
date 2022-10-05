
############################ Try first partition ############################
#Make a copy of the daily data
dat_pre1 <- dat_d

#For continuous covariate LS1, create a copy without LOCF (we want to aggregate observed L)
dat_pre2 <- dat_pre1 %>%
  mutate(LM1 = ifelse(V==1, LS1, NA))%>%
  dplyr::select(id,t0,V,LS1,LM1,A,Y,L1b)

#Assign daily timepoints (t0) to intervals (t0_2) and identify the position within an interval of V and A
dat_pre3 <- dat_pre2 %>%
  mutate(Alag=lag(A),
         t0_2 = (t0 %/% 8),
         position=t0-t0_2*8,
         pos_V=ifelse(V==1,position,NA),
         pos_A_pre=ifelse(A==1 & (Alag==0 | is.na(Alag)==T),position,NA))

#Need the position of A for a whole interval, to compare with each V
dat_pre4 <- dat_pre3 %>%
  group_by(id,t0_2)%>%
  dplyr::summarize(pos_A=ifelse(max(pos_A_pre,na.rm=T)!=-Inf,max(pos_A_pre,na.rm=T),NA))

#Merge with daily data in order to compare position of V with A for whole interval
dat_pre5 <-merge(dat_pre3,
                 dat_pre4,
                 by=intersect(names(dat_pre3), names(dat_pre4)))

#Drop unnecesary variables
dat_pre6 <- select(dat_pre5,-c(Alag,pos_A_pre,position))

#If treatment initiated in interval s and there is a post-treatment covariate measurement
#in interval s, push that covariate measurement to interval s+1
dat_pre7 <- dat_pre6 %>%
  mutate(t0_2L=ifelse(V==1 & pos_V>pos_A & is.na(pos_A)==F,t0_2+1,t0_2))

#Take the mean of LM1 and the max of A, V, and Y within each interval
sum2 <- dat_pre7 %>%
  group_by(id, t0_2L) %>%
  dplyr::summarize(sum_V = sum(V))

#hist(sum3$max_V)
#Proportion of patients with 1 visit per interval
prop1<-sum((count(sum2$sum_V)[1,2]+count(sum2$sum_V)[2,2])/sum(count(sum2$sum_V)[,2]))
print(prop1)
ifelse(prop1<0.75,print("make interval size narrower"),print("make interval size wider"))


hist(sum3$max_V)
