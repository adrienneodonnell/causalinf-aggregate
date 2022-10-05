# Load packages
library(data.table)
library(plyr)
library(ggplot2)
library(dplyr)
library(stringr)

########################################################
# 1. Build function to generate simulated daily EHR data
########################################################
datagenerate <- function(i){
  # Number of timepoints
  K <- 60
  # Patient ID number
  id <- as.numeric(i) 
  # Baseline time
  t0 <- 0
  # Baseline data
    # Visit process for confounder L
    V <- 1
    # Continuous confounder L
    L <- rnorm(1, 0, sd=1)
    Llag <- L
    # Measured/observed continuous covariate LS: When V=1, LS=L
    LS <- L
    LSlag <- Llag
    # Binary treatment exposure indicator A
    A <- rbinom(1, 1, plogis(-5-3*LS))
    Alag <- A
    # Binary outcome indicator Y(t+1)
    Y <- rbinom(1, 1, plogis(-10-4*LS+1*A))
  # Generate vectors to build data frame
  id_     <- c(id)
  t0_     <- c(t0)
  V_      <- c(V)
  L_     <- c(L)
  LS_    <- c(LS)
  A_      <- c(A)
  Y_      <- c(Y)
  # Post-baseline data
  if ((K > 1) && (Y==0)){
    for (j in 2:K){
      # Time interval
      t0 <- j-1
      # Visit process for confounder L
      Vstar   <- rbinom(1, 1, plogis(-2))
      temp_V <- c(V_, Vstar)
      # Confounder L
      Lstar  <- rnorm(1, mean= -0.5*t0, sd=1)
      temp_L <- c(L_, Lstar)
      # Measured/observed confounder LS with LOCF when not measured
      if (Vstar == 1){
        LSstar <- Lstar
      }
      else{
        LSstar <- LS_[j-1]
      }
      temp_LS <- c(LS_, LSstar)
      # Treatment A       
      if (A_[j-1] == 1){
        Astar <- 1
      }
      else{
        Astar <- rbinom(1, 1, plogis(-5-3*LSstar))
      }
      temp_A <- c(A_, Astar)
      # TTE Outcome Y
      Ystar <- rbinom(1, 1, plogis(-10-4*LSstar+1*A))
      # Finalize new data to add to longitudinal data frame
      id_[j]        <- id
      t0_[j]        <- t0
      V_[j]         <- Vstar
      L_[j]        <- Lstar
      LS_[j]       <- LSstar
      A_[j]         <- Astar
      Y_[j]         <- Ystar
      
      # If outcome event realized, discontinue generating data for indiv.
      if (Ystar==1){
        break
      }
    }
  }  
  # Consolidate data in a single data frame
  temp_data <- data.frame(id = id_,
                          t0 = t0_,
                          V = V_,
                          L = L_,
                          LS = LS_,
                          A = A_,
                          Y = Y_)
  return(temp_data)
}


########################################################
# 2. Run function to generate simulated daily EHR data
########################################################

set.seed(50)
# Set sample size
n<-2000
# Call function
dat_d_l <- lapply(as.list(1:n), FUN=function(ind){
  datagenerate(ind)
})
# Save variables as a dataset
dat_d <- rbindlist(dat_d_l)

########################################################
# 3. Aggregate daily data into weekly data
########################################################

# Make a copy of the daily data
dat_pre1 <- dat_d

# For continuous covariate LS, create a copy without LOCF (we only want to aggregate observed LS)
# Call this LM
dat_pre2 <- dat_pre1 %>%
  mutate(LM = ifelse(V==1, LS, NA))%>%
  dplyr::select(id,t0,V,LS,LM,A,Y)

# Assign daily timepoints (t0) to intervals (t0_2) and identify the position within an interval of V and A
dat_pre3 <- dat_pre2 %>%
  mutate(Alag=lag(A),
         t0_2 = (t0 %/% 7),
         position=t0-t0_2*7,
         pos_V=ifelse(V==1,position,NA),
         pos_A_pre=ifelse(A==1 & (Alag==0 | is.na(Alag)==T),position,NA))

# Need the position of A for a whole interval, to compare with each V
dat_pre4 <- dat_pre3 %>%
  group_by(id,t0_2)%>%
  dplyr::summarize(pos_A=ifelse(max(pos_A_pre,na.rm=T)!=-Inf,max(pos_A_pre,na.rm=T),NA))

# Merge with daily data in order to compare position of V with A for whole interval
dat_pre5 <-merge(dat_pre3,
                 dat_pre4,
                 by=intersect(names(dat_pre3), names(dat_pre4)))

# Drop unnecessary variables
dat_pre6 <- select(dat_pre5,-c(Alag,pos_A_pre,position))

# If treatment initiated in interval s and there is a post-treatment covariate measurement
# in interval s, push that covariate measurement to interval s+1
dat_pre7 <- dat_pre6 %>%
  mutate(t0_2L=ifelse(V==1 & pos_V>pos_A & is.na(pos_A)==F,t0_2+1,t0_2))

# Take the mean of LM and the max of A, V, and Y within each interval
empty_frame2 <- dat_pre7 %>%
  group_by(id, t0_2L) %>%
  dplyr::summarize(ave_LM = ifelse(any(!is.na(LM)), mean(LM, na.rm=T), NA),
                   max_V = max(V), 
                   max_A = max(A),
                   max_Y = max(Y),
                   .groups="drop") %>%
  tidyr::fill(ave_LM) #Carry forward continuous covariate information

# Rename variables
frame3<-dplyr::rename(empty_frame2,t0=t0_2L,L=ave_LM,V=max_V,A=max_A,Yn=max_Y)

# When we aggregated  (using dplyr::summarize), if Y(s)=1 we need to move Y to the (s-1)th row, since 
# the column is Y(t+1), not Y(t) 
# Save the ID and time which has Y(t+1)=1
tofY2 <- frame3[which(frame3$Yn==1),] #select the rows with Y=1
tofY2$tY <- tofY2$t0-1 #create variable tY that shows the time Y=1 occurred
tofY3 <- tofY2[c("id","tY","Yn")] #keep only the necesary variables
tofY4<- tofY3 %>% mutate(tY = ifelse(tY==-1,0,tY)) #if Y occurred at (t+1)=0, force it to (t+1)=1
tofY5<-dplyr::rename(tofY4,t0=tY,Y=Yn) #rename variables to be compatible with merging back

# merge the dataset with id and Y(t+1) with the complete dataset
frame4<- full_join(frame3,tofY5,by=c("id","t0"))

# Set Y=NA to Y=0
frame5 <- frame4 %>% mutate(Y = ifelse(is.na(Y), 0, Y))

# Save the location that Y=0
loc<-dplyr::rename(tofY5,Ypos1=t0)
loc2<-loc %>% mutate(Ypos=ifelse(Ypos1>16,16,Ypos1))

# merge location of Y with complete dataset and delete time periods after Y=1
frame5a<- full_join(frame5,loc2[, c("id","Ypos")],by=c("id"))
frame6 <- frame5a %>% filter(Ypos>=t0 | (is.na(Ypos)==T & t0<=floor(119/7)-1)) #keep if any of these

# Order the columns
col_order <- c("id","t0","V","L","A","Y")
dat_w <- frame6[, col_order]
