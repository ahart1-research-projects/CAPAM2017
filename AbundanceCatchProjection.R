AgeBasedMSY <- function(data, ZData=NULL, NAAData=NULL, NAgeClass=NULL, WeightAtAge=NULL, MaturityAtAge=NULL, CAAData=NULL, Nyr=NULL, RecruitOption="Option1", StdDev=NULL, DegFreedom=NULL){
  # This function projects Abundance-at-age and Catch-at-age forward in time
  # Args:
       # ZData: Matrix of mortality-at-age
       # NAAData: Matrix of abundance-at-age
       # NAgeClass: Number of age classes
       # WeightAtAge: Vector of weight-at-age 
       # MaturityAtAge: Vector of maturity-at-age
       # CAAData: Matrix of catch-at-age data
       # Nyr: Number of years to project forward
       # RecruitOption: Chooses method used for calculating recruitment
          # RecruitOption=="Option1": Random walk
          # RecruitOption=="Option2": Random about the mean
          # RecruitOption=="Option3": Student-T 
       # StdDev: Standard deviation required when RecruitOption 1 and 2 chosen
       # DegFreedom: Degrees of freedom required for when RecruitOption 3 chosen
  # Return:
       # A list containing a matrix of NAAData with appended projection Nyr long, a matrix of CAAData with appended projection Nyr long, and
       # a vector of projected recruitment
 
  NAA <- NAAData
  CAA <- CAAData
  ZMortality <- ZData
  # Set up storage 
  NAAProjection <- rep(NA, NAgeClass)
  SSBProjection <- rep(NA, Nyr)
  Recruitment <- rep(NA, Nyr)
  CAAProjection <- rep(NA, NAgeClass)
  
  ########## Forward Projection ##########
  for(t in 1:Nyr){
    #Define SSB for each year(t) and age(a)
  #   if(RecruitOption=="Option"){
  #     for(a in 1:Nyr){
  #       SSBProjection[t+1] <- NAAProjection[t,a]*MaturityAtAge[t, a]*WeightAtAge[t, a]
  #     }
  #   Recruitment[t+1] <- b1_par*SSB[t]/(SSB[t]+b2_par)
  # } else if(RecruitOption=="Option2"){} #??????????
  
    if(RecruitOption=="Option1"){
      Recruitment[t] <- rnorm(1, mean = NAA[nrow(NAA),1], sd=StdDev) # Random walk
    } else if(RecruitOption=="Option2"){
      Recruitment[t] <- rnorm(1, mean = mean(NAA[ ,1]), sd=StdDev) # Random about the mean
    } else if(RecruitOption=="Option3"){
      Recruitment[t] <- rt(1, df=DegFreedom) # Student-T
    }
    
    #Set age 1 abundance equal to recruitment based on RecruitmentOption chosen
    NAAProjection[t, 1] <- Recruitment[t]
    
    # NEEd to update ZMortality also!!!???????
    UpdateZ <- # Calculation here ???????
    ZMortality <- rbind(ZMortality, UpdateZ)
    
    #Calculate abundance-at-age for year 2-Nyr
    for(a in 2:NAgeClass){
      NAAProjection[t, a]<-NAA[nrow(NAA), a-1]*exp(-ZMortality[nrow(ZMortality), a-1])
    }
    NAA <- rbind(NAA, NAAProjection)
    
    # Catch-at-age projection
    for(a in 1:NAgeClass){
      CAAProjection[t,a] <- NAAProjection[t, a]*F_par[t, a]*(1-exp(-ZMortality[t, a]))/ZMortality[t,a]
    }
    CAA <- rbind(CAA, CAAProjection)
  }
  
  return(list(NAAWithProjection=NAA, Recruitment=Recruitment, CAAWithProjection=CAA))
}




################### Run this function ##########
ZDataFile <- read.csv("Z.csv")
NAADataFile <- read.csv("NAA_mat.csv")
