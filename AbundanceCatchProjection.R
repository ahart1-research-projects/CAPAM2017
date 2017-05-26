AgeBasedMSY <- function(data, ZData=NULL, NAAData=NULL, NAgeClass=NULL, WeightAtAge=NULL, MaturityAtAge=NULL, MortalityAtAge=NULL, CAAData=NULL, Nyr=NULL, RecruitOption="Option1", StdDev=NULL, DegFreedom=NULL){
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
  names(NAAProjection) <- colnames(NAA)
  SSBProjection <- rep(NA, Nyr)
  Recruitment <- rep(NA, Nyr)
  CAAProjection <- rep(NA, NAgeClass) 
  names(CAAProjection) <- colnames(NAA)
  
  ########## Forward Projection ##########
  for(t in 1:Nyr){
    if(RecruitOption=="Option1"){
      Recruitment[t] <- rnorm(1, mean = NAA[nrow(NAA),1], sd=StdDev) # Random walk
    } else if(RecruitOption=="Option2"){
      Recruitment[t] <- rnorm(1, mean = mean(NAA[ ,1]), sd=StdDev) # Random about the mean
    } else if(RecruitOption=="Option3"){
      Recruitment[t] <- rt(1, df=DegFreedom) # Student-T
    }
    
    #Set age 1 abundance equal to recruitment based on RecruitmentOption chosen
    NAAProjection[1] <- Recruitment[t]
    
    # NEEd to update ZMortality also!!!???????
    UpdateZ <- rep(0.2, NAgeClass)
    #UpdateZ <- # Calculation here ???????
    ZMortality <- rbind(ZMortality, UpdateZ)
    
    #Calculate abundance-at-age for year 2-Nyr
    for(a in 2:NAgeClass){
      NAAProjection[a]<-NAA[nrow(NAA), a-1]*exp(-ZMortality[nrow(ZMortality), a-1])
    }
    NAA <- rbind(NAA, NAAProjection)
    
    # Catch-at-age projection
    for(a in 1:NAgeClass){
      CAAProjection[a] <- NAAProjection[a]*MortalityAtAge[t, a]*(1-exp(-ZMortality[t, a]))/ZMortality[t,a]
    }
    CAA <- rbind(CAA, CAAProjection)
  }
  
  return(list(NAAWithProjection=NAA, Recruitment=Recruitment, CAAWithProjection=CAA))
}




################### Run this function ##########
ZDataFile <- read.csv("Z.csv")
ZDataFile <- ZDataFile[,-1]
NAADataFile <- read.csv("NAA_mat.csv")
NAADataFile <- NAADataFile[,-1]
MaturityFile <- read.csv("maturity.csv")
MaturityFile <- MaturityFile[,-1]
MortalityFile <- read.csv("mort.csv")
MortalityFile <- MortalityFile[,-1]
WeightFile <- read.csv("waa.csv")
WeightFile <- WeightFile[,-1]
CAADataFile <- read.csv("CAA_mat.csv")
CAADataFile <- CAADataFile[,-1]

TestOption2 <- AgeBasedMSY(ZData = ZDataFile, NAAData = NAADataFile, NAgeClass = 9, WeightAtAge = WeightFile, MaturityAtAge = MaturityFile, CAAData = CAADataFile, Nyr = 5, RecruitOption = "Option2", StdDev=0.6160962, MortalityAtAge=MortalityFile)
