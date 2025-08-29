## Estimate FE for Pop-Elches (2013)
# with matlab it is out of memory...

rm( list = ls())
library(tidyverse)
#install.packages("plm")
library(plm)
library(haven)

## Get the data
path <- "/Users/areguly6/Library/CloudStorage/OneDrive-CorvinusUniversityofBudapest/Research/RDD/develop/"
data_AER_1 <- read_dta(paste0(path,"data/Pop_Elches/data-AER-1.dta"))

# Use 0.1 bandwidth
df <- data_AER_1#filter( data_AER_1 , dzag != 0) #dzag >= -1 ,  , dzag <= 1 )

## Estimate different FE models to get outcome without FEs
# School transition score
agus_fe <- plm( agus ~ dga + dzag + dzag_after , data = df,
               index = c("uazY"), model = "within")
# Taking BA
bct_fe <- plm( bct ~ dga + dzag + dzag_after , data = df,
                    index = c("uazY"), model = "within")
# BA grade
bcg_fe <- plm( bcg ~ dga + dzag + dzag_after , data = df,
               index = c("uazY"), model = "within")


# Model output to check
summary(agus_fe)
summary(bct_fe)
summary(bcg_fe)

# Get the fixed effects and save them
f_agus <- fixef( agus_fe )
f_agus <- data.frame( uazY=names( f_agus ) , FE_agus = f_agus , row.names=NULL )
f_bct  <- fixef( bct_fe )
f_bct  <- data.frame( uazY=names( f_bct ) , FE_bct = f_bct , row.names=NULL )
f_bcg  <- fixef( bcg_fe )
f_bcg  <- data.frame( uazY=names( f_bcg ) , FE_bcg = f_bcg , row.names=NULL )

# Add FEs to the data
df <- left_join( df , f_agus , by = "uazY" )
df <- left_join( df , f_bct , by = "uazY" )
df <- left_join( df , f_bcg , by = "uazY" )

# Create outcome without FEs
df <- mutate( df , agus_dm = agus - FE_agus )
df <- mutate( df , bct_dm = bct - FE_bct )
df <- mutate( df , bcg_dm = bcg - FE_bcg )

# Checks
chck_1 <- lm( formula = agus_dm ~ dga + dzag + dzag_after +0 , data = df )
summary( chck_1 )
chck_2 <- lm( formula = bct_dm ~ dga + dzag + dzag_after +0 , data = df )
summary( chck_2 )
chck_3 <- lm( formula = bcg_dm ~ dga + dzag + dzag_after +0 , data = df )
summary( chck_3 )

## Save output
write_csv(df , paste0( path , "data/Pop_Elches/dataAER1_R_bw01.csv"))



# IV for BA grade with peer effects as well
#df <- mutate( df , uazYprsp = paste0(uazY,prsp))
#agus_fe_IV <- plm( agus ~ dga + dzag + dzag_after , data = df,
#                   index = c("uazYprsp"), model = "within")
#df <- mutate( df , agus_hat = predict(agus_fe_IV,data=df) )
              
