rm(list=ls())
library(tidyverse)
library(rdrobust)


path <- '/Users/areguly6/Library/CloudStorage/OneDrive-CorvinusUniversityofBudapest/Research/RDD/develop/'
df = read_csv( paste0( path, 'data/Pop_Elches/dataAER1_R_bw01.csv' ) )
df <- df %>% 
  filter( dzag >= -1, dzag <= 1, dzag != 0)

df_bt <- df %>% 
  filter( zga < 6.77 )
df_ut <- df %>% 
  filter( zga > 7.74 )
df_three <- df %>% 
  filter( nusua == 3 )
df_two <- df %>% 
  filter( nusua == 2 )
df_four <- df %>% 
  filter( nusua > 3 )
# School avg transition scores
#res_all_a <- rdrobust( df$agus_dm, df$dzag, c = 0, cluster = df$sid2 )
res_all_b <- rdrobust( df$bct_dm, df$dzag, c = 0, cluster = df$sid2 )
res_all_c <- rdrobust( df$bcg_dm, df$dzag, c = 0, cluster = df$sid2 )
#summary(res_all)

#res_bottom_a <- rdrobust( df_bt$agus_dm, df_bt$dzag, c = 0, cluster = df_bt$sid2 )
res_bottom_b <- rdrobust( df_bt$bct_dm, df_bt$dzag, c = 0, cluster = df_bt$sid2 )
res_bottom_c <- rdrobust( df_bt$bcg_dm, df_bt$dzag, c = 0, cluster = df_bt$sid2 )
#summary(res_bottom)

#res_upper_a <- rdrobust( df_ut$agus_dm, df_ut$dzag, c = 0, cluster = df_ut$sid2 )
res_upper_b <- rdrobust( df_ut$bct_dm, df_ut$dzag, c = 0, cluster = df_ut$sid2 )
res_upper_c <- rdrobust( df_ut$bcg_dm, df_ut$dzag, c = 0, cluster = df_ut$sid2 )
#summary(res_upper_a)

# Number of schools in town

# Two schools
res_two_a <- rdrobust( df_two$agus_dm, df_two$dzag, c = 0, cluster = df_two$sid2 )
res_two_b <- rdrobust( df_two$bct_dm, df_two$dzag, c = 0, cluster = df_two$sid2 )
res_two_c <- rdrobust( df_two$bcg_dm, df_two$dzag, c = 0, cluster = df_two$sid2 )
#summary(res_two_a)

# Three schools
res_three_a <- rdrobust( df_three$agus_dm, df_three$dzag, c = 0, cluster = df_three$sid2 )
res_three_b <- rdrobust( df_three$bct_dm, df_three$dzag, c = 0, cluster = df_three$sid2 )
res_three_c <- rdrobust( df_three$bcg_dm, df_three$dzag, c = 0, cluster = df_three$sid2 )
#summary(res_two_a)

# Four or more
res_four_a <- rdrobust( df_four$agus_dm, df_four$dzag, c = 0, cluster = df_four$sid2 )
res_four_b <- rdrobust( df_four$bct_dm, df_four$dzag, c = 0, cluster = df_four$sid2 )
res_four_c <- rdrobust( df_four$bcg_dm, df_four$dzag, c = 0, cluster = df_four$sid2 )
#summary(res_two)