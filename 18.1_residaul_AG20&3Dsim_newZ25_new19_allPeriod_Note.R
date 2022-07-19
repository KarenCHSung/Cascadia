#**************************************************************************************#
#                                                                                      #
# 2021/1/22                                                                            #
# calculate the total residual from AG20 model                                         #
# dataset is from the 3D simulation                                                    #
# 2021/03/06                                                                           #
# turn off the Basin-Depth Scaling term                                                #
# 2021/03/21                                                                           #
# add the new Basin-Depth dscaling term                                                #
#                                                                                      #
#**************************************************************************************#

#--------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------#
#0. read data -------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------#

#===read 3D simulation ====
datapath3 <- paste("C:/Users/karen/OneDrive/桌面/Berkeley/SUB_CAS/data/final_data/M9_SA_data_fianl.csv",seq="")
sub_data <- read.table(datapath3, stringsAsFactors=FALSE, header=T, sep=",")
dim(sub_data)
names(sub_data)

#===read coeff for AG20 ===
datapath3 <- paste("C:/Users/karen/OneDrive/桌面/Berkeley/SUB_CAS/data/Coeff_AG20/AG20_coeff_with_inter_allPeriod.csv",seq="")
coeff_AG20 <- read.table(datapath3, stringsAsFactors=FALSE, header=T, sep=",")
names(coeff_AG20)


#--------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------#
#1. set up the adjusted AG20 model (Ergodic model for Cascadia) -----------------------#
#--------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------#

AG20_fun = function(M,R,Fualt,Vs30,Z25,Ztor,FAS,i,PGA1000,CASS){
  
  a6 = coeff_AG20$a6[i] 
  a12 = coeff_AG20$a12[i] 
  a1 = coeff_AG20$a1[i] 
  C4 = coeff_AG20$C4[i]
  a9 = coeff_AG20$a9[i]
  C1i = coeff_AG20$c1i[i]
  a4 = coeff_AG20$a4[i]
  a45 = coeff_AG20$a45[i]
  a13 = coeff_AG20$a13[i]
  a5 = coeff_AG20$a5[i]
  a8 = coeff_AG20$a8[i]
  a11 = coeff_AG20$a11[i]
  Vlin = coeff_AG20$Vlin[i]
  b = coeff_AG20$b[i]
  c = 1.88
  n = 1.18
  a39 =  coeff_AG20$a39[i] 
  a15 =  coeff_AG20$a15[i]
  a25 =  coeff_AG20$a25[i] 
  a18 =  coeff_AG20$a18[i] 
  a32 =  coeff_AG20$a32[i]
  a10 =  coeff_AG20$a10[i]
  a14 =  coeff_AG20$a14[i]
  adj_CAS =  coeff_AG20$Cascadia[i]
  a2 = coeff_AG20$a2[i]
  a3 = coeff_AG20$a3[i]
  
  #==== new coeff for Basin-Depth scaling ===
  a39_b22 = coeff_AG20$CSIM[i] #CSIM
  a39_b1 = coeff_AG20$b1[i] #b1
  a39_b2 = coeff_AG20$b2[i] #b2
  
  #Regionalization of the Base Model for Cascadia (CASS == 1)
  if(CASS == 1){  
    #large distance (linear R) scaling 
    a_6 = a6 + a25
    #linear site-amplification scaling
    a_12 = a12 + a18
    #the constant term
    a_1 =  a32
  }
  
  #no Regionalization for the Base Model (CASS == 0)
  #if need to use into other regions, still need to add another Regionalization
  if(CASS == 0){  
    #large distance (linear R) scaling
    a_6 = a6 
    #linear site-amplification scaling
    a_12 = a12 
    #the constant term
    a_1 =  a1 
  }
  
  #finite-fault term
  HFF = C4 * exp(a9 * (M - 6))
  
  #Magnitude Scaling
  C1s = 7.1
  C1 = Fualt * C1s + (1 - Fualt) * C1i
  
  if(M <= C1){
    fmag = (a4 + Fualt * a45) + a13 * (10 - M)^2
  }
  
  if(M > C1){
    fmag = a5 * (M - C1) + a13 * (10 - M)^2
  }
  
  #Intraslab Scaling
  deltaC1s = C1s - 7.5
  fslab = a10 + a4 * (C1s - deltaC1s - 7.5) + a14 * log(R + HFF)
  
  #Depth Scaling
  if(Ztor <= 50){
    fztor = a8 * (Ztor - 50) * Fualt
  }
  if(Ztor > 50 & Ztor <= 200){
    fztor = a11 * (Ztor - 50) * Fualt
  }
  if(Ztor > 200){
    fztor = a11 * (150) * Fualt
  }
  
  #Site-Response Scaling
  if(Vs30 > 1500) {Vs = 1500}
  if(Vs30 <= 1500) {Vs = Vs30}
  
  if(Vs30 < Vlin){
    fsite = a_12 * log(Vs/Vlin) - b * log(PGA1000 + c) + b * log(PGA1000 + c * (Vs/Vlin)^n)
  }
  if(Vs30 >= Vlin){
    fsite = (a_12 + b * n) * log(Vs/Vlin)
  }
  
  #Basin-Depth Scaling 
  if(Vs30 <= 200){lnZ25_CAS = 8.52}
  if(Vs30 > 200 & Vs30 < 570) {lnZ25_CAS = 8.52 - 0.88 * log(Vs30/200)}
  if(Vs30 >=  570) {lnZ25_CAS = 7.6}   
  Z25_ref_CAS = exp(lnZ25_CAS)
  Z_25 = (Z25 + 50) / (Z25_ref_CAS + 50)
  
  #the basin term from 3D simulation
  Zx_AG20 = 1
  if(Z_25 < Zx){ fbasin_CAS = b2 } #b2
  if(Z_25 >= Zx) fbasin_CAS = b22 + b1 * log(Z_25/Zx_AG20) #b1*log(Z_25/Zx) + b22
  
  #original basin term for AG20 
  #if(log(Z_25) <= -0.58){ fbasin_CAS = 0.0 }
  #if(log(Z_25) > -0.58) fbasin_CAS = fbasin_CAS = a39 * log(Z_25)

  
  #adjustment for CAS
  adj_CAS = adj_CAS
  
  #PSA--------------------------------------------------------------------------------------------------------------------------------- 
  #Cascadia
  if(CASS == 1){
    #intra
    if(Fualt == 1){ 
      lnPSA =  a_1 + (a2 + a3*(M - 7))*log(R + HFF) + a_6*R + fmag  + Fualt*fztor + fsite + Fualt*fslab + fbasin_CAS  + a15*FAS  + adj_CAS
    }
    #interface
    if(Fualt == 0){ 
      lnPSA =  a_1 + (a2 + a3*(M - 7))*log(R + HFF) + a_6*R + fmag  + Fualt*fztor + fsite + Fualt*fslab + fbasin_CAS   + a15*FAS  + adj_CAS
    }
  }
  
  #not Cascadia
  if(CASS == 0){
    #intra
    if(Fualt == 1){ 
      lnPSA =  a_1 + (a2 + a3 * (M - 7)) * log(R + HFF) + a_6 * R + fmag + fsite + Fualt*fztor + fbasin_CAS + Fualt*fslab   + a15 * FAS 
    }
    #interface
    if(Fualt == 0){ 
      lnPSA =  a_1 + (a2 + a3 * (M - 7)) * log(R + HFF) + a_6 * R + fmag + fsite + Fualt*fztor + fbasin_CAS + Fualt*fslab + a15 * FAS 
    }
    
  }   
  
  PSA = exp(lnPSA) 
  PSA
}



#--------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------#
#2. calculate the total residuals                               -----------------------#
#--------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------#
#==== 30 sources ====
for(m in 1:30){

 #chose the data from one source
  sub_data2 = sub_data[sub_data$EQ_ID == m,]


 #calculate the PSA from AG20 using the parameter from 3D-simulation 
 M=9
 Fualt=0
 Vs30=600
 FAS=0
 Ztor=0
 CASS=1
 
 for(j in 1:dim(sub_data2)[1]){
     PGA1000 = AG20_fun(M,sub_data2$Rrup[j],Fualt,1000,sub_data2$Z25[j],Ztor,FAS,1,0,0)
     sub_data2$AG20_T_0.01[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,1,PGA1000,CASS)
     sub_data2$AG20_T_0.02[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,2,PGA1000,CASS)
     sub_data2$AG20_T_0.03[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,3,PGA1000,CASS)
     sub_data2$AG20_T_0.05[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,4,PGA1000,CASS)
     sub_data2$AG20_T_0.075[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,5,PGA1000,CASS)
     sub_data2$AG20_T_0.1[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,6,PGA1000,CASS)
     sub_data2$AG20_T_0.15[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,7,PGA1000,CASS)
     sub_data2$AG20_T_0.2[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,8,PGA1000,CASS)
     sub_data2$AG20_T_0.25[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,9,PGA1000,CASS)
     sub_data2$AG20_T_0.3[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,10,PGA1000,CASS)
     sub_data2$AG20_T_0.4[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,11,PGA1000,CASS)
     sub_data2$AG20_T_0.5[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,12,PGA1000,CASS)
     sub_data2$AG20_T_0.6[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,13,PGA1000,CASS)
     sub_data2$AG20_T_0.75[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,14,PGA1000,CASS)
     sub_data2$AG20_T_1[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,15,PGA1000,CASS)
     sub_data2$AG20_T_1.5[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,16,PGA1000,CASS)
     sub_data2$AG20_T_2[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,17,PGA1000,CASS)
     sub_data2$AG20_T_2.5[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,18,PGA1000,CASS)
     sub_data2$AG20_T_3[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,19,PGA1000,CASS)
     sub_data2$AG20_T_4[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,20,PGA1000,CASS)
     sub_data2$AG20_T_5[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,21,PGA1000,CASS)
     sub_data2$AG20_T_6[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,22,PGA1000,CASS)
     sub_data2$AG20_T_7.5[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,23,PGA1000,CASS)
     sub_data2$AG20_T_10[j] = AG20_fun(M,sub_data2$Rrup[j],Fualt,Vs30,sub_data2$Z25[j],Ztor,FAS,24,PGA1000,CASS)
     
     sub_data2$resid_T0.01[j] = log(sub_data2$T_0.01[j]) - log(sub_data2$AG20_T_0.01[j])
     sub_data2$resid_T0.02[j] = log(sub_data2$T_0.02[j]) - log(sub_data2$AG20_T_0.02[j])
     sub_data2$resid_T0.03[j] = log(sub_data2$T_0.03[j]) - log(sub_data2$AG20_T_0.03[j])
     sub_data2$resid_T0.05[j] = log(sub_data2$T_0.05[j]) - log(sub_data2$AG20_T_0.05[j])
     sub_data2$resid_T0.075[j] = log(sub_data2$T_0.75[j]) - log(sub_data2$AG20_T_0.075[j])
     sub_data2$resid_T0.1[j] = log(sub_data2$T_0.1[j]) - log(sub_data2$AG20_T_0.1[j])
     sub_data2$resid_T0.15[j] = log(sub_data2$T_0.15[j]) - log(sub_data2$AG20_T_0.15[j])
     sub_data2$resid_T0.2[j] = log(sub_data2$T_0.2[j]) - log(sub_data2$AG20_T_0.2[j])
     sub_data2$resid_T0.25[j] = log(sub_data2$T_0.25[j]) - log(sub_data2$AG20_T_0.25[j])
     sub_data2$resid_T0.3[j] = log(sub_data2$T_0.3[j]) - log(sub_data2$AG20_T_0.3[j])
     sub_data2$resid_T0.4[j] = log(sub_data2$T_0.4[j]) - log(sub_data2$AG20_T_0.4[j])
     sub_data2$resid_T0.5[j] = log(sub_data2$T_0.5[j]) - log(sub_data2$AG20_T_0.5[j])
     sub_data2$resid_T0.6[j] = log(sub_data2$T_0.6[j]) - log(sub_data2$AG20_T_0.6[j])
     sub_data2$resid_T0.7[j] = log(sub_data2$T_0.7[j]) - log(sub_data2$AG20_T_0.75[j])
     sub_data2$resid_T1[j] = log(sub_data2$T_1[j]) - log(sub_data2$AG20_T_1[j])
     sub_data2$resid_T1.5[j] = log(sub_data2$T_1.5[j]) - log(sub_data2$AG20_T_1.5[j])
     sub_data2$resid_T2[j] = log(sub_data2$T_2[j]) - log(sub_data2$AG20_T_2[j])
     sub_data2$resid_T2.5[j] = log(sub_data2$T_2.5[j]) - log(sub_data2$AG20_T_2.5[j])
     sub_data2$resid_T3[j] = log(sub_data2$T_3[j]) - log(sub_data2$AG20_T_3[j])
     sub_data2$resid_T4[j] = log(sub_data2$T_4[j]) - log(sub_data2$AG20_T_4[j])
     sub_data2$resid_T5[j] = log(sub_data2$T_5[j]) - log(sub_data2$AG20_T_5[j])
     sub_data2$resid_T6[j] = log(sub_data2$T_6[j]) - log(sub_data2$AG20_T_6[j])
     sub_data2$resid_T7.5[j] = log(sub_data2$T_7.5[j]) - log(sub_data2$AG20_T_7.5[j])
     sub_data2$resid_T10[j] = log(sub_data2$T_10[j]) - log(sub_data2$AG20_T_10[j])
     #print(j)
 }
 

 path3 <- paste0("C:/Users/karen/OneDrive/桌面/Berkeley/SUB_CAS/output/03_Residual_GMM_newZ2.5/M9_SA_data_Residual_EQID",m,"_new19_allPeriod.csv",seq="")
 write.table(sub_data2, file= path3,sep = ",",row.names = F)
 
 print(m)
}
 

#--------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------#
# Combine the files--------------------------------------------------------------------#
#--------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------#

#5.read EQID
EQ <- read.csv(paste0("C:/Users/karen/OneDrive/桌面/Berkeley/SUB_CAS/data/source/source_data_detail_latlon.csv",seq=""), header=T)
names(EQ)

EQ_NORTH2 = EQ

#6. combine EQs
u = 1
file0 =  read.csv(paste0("C:/Users/karen/OneDrive/桌面/Berkeley/SUB_CAS/output/03_Residual_GMM_newZ2.5/M9_SA_data_Residual_EQID",u,"_new19_allPeriod.csv",seq=""), header=T,sep=",")
dim(file0)
#file0$NEQ = EQ_NORTH2$EQID[1]
file0$EQ_lon = EQ_NORTH2$lon[1]
file0$EQ_lat = EQ_NORTH2$lat[1]

for(j in 2:dim(EQ_NORTH2)[1]){
  
  u = j
  
  file2 =  read.csv(paste0("C:/Users/karen/OneDrive/桌面/Berkeley/SUB_CAS/output/03_Residual_GMM_newZ2.5/M9_SA_data_Residual_EQID",u,"_new19_allPeriod.csv",seq=""), header=T)
  #file2$NEQ = EQ_NORTH2$EQID[j]
  file2$EQ_lon = EQ_NORTH2$lon[j]
  file2$EQ_lat = EQ_NORTH2$lat[j]
  
  file0 = rbind(file0,file2)
  
  print(dim(file0))
  print(j)
  
}

#7. output
path4 <- paste0("C:/Users/karen/OneDrive/桌面/Berkeley/SUB_CAS/output/03_Residual_GMM_newZ2.5/M9_SA_data_Resid_All_w.o._Z2.5_new19_allPeriod.csv",seq="")
write.table(file0, file= path4,sep = ",",row.names = F)



#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#
#---calculate between-site residuals for ergodic -------------------------------#
#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#

library(nlme)

#0.set up the period
#pr =  c(1, 1.5, 2, 2.5, 3, 4, 5, 7.5, 10)
pr =  c(0.01, 0.02, 0.03, 0.05, 0.075, 0.1,
        0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75,
        1, 1.5, 2, 2.5, 3, 4, 5, 6, 7.5, 10)

#8.read data: sampled recordings
datapath3 <- paste0("C:/Users/karen/OneDrive/桌面/Berkeley/SUB_CAS/output/03_Residual_GMM_newZ2.5/M9_SA_data_Resid_All_w.o._Z2.5_new19_allPeriod.csv",seq="")
sub_data <- read.table(datapath3, stringsAsFactors=FALSE, header=T, sep=",")
dim(sub_data)
names(sub_data)

for(i in 1:length(pr)[1]){
  
  #total residual
  sub_data$CAZ_Tstd_erg = round(sd(sub_data[,(59+i)]),digits = 3)
 
  #single-site residual (partail)
  TotalResid = sub_data[,(59+i)]
  CAZ_LME.site_pt<-lme(TotalResid~ 1, random=~1|Sta_ID, data=sub_data,control=lmeControl(opt = "optim"))
  CAZ_Ranef.site_pt <-ranef(CAZ_LME.site_pt,augFrame=T,which=c("Sta_ID","Sta_lat","Sta_long","Sta_ID","Z25"))
  InterSitepath <-paste0("C:/Users/karen/OneDrive/桌面/Berkeley/SUB_CAS/output/04_new_Z25/Residual&std_erg/Between_site_T",pr[i],"_SFirst19_allPeriod.csv",seq="")
  write.table(CAZ_Ranef.site_pt, file= InterSitepath,sep = ",",row.names = F)
  
  print(CAZ_LME.site_pt$coefficients$fixed)
  
  cc = VarCorr(CAZ_LME.site_pt)
  InterSitepath <- paste0("C:/Users/karen/OneDrive/桌面/Berkeley/SUB_CAS/output/04_new_Z25/Residual&std_erg/Between_site_T",pr[i],"_SFirst19_allPeriod.txt",seq="")
  write.table(cc, file= InterSitepath,sep = ",",row.names = F)
  
  
  CAZ_pt_RSD <-residuals(CAZ_LME.site_pt,level=1)
  sub_data$CAZ_pt_RSD <- CAZ_pt_RSD
  #sub_data$Ptstd_pt <- round(sd(sub_data$CAZ_pt_RSD),digits = 3)
  
  
  #output-------------------------------------------------------------------------------------------------------------------
  new.file = matrix(0, ncol=16, nrow=dim(sub_data)[1])
  
  #genertal info
  new.file[,1] = sub_data$sim_EQ_ID
  new.file[,2] = sub_data$sim_Sta_ID
  new.file[,3] = sub_data$Sta_lat
  new.file[,4] = sub_data$Sta_long
  new.file[,5] = sub_data$Rrup
  new.file[,6] = sub_data$Z25
  new.file[,7] = sub_data[,(59+i)]
  new.file[,8] = sub_data$Sta_ID
  new.file[,9] = sub_data$EQ_ID
  new.file[,10] = sub_data$EQ_lon
  new.file[,11] = sub_data$EQ_lat
  #resid&std_ergodic
  new.file[,12] = sub_data$CAZ_Tstd_erg
  #new.file[,13] = sub_data$CAZ_WE_RSD_erg
  #new.file[,14] = sub_data$CAZ_WS_RSD_erg
  new.file[,13] = sub_data$CAZ_pt_RSD
  new.file[,14] = sub_data[,(11+i)]
  new.file[,15] = sub_data[,(35+i)]
  new.file[,16] = CAZ_LME.site_pt$coefficients$fixed[1]
  
  path3 <- paste0("C:/Users/karen/OneDrive/桌面/Berkeley/SUB_CAS/output/04_new_Z25/Residual&std_erg/resid&std_T",pr[i],"_erg19_allPeriod.csv",seq="")
  write.table(new.file, file= path3,sep = ",",row.names = F,col.names = c("sim_EQ_ID","sim_Sta_ID","Sta_lat","Sta_long","Rrup","Z25",
                                                                          "Resid", "NSTA","NEQ","EQ_lon","EQ_lat",
                                                                          "Total_std_erg",#"intra_RSD_erg","record_RSD_erg",
                                                                          "singleST_RSD_erg",
                                                                          "sim_PSA","AG20_PSA","intercept"))
  print(i)
  
}




 