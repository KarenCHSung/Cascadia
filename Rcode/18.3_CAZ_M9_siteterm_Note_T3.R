#***********************************************************************************************#
#                                                                                               #
#  Varying coefficient model using INLA                                                         #
#                                                                                               #
#***********************************************************************************************#

#-----------------------------------------------------------------------------------------------#
#0.load packages--------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

  library(INLA)
  library(fields)
  library(viridisLite)
  library(gstat)
  library(sp)
  library(rgdal)

#-----------------------------------------------------------------------------------------------#    
#1. function definitions------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
  
#=== 1.1 plot figures in R ===  
 local.plot.field = function(field, mesh, xlim=c(0,10), ylim=c(0,10), ...){
        stopifnot(length(field) == mesh$n)
        proj = inla.mesh.projector(mesh, xlim = xlim, 
                             ylim = ylim, dims=c(300, 300))
        field.proj = inla.mesh.project(proj, field)
        n.col = 20
        image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
             xlim = xlim, ylim = ylim, col = plasma(n.col), nlevel=n.col+1)
  }

#=== here separate to three function 1.2 - 1.4 to avoid mixing======  
#=== 1.2 output the Lon ===
 local.plot.field2 = function(field, mesh, xlim=c(0,10), ylim=c(0,10), ...){
        stopifnot(length(field) == mesh$n)
        proj = inla.mesh.projector(mesh, xlim = xlim, 
                             ylim = ylim, dims=c(200, 200))
        return(proj$x)
 }
 
#=== 1.3 output the Lat ===
 local.plot.field3 = function(field, mesh, xlim=c(0,10), ylim=c(0,10), ...){
       stopifnot(length(field) == mesh$n)
       proj = inla.mesh.projector(mesh, xlim = xlim, 
                             ylim = ylim, dims=c(200, 200))
       return(proj$y)
 }

#=== 1.4 output the mean or std ===
 local.plot.field4 = function(field, mesh, xlim=c(0,10), ylim=c(0,10), ...){
       stopifnot(length(field) == mesh$n)
       proj = inla.mesh.projector(mesh, xlim = xlim, 
                             ylim = ylim, dims=c(200, 200))
       field.proj = inla.mesh.project(proj, field)
       return(field.proj)
 } 

 
#-----------------------------------------------------------------------------------------------#    
#2. set up the periods  ------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------# 
 
 pr =  c(1, 1.5, 2, 2.5, 3, 4, 5, 6, 7.5, 10)  

#-----------------------------------------------------------------------------------------------#   
#3.INLA      -----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------# 
#for(i in 1:length(pr)){
 
 i=5

#=== 3.1 read data ====  
 dat_stan <- read.csv(paste0("C:/Users/karen/OneDrive/орн▒/Berkeley/SUB_CAS/output/04_new_Z25/Residual&std_erg/UTM/Between_site_T",pr[i],"_SFirst19_allPeriod_UTM.csv",seq=""), header=T)
 names(dat_stan)
 
#=== 3.2 semi-variance ====
#--- use the semivariogram to pick up the correlation length ---#
#--- instead of using the INLA to estimate (save more time)  ---#
#--- Sta_lat/Sta_long: latitude/longitude for sites          ---#
#--- X.Intercept.: Between-site residual for ergodic model   ---#
#--- UTME/UTMN: UTM east/UTM north for sites                 ---# 
#--- fitVar$range[2]: correlation length                     ---#
 
 corrL = dat_stan[c("Sta_lat","Sta_long","X.Intercept.", "UTME","UTMN")]
 coordinates(corrL) = ~ UTMN +UTME 
 SemiVar = variogram(X.Intercept.~1, corrL)
 fitVar = fit.variogram(SemiVar, vgm("Sph"))
 
 
#=== 3.3 create the file (df.covar) ====
#--- Sta_lat/Sta_long: latitude/longitude for sites          ---#
#--- X.Intercept.: Between-site residual for ergodic model   ---#
#--- Sta_ID.1: site ID from 3D simulaiton                    ---# 
 coords_stat <- cbind(dat_stan$Sta_lat,dat_stan$Sta_long)
 statid <- dat_stan$Sta_ID.1
 resid <- dat_stan$X.Intercept.
 
 df.covar <- data.frame(intercept = 1,
                        stat = statid,
                        resid = resid)

#=== 3.4 spatial model for event and station terms terms ==== 
#3.4.1 use UTME/UTMN to build mesh, here unit is meter   
location <- cbind(dat_stan$UTME,dat_stan$UTMN)

#3.4.2 decide the maxedge (you can set up by yourself), here unit is km
exp.range = round(as.numeric(quantile(dist(location, upper=T, diag=F)/1000,0.1)),0)
MaxEdge = exp.range/5

#3.4.3 set up the mesh
#loc: location coordinates that are used as initial vertices
#max.edge: Maximum triangle edge length in the c(inner,outer) segment
#cutoff: Minimum distance between two distinct points (the unit should be equal to loc)
#offset: distance that specifies the size of the inner and outer extensions around the
#        data locations (here I turn off)
#can use "mesh1$n" to check the number of mesh
mesh1 = inla.mesh.2d(loc=location,
                     max.edge = c(1,5)*(MaxEdge*1000),
                     cutoff = (MaxEdge*1000)/5)#,
                     #offset = c(5 * max.edge, bound.outer))


#3.4.4 plot the mesh to check
#X1_s <- coords_stat[,2]
#X2_s <- coords_stat[,1]
#coords <- unique(cbind(c(X1_s), c(X2_s)))
#co_stat <- cbind(X1_s,X2_s)
#plot(mesh1, main="1st attempt");
#points(co_stat[,1], co_stat[,2], col="blue")
#axis(1); axis(2)


#=== 3.5 set up the prior and define model for stat terms  ===#
spde_stat <- inla.spde2.pcmatern(
  # Mesh and smoothness parameter(alpha)
  mesh = mesh1, alpha = 2,
  # P(practic.range < 0.3) = 0.5
  prior.range = c(10, 0.9),
  # P(sigma > 1) = 0.01
  prior.sigma = c(1, 0.01)) 

A_stat <- inla.spde.make.A(mesh1, loc = location)
idx.stat <- inla.spde.make.index("idx.stat",spde_stat$n.spde)


#=== 3.6 create the stack ===#
#--- information is from df.covar file ---#
stk1 <- inla.stack(
  data = list(y = resid),
  A = list(A_stat, 1),
  effects = list(idx.stat = idx.stat,
                 df.covar), tag = 'model_eqstat')

#=== 3.7 set up the formula ===#
#--- use the correlation length from semivariogram ---#
#--- must be use the log unit in INLA                 #
fix_rho = log(fitVar$range[2])

#--- if user dont want to fix the correlation length -#
#--- just need to remove hyper=list()                 #
form <- y ~ -1 +
  f(stat, model = "iid") +
  f(idx.stat, model = spde_stat, hyper=list(theta1=list(initial=fix_rho, fixed=TRUE)))

#=== 3.8 run inla ===#
#--- dic: the model deviance information criterion  ---#
#--- waic: the Watanabe-Akaike information criterion --#
#--- cpo: conditional predictive ordinates          ---#
fit_inla_spatial <- inla(form, 
                         data = inla.stack.data(stk1), 
                         control.predictor = list(A = inla.stack.A(stk1)),
                         control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE)
)


#=== 3.9 out put the hyper parameter ===#
hyper = fit_inla_spatial$summary.hyperpar


#=== 3.10 out put the map for UTME/UTMN ===#
ST_LON = local.plot.field2(fit_inla_spatial$summary.random[['idx.stat']][['mean']], mesh1, xlim = c(450000,650000), ylim = c(5150000,5400000))
ST_LAT = local.plot.field3(fit_inla_spatial$summary.random[['idx.stat']][['mean']], mesh1,xlim = c(450000,650000), ylim = c(5150000,5400000))

#=== 3.11 out put the prediction for nonergodic site ===#
ST_MEAN = local.plot.field4(fit_inla_spatial$summary.random[['idx.stat']][['mean']], mesh1, xlim = c(450000,650000), ylim = c(5150000,5400000))

#=== 3.12 out put the epistemic for nonergodic site ===#
ST_SD = local.plot.field4(fit_inla_spatial$summary.random[['idx.stat']][['sd']], mesh1,  xlim = c(450000,650000), ylim = c(5150000,5400000))

#=== 3.13 out put the file ===#
ST_infor = matrix(0,nrow = length(ST_LON)*length(ST_LAT), ncol= 4)

for(m in 1:length(ST_LON)){
  for(f in 1:length(ST_LAT)){
    
    Poo = length(ST_LON)*(m-1) + f
    ST_infor[Poo,1] =  ST_LON[m]
    ST_infor[Poo,2] =  ST_LAT[f]
    ST_infor[Poo,3] =  ST_MEAN[m,f]
    ST_infor[Poo,4] =  ST_SD[m,f]
  }}

path4 <- paste0("C:/Users/karen/OneDrive/орн▒/Berkeley/SUB_CAS/output/06_nonerg_adj/STA_const_T",pr[i],"_mu&std_nointer19_allPeriod_UTM3.csv",seq="")
write.table(ST_infor, file= path4,sep = ",",row.names = F,col.names=c("UTME","UTMN","B0","sd"))

#-----------------------------------------------------------------------------------------------#   
# 4. put the non-ergodic site term into the original file --------------------------------------#
#-----------------------------------------------------------------------------------------------# 

#nonergodic term
for (z in 1 : dim(dat_stan)[1]){

  #adj_STconst
  del_STLon =  dat_stan$UTME[z] - ST_LON
  LonST_ID = which.min(abs(del_STLon))
  dat_stan$adj_STlon[z] = ST_LON[LonST_ID]
  
  del_STLat =  dat_stan$UTMN[z] - ST_LAT
  LatST_ID = which.min(abs(del_STLat))
  dat_stan$adj_STlat[z] = ST_LAT[LatST_ID]
  
  dat_stan$adj_STconst[z] = ST_MEAN[LonST_ID,LatST_ID]
  
}


#New Residual
sub_data = dat_stan

  sub_data$newResid = sub_data$X.Intercept. - (sub_data$adj_STconst)
  sub_data$SemiVar_range = fitVar$range[2]
  sub_data$SemiVar_sill = fitVar$psill[2]


#plot Z2.5 versus Residual (check results)
windows(width = 6, height = 6,bg="white")
par(mar=c(3.5,5,1,1))
plot(sub_data$Rrup,sub_data$T_0.2 ,type="n",log="x",ylim=c(-2.5,2.5),xlim=c(2,9000),cex=0.5,pch=24,ann=FALSE ,axes=FALSE,cex.lab=2.3) 
x.at <- c(0.1,1,5,10,20,50,100,200,300,400,500,1000,10000)
y.at <- c(seq(-4,4,by=0.5))
axis(1, at=x.at,las=1, labels = formatC(x.at, format="fg"), cex.axis=1.4,lwd=2,tck=0.01,line=-0.4,col="white")
axis(2, at=y.at,las=2, labels = formatC(y.at, format="fg"), cex.axis=1.4,lwd=2,tck=0.01,line=-0.4,col="white")
tic<-c(seq(-4,4,by=0.5))
tic2<-c(seq(0.1,0.9,by=0.1),seq(1,10,by=1),seq(10,100,by=10),seq(100,1000,by=100),seq(1000,10000,by=1000))
abline(h=tic, col="grey")
abline(v=tic2, col="grey")
mtext("       Z2.5", side=1, line=-1.5,outer = T,cex=1.8)
mtext("     Between-site residual (ln unit)", side=2, line=-1.5,outer = T,cex=1.8)
axis(3,at=x.at, labels = F,lwd=1.5,tck=0)
axis(4,at=y.at, labels = F,lwd=1.5,tck=0)   
axis(1,at=x.at, labels = F,lwd=1.5,tck=0.01)
axis(2,at=y.at, labels = F,lwd=1.5,tck=0.01) 

points(sub_data$Z25,sub_data$X.Intercept. , cex=1)
points(sub_data$Z25,sub_data$newResid , cex=1, col="red")
text(50,-2,paste("T=",pr[i],sep=""),cex=2)
dev.off()

print(i)

#}
