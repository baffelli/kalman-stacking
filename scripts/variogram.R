library(geoR)
library(akima)
#Read data
unw_stable_file<-read.table('/home/baffelli/PhD/trunk/Bisgletscher_pol/ipt/20150803_120249_AAAl_20150803_120519_AAAl.punw_stable.csv', sep=',', head=T)
unw_grid_file<-read.table('/home/baffelli/PhD/trunk/Bisgletscher_pol/ipt/20150803_120249_AAAl_20150803_120519_AAAl.punw_grid.csv', sep=',', head=T)

#Read heights for corresponding points
hgt_stable_file<-read.table('/home/baffelli/PhD/trunk/Bisgletscher_pol/ipt/20150803_120249_AAAl_20150803_120519_AAAl.phgt_stable.csv', sep=',', head=T)
hgt_grid_file<-read.table('/home/baffelli/PhD/trunk/Bisgletscher_pol/ipt/20150803_120249_AAAl_20150803_120519_AAAl.phgt_grid.csv', sep=',', head=T)


#Create convenient data frames
hgt_grid = data.frame(x=hgt_grid_file$x, y=hgt_grid_file$y,  r=hgt_grid_file$ridx, z=hgt_grid_file$record_0)
hgt_stable = data.frame(x=hgt_stable_file$x, y=hgt_stable_file$y, r=hgt_stable_file$ridx, az=hgt_stable_file$azidx, z=hgt_stable_file$record_0)
ifgram_stable = data.frame(x=unw_stable_file$x, y=unw_stable_file$y, r=unw_stable_file$ridx, az=unw_stable_file$azidx, phase=unw_stable_file$record_0, z=hgt_stable$z)
ifgram_grid = data.frame(x=unw_grid_file$x, y=unw_grid_file$y, r=unw_grid_file$ridx, az=unw_grid_file$azidx, phase=unw_grid_file$record_0, z=hgt_grid$z)




#Color ramp for plotting
cols <-  colorRampPalette(c("blue", 'white', "red"))(256)

#Model as a function of location
model_loc_eq <- phase ~ x + y + I(x ** 2) + I(y ** 2)  + I(y * x)
model_loc <- lm(model_loc_eq, data = ifgram_stable)
predicted_aps_stable.first_fit = predict(model_loc, se.fit = TRUE)

#Model as a function of height
model_h_eq <- phase ~ z + I(z ** 2)
model_h <- lm(model_h_eq, data = ifgram_stable)
predicted_aps_stable.second_fit = predict(model_h, se.fit = TRUE)

#Mixed models as a function of height and location
model_mixed_eq <-
  phase ~ z + I(z ** 2)  + x + y + I(x ** 2) + I(y ** 2)  + I(y * x)
model_mixed <- lm(model_mixed_eq, data = ifgram_stable)
predicted_aps_stable.third_fit = predict(model_mixed, se.fit = TRUE)

#Mixed models as a function of range and azimuth
model_mixed_pol_eq <- phase ~ r + I(r ** 2)
model_mixed_pol <- lm(model_mixed_pol_eq, data = ifgram_stable)
predicted_aps_stable.fourth_fit = predict(model_mixed_pol, se.fit = TRUE)

#Function to scale colors
scale_colors <-
  function(data,
           color_gradient ,
           lmin = -2 * pi,
           lmax = 2 * pi) {
    color_gradient[cut(data, breaks = seq(lmin, lmax, l = length(color_gradient)))]
  }

par(mfrow = c(4, 3))
lmin <- -10
lmax <- 10
sf <- 2
pch <- 1
r <- ifgram_stable$r
az <- ifgram_stable$az
#Unwrapped phase
plot(
  r,
  az,
  col = scale_colors(ifgram_stable$phase, cols , lmin = lmin, lmax = lmax),
  pch = pch,
  main = "Unwrapped point scatterers phase"
)
#First model
plot(
  r,
  az,
  col = scale_colors(
    predicted_aps_stable.first_fit$fit,
    cols,
    lmin = lmin,
    lmax = lmax
  ),
  pch = pch,
  main = "Modeled phase (location)",
  cex = predicted_aps_stable.first_fit$se.fit * sf
)
#First differential phase
plot(
  r,
  az,
  col = scale_colors(
    ifgram_stable$phase - predicted_aps_stable.first_fit$fit,
    cols,
    lmin = lmin,
    lmax = lmax
  ),
  pch = pch,
  main = "Differential phase"
)
#Phase again
plot(
  r,
  az,
  col = scale_colors(ifgram_stable$phase, cols , lmin = lmin, lmax = lmax),
  pch = pch,
  main = "Unwrapped point scatterers phase"
)
#Second model
plot(
  r,
  az,
  col = scale_colors(
    predicted_aps_stable.second_fit$fit,
    cols,
    lmin = lmin,
    lmax = lmax
  ),
  pch = pch,
  main = "Modeled phase (height)",
  cex = predicted_aps_stable.second_fit$se.fit * sf
)
#Second differential phase
plot(
  r,
  az,
  col = scale_colors(
    ifgram_stable$phase - predicted_aps_stable.second_fit$fit,
    cols,
    lmin = lmin,
    lmax = lmax
  ),
  pch = pch,
  main = "Differential phase"
)
#Phase again
plot(
  r,
  az,
  col = scale_colors(ifgram_stable$phase, cols , lmin = lmin, lmax = lmax),
  pch = pch,
  main = "Unwrapped point scatterers phase"
)
#Third model
plot(
  r,
  az,
  col = scale_colors(
    predicted_aps_stable.third_fit$fit,
    cols,
    lmin = lmin,
    lmax = lmax
  ),
  pch = pch,
  main = "Modeled phase (mixed)",
  cex = predicted_aps_stable.third_fit$se.fit * sf
)
#Third differential phase
plot(
  r,
  az,
  col = scale_colors(
    ifgram_stable$phase - predicted_aps_stable.third_fit$fit,
    cols,
    lmin = lmin,
    lmax = lmax
  ),
  pch = pch,
  main = "Differential phase"
)
#Phase again
plot(
  r,
  az,
  col = scale_colors(ifgram_stable$phase, cols , lmin = lmin, lmax = lmax),
  pch = pch,
  main = "Unwrapped point scatterers phase"
)
#Third model
plot(
  r,
  az,
  col = scale_colors(
    predicted_aps_stable.fourth_fit$fit,
    cols,
    lmin = lmin,
    lmax = lmax
  ),
  pch = pch,
  main = "Modeled phase (polar)",
  cex = predicted_aps_stable.third_fit$se.fit * sf
)
#Third differential phase
plot(
  r,
  az,
  col = scale_colors(
    ifgram_stable$phase - predicted_aps_stable.fourth_fit$fit,
    cols,
    lmin = lmin,
    lmax = lmax
  ),
  pch = pch,
  main = "Differential phase"
)

#Diff ifgram
res<-data.frame(x=ifgram$x,y=ifgram$y, phase=ifgram$phase - predicted_hgt$fit)

#Create dataframe
unw.geo <- as.geodata(ifgram, coords.col = 1:2, data.col = 3)
res.geo <- as.geodata(res, coords.col = 1:2)
#Compute variogram
lags = seq(30,5000, l=50)
unw.var <- variog(unw.geo, breaks=lags)
res.var <- variog(res.geo, breaks=lags)



#Fit model
unw.var.fit <- variofit(unw.var, ini.cov.pars = c(0.5,1000),cov.model ="exponential",fix.nugget=FALSE, nugget=0.5)
#Fit model
res.var.fit <- variofit(res.var, ini.cov.pars = c(0.5,2000),cov.model ="exponential",fix.nugget=FALSE, nugget=0.5)

#Plot models
plot(unw.var,ylim=c(0,2))
lines(res.var, col='red')
lines.variomodel(unw.var.fit)
lines.variomodel(res.var.fit,col='red')


#Create trend
aps_trend <-trend.spatial(model_eq,aps.geo)
#Kriging object
k <- krige.control(cov.model=res.var.fit$cov.model, cov.pars=res.var.fit$cov.pars, type.krige = 'OK', beta=0.2, dist.epsilon = 20)


#Locations of prediction
glacier_pg <- data.frame(x=glacier_hgt$x, y=glacier_hgt$y)

#Krige 
aps.kriged <- krige.conv(res.geo,loc=glacier_pg,krige=k)


#Plot variogram and fit
plot(unw.var)
lines(res.var, col='red')
lines.variomodel(unw.var.fit)




plot(res.geo$coords, col=cols[cut(res.geo$data, breaks = seq(-6,6,l=256))])
points(glacier_pg, col=cols[cut(aps.kriged$predict, breaks = seq(-6,6,l=256))], pch=0)

#Plot trend with height
plot(hgt$record_0, aps.geo$data)
