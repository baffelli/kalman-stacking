library(geoR)
library(akima)
#Read data
unw<-read.table('/home/baffelli/PhD/trunk/Bisgletscher_pol/ipt/20150803_060519_stack_20_AAAl.aps.csv', sep=',', head=T)
#Read heights for corresponding points
hgt_stable<-read.table('/home/baffelli/PhD/trunk/Bisgletscher_pol/ipt/20150803_060519_stack_20_AAAl.hgt_stable.csv', sep=',', head=T)
hgt_grid<-read.table('/home/baffelli/PhD/trunk/Bisgletscher_pol/ipt/20150803_060519_stack_20_AAAl.hgt_grid.csv', sep=',', head=T)


grd_hgt = data.frame(x=hgt_grid$x, y=hgt_grid$y, z=hgt_grid$record_0)
pt_hgt = data.frame(x=hgt_stable$x, y=hgt_stable$y, z=hgt_stable$record_0)
ifgram = data.frame(x=unw$x, y=unw$y, phase=unw$record_5, z=pt_hgt$z)




#Color ramp for plotting
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
#Model as a function of elevation
model_eq <- phase ~ x + y + I(x^2) + I(y^2) + I(x*y)  + I(z) + I(z^2)
model <- lm(model_eq, data=ifgram)
predicted_hgt = predict(model, se.fit=TRUE)
predicted_glacier = predict(model,grd_hgt, se.fit=TRUE)

#diagnostic plot
#plot(hgt$x, hgt$y, col=cols[cut(unw$record_0, breaks = seq(-5,5,l=256))])
plot(ifgram$x, ifgram$y, col=cols[cut(ifgram$phase- predicted_hgt$fit, breaks = seq(-10,10,l=256))], pch=12, cex=(predicted_hgt$se.fit)*20)
points(grd_hgt$x, grd_hgt$y, col=cols[cut(predicted_glacier$fit, breaks = seq(-2,2,l=256))], pch=12, cex=2)


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
k <- krige.control(cov.model=unw.var.fit$cov.model, cov.pars=unw.var.fit$cov.pars, type.krige = 'OK', beta=0.2, dist.epsilon = 20)


#Locations of prediction
glacier_pg <- c(grd_hgt$x, grd_hgt$y)

#Krige 
aps.kriged <- krige.conv(res.geo,loc=glacier_pg,krige=k)


#Plot variogram and fit
plot(unw.var)
lines(res.var, col='red')
lines.variomodel(unw.var.fit)




plot(aps.geo$coords[,1], aps.geo$coords[,2], col=cols[cut(aps.geo$data, breaks = seq(-6,6,l=256))])

points(glacier_pg, col=cols[cut(aps.kriged$predict, breaks = seq(-6,6,l=256))], pch=0)

#Plot trend with height
plot(hgt$record_0, aps.geo$data)
