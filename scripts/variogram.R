library(geoR)
library(akima)
#Read data
unw<-read.table('/home/baffelli/PhD/trunk/Bisgletscher_pol/ipt/20150803_060519_stack_20_AAAl.aps.csv', sep=',', head=T)
#Read heights 
hgt<-read.table('/home/baffelli/PhD/trunk/Bisgletscher_pol/ipt/20150803_060519_stack_20_AAAl.hgt.csv', sep=',', head=T)


#Color ramp for plotting
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
#Model as a function of elevation
model_eq <- unw$record_0 ~ hgt$x + hgt$y + I(hgt$x^2) + I(hgt$y^2) + hgt$x*hgt$y + hgt$record_0 + I(hgt$record_0^2) + I(hgt$y^2*hgt$x^2)
model <- lm(model_eq)
predicted_hgt = predict(model)

#diagnostic plot
plot(hgt$x, hgt$y, col=cols[cut(unw$record_0, breaks = seq(-5,5,l=256))])
points(hgt$x, hgt$y, col=cols[cut(predicted_hgt, breaks = seq(-5,5,l=256))], pch=12, cex=(unw$record_0 - predicted_hgt)**2/5)



#Create dataframe
unw.geo <- as.geodata(unw, coords.col = 1:2, data.col=16)
res.geo <- as.geodata(unw - predicted_hgt, coords.col = 1:2, data.col=16)
#Compute variogram
unw.var <- variog(unw.geo, breaks=seq(30, 3000, l=50))
res.var <- variog(res.geo, breaks=seq(30, 3000, l=50))



#Fit model
unw.var.fit <- variofit(unw.var, ini.cov.pars = c(0.5,1000),cov.model ="exponential",fix.nugget=FALSE, nugget=0.5)
#Fit model
res.var.fit <- variofit(res.var, ini.cov.pars = c(0.5,1000),cov.model ="exponential",fix.nugget=FALSE, nugget=0.5)

#Plot models
plot(unw.var,ylim=c(0,10))
lines(res.var, col='red')
lines.variomodel(unw.var.fit)
lines.variomodel(res.var.fit,col='red')


#Create trend
aps_trend <-trend.spatial(model_eq,aps.geo)
#Kriging object
k <- krige.control(cov.model=unw.var.fit$cov.model, cov.pars=unw.var.fit$cov.pars, type.krige = 'OK', beta=0.2, dist.epsilon = 20)


#Locations of prediction
glacier_pg <- read.table('/home/baffelli/PhD/trunk/Bisgletscher_pol/ipt/20150803_060519_stack_20_AAAl.glacier.csv', sep=',', head=T)

#Krige 
aps.kriged <- krige.conv(aps.geo,loc=glacier_pg,krige=k)


#Plot variogram and fit
plot(unw.var)
lines(res.var, col='red')
lines.variomodel(unw.var.fit)




plot(aps.geo$coords[,1], aps.geo$coords[,2], col=cols[cut(aps.geo$data, breaks = seq(-6,6,l=256))])

points(glacier_pg, col=cols[cut(aps.kriged$predict, breaks = seq(-6,6,l=256))], pch=0)

#Plot trend with height
plot(hgt$record_0, aps.geo$data)
