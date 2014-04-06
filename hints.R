library(RSEIS)
library(GEOmap)
library(akima)
vv = getpfile("03092106124p" )

####  the station and arrival information is stored
###  in  a list of station stuff:  STAS

# data.frame(vv$STAS)


# ####  the earthquake location is stored in a list LOC
# vv$LOC


####   the arrival times are seconds after the event minute

# vv$LOC$yr
# vv$LOC$jd
# vv$LOC$hr
# vv$LOC$mi

#####

##### The station locations are stored in two files
##  I am not sure which is best -
##   maybe either
##wash2.sta.now
##wash2.sta


sta = setstas("wash2.sta.now")

vel <- Get1Dvel("PNW3.vel", PLOT = FALSE)

# names(sta)

####   from this information you should be able to locate the earthquakes

## find the stations that both in the pick file and the station file
indp=which(vv$STAS$phase=='P')  ## only use P wave here
name.stations=intersect(vv$STAS$name[indp],sta$name)
ind1=match(name.stations,vv$STAS$name)
tt=vv$STAS$sec[ind1]
err=vv$STAS$err[ind1]
ind2=match(name.stations,sta$name)
STA=cbind(sta$lat[ind2],sta$lon[ind2],sta$z[ind2],tt,err)
rownames(STA)=name.stations
colnames(STA)=c('lat','lon','z','sec','err')

## then convert the lat/lon to cartisian coordinates
orglat=median(STA[,'lat'])
orglon=median(STA[,'lon'])
proj=setPROJ(type=2,orglat,orglon)
xy=GLOB.XY(STA[,1],STA[,2],proj)
STAxy=cbind(xy$x,xy$y,STA[,'sec'],STA[,'err'])
rownames(STAxy)=rownames(STA)
colnames(STAxy)=c('x','y','sec','err')

## plot the station locations colored by the arrival time
x=STAxy[,1]
y=STAxy[,2]
tt=STAxy[,3]
cols=(tt-min(tt))/(max(tt)-min(tt))
dev.new()
plot(x,y,col=gray(cols),pch=20)
points(x,y)
title(main='Station locations colored by arrival time',sub='The darker the smaller the arrival time')

## interpolate the timing and contour plot
## contour the time arrival, guess the initial lcoation
xym=interp(x,y,tt)
ind=which(xym[[3]]==min(xym[[3]],na.rm=T),arr.ind=T)
image(xym[[1]],xym[[2]],xym[[3]],col=gray(seq(0.3,0.9,length=100)),ann=F)
points(x,y,pch=20,col=gray(1,alpha=0.5))
points(x,y)
points(xym[[1]][ind[1]],xym[[2]][ind[2]],pch=10,col='red')
title(xlab='X /km',ylab='Y /km',main='Station locations and arrival time',sub='The darker the smaller the arrival time')

## set up initial guess
EQ=list(x=xym[[1]][ind[1]],y=xym[[2]][ind[2]],z=10,t=xym[[3]][ind[1],ind[2]])
#EQ=list(x=50,y=-50,z=10,t=xym[[3]][ind[1],ind[2]])
EQini=EQ
## assign tolerances
xtol=0.001
ytol=0.001
ztol=0.01
lambdareg=100

for (i in 1:10000) {
	delx=EQ$x-x
	dely=EQ$y-y
	deltadis=sqrt(delx^2+dely^2) ## distance from initial guess to each station
	temp=GETpsTT(rep('P',length(tt)),eqz=EQ$z,staz=0,delx=delx,dely=dely,deltadis=deltadis,vel=vel) ## calculate the travel time and derivatives
	G=cbind(rep(1,nrow(temp$Derivs)),temp$Derivs) ## G matrix to Ax=b (b is the residules, x is perturbation that is solved for)

	observed=tt ## observed arrival time
	## create weighting matrix according to the distance
	distwt=10
	wts = DistWeightXY(xy$x, xy$y, EQ$x, EQ$y, STAxy[,'err'], distwt)
	predictedTT=EQ$t+temp$TT
	## cors
	weights=wts

	## solve the linear equation
	S=svd(G)
	RHS=weights*(observed-predictedTT)
	LAM=diag(S$d/(S$d^2+lambdareg^2))
	per=S$v %*% LAM %*% t(S$u) %*% RHS
	
	## test the tolerance
	if (abs(per[2])<xtol & abs(per[3])<ytol & abs(per[4])<ztol) break
	
	## update the earthquake location using the perturbaton
	EQ$x=EQ$x+per[2]
	EQ$y=EQ$y+per[3]
	EQ$z=EQ$z+per[4]
	EQ$t=EQ$t+per[1]

}

print(paste('tolerance reached at step',i))  ## test if the result converge
points(EQ$x,EQ$y,pch=20,col='blue')  ## 
