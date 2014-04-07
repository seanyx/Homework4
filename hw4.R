## This code extract pick time for P wave from "03092106124p" 
## and seismic station information from "wash2.sta.now", 
## then solve for the hypocenter location using 
## velocity model stored in file "PNW3.vel"
library(RSEIS) 
## used for extract arrival time, velocity model, 
## and station information from files 
library(GEOmap) ## used for convert lat/lon to UTM coordinates
library(akima) ## used for 2d interpolation contour plot for arrival time
library(Rquake) ## calculate the derivative matrix and arrival time


vv = getpfile("03092106124p" ) ## get arrival time

###  the station and arrival information is stored
###  in  a list of station stuff:  STAS

# data.frame(vv$STAS)

# ####  the earthquake location is stored in a list LOC
# vv$LOC

####   the arrival times are seconds after the event minute

# vv$LOC$yr
# vv$LOC$jd
# vv$LOC$hr
# vv$LOC$mi

sta = setstas("wash2.sta.now") ## get station information

vel = Get1Dvel("PNW3.vel", PLOT = FALSE) ## get velocity model

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
colnames(STA)=c('lat','lon','staz','sec','err')

## then convert the lat/lon to cartisian coordinates
## set the central latitude for coordinates conversion
orglat=median(STA[,'lat']) 
## set the central longitude for coordinates conversion
orglon=median(STA[,'lon']) 
## set the target projection to utm.sphr
proj=setPROJ(type=2,orglat,orglon) 
xy=GLOB.XY(STA[,1],STA[,2],proj)
STAxy=cbind(xy$x,xy$y,STA[,'sec'],STA[,'staz'],STA[,'err'])
rownames(STAxy)=rownames(STA)
colnames(STAxy)=c('x','y','sec','staz','err')

## plot the station locations colored by the arrival time
x=STAxy[,1]
y=STAxy[,2]
tt=STAxy[,3]
err=STAxy[,5]
staz=STAxy[,4]
sta.names=rownames(STAxy)
cols=(tt-min(tt))/(max(tt)-min(tt)) ## scale for color map
dev.new()
plot(x,y,col=gray(cols),pch=20)
points(x,y)
title(main='Station locations colored by arrival time',
      sub='The darker the smaller the arrival time')

## interpolate the timing and contour plot
## contour the time arrival, guess the initial lcoation
cex.staz=(staz-min(staz))/(max(staz)-min(staz))*3+1
xym=interp(x,y,tt)
ind=which(xym[[3]]==min(xym[[3]],na.rm=T),arr.ind=T)
image(xym[[1]],xym[[2]],xym[[3]],col=gray(seq(0.3,0.9,length=100)),ann=F)
points(x,y,pch=21,bg=gray(1,alpha=0.5),cex=cex.staz,col='black')
points(xym[[1]][ind[1]],xym[[2]][ind[2]],pch=10,col='red')
title(xlab='X (E-W) /km',ylab='Y (N-S) /km',
      main='Station locations and arrival time',
      sub='The darker the region the smaller the arrival time')

## set up initial guess
EQ=list(x=xym[[1]][ind[1]],
	y=xym[[2]][ind[2]],
	z=10,
	t=xym[[3]][ind[1],ind[2]])
#EQ=list(x=50,y=-50,z=10,t=xym[[3]][ind[1],ind[2]])
EQini=EQ
## assign tolerances
xtol=0.001
ytol=0.001
ztol=0.01
lambdareg=100 ## regularization factor
distwt=10 ## distance weighting factor

for (i in 1:10000) {
	delx=EQ$x-x
	dely=EQ$y-y
	## distance from previous EQ location guess to each station
	deltadis=sqrt(delx^2+dely^2) 
	## calculate the travel time and derivatives
	temp=GETpsTT(rep('P',length(tt)),eqz=EQ$z,staz=-staz,
		     delx=delx,dely=dely,deltadis=deltadis,vel=vel) 
	## G matrix to Ax=b (b is the residules, 
	## x is perturbation that needs to be solved)
	G=cbind(rep(1,nrow(temp$Derivs)),temp$Derivs) 

	observed=tt ## observed arrival time
	## create weighting matrix according to the distance
	wts = DistWeightXY(xy$x, xy$y, EQ$x, EQ$y, err, distwt)
	predictedTT=EQ$t+temp$TT
	## cors is the station corrections
	weights=wts

	## solve the linear equation
	S=svd(G)
	RHS=weights*(observed-predictedTT)
	LAM=diag(S$d/(S$d^2+lambdareg^2))
	per=S$v %*% LAM %*% t(S$u) %*% RHS
	
	## test the tolerance
	if (abs(per[2])<xtol 
	    & abs(per[3])<ytol 
	    & abs(per[4])<ztol) break
	
	## update the earthquake location using the perturbaton
	EQ$x=EQ$x+per[2]
	EQ$y=EQ$y+per[3]
	EQ$z=EQ$z+per[4]
	EQ$t=EQ$t+per[1]

}

print(paste('tolerance reached at step',i))  ## test if the result converge
points(EQ$x,EQ$y,pch=4,col='red')  ## plot the calculated earthquake location

## convert the earthquake location stored in the pick file to UTM projection
refxy=GLOB.XY(vv$LOC$lat,vv$LOC$lon,proj) 
refxy['z']=vv$LOC$z
refxy['sec']=vv$LOC$sec

points(refxy$x,refxy$y,pch=2,col='green')
legend('topright',legend=
       c('Station location: size by elevation','Initial guess',
	 'Calculated location','Location from pickfile'),
       pch=c(21,10,4,2),col=c('black','red','red','green'))

## Convert the EQ location back to LL format
EQLL=XY.GLOB(EQ$x,EQ$y,proj)
EQ$lat=EQLL$lat; EQ$lon=EQLL$lon-360
print(paste('The EQ is located at latitude',format(EQ$lat,digits=4),
	    ', longitude',format(EQ$lon,digits=4),', and depth',
	    format(EQ$z,digits=4),'km'))

## Here to calculate the timing err for final EQ locating
## EQ locating is an iterative process. Here I only used the last 
## iteration to calculate the locating error, thus the results only
## represents the stability of the last iteration, not the whole process.

delx=EQ$x-x
dely=EQ$y-y
## distance from previous EQ location guess to each station
deltadis=sqrt(delx^2+dely^2) 
## calculate the travel time and derivatives
temp=GETpsTT(rep('P',length(tt)),eqz=EQ$z,staz=-staz,
	     delx=delx,dely=dely,deltadis=deltadis,vel=vel) 
## G matrix to Ax=b (b is the residules, 
## x is perturbation that needs to be solved)
G=cbind(rep(1,nrow(temp$Derivs)),temp$Derivs) 

observed=tt ## observed arrival time
## create weighting matrix according to the distance
wts = DistWeightXY(xy$x, xy$y, EQ$x, EQ$y, err, distwt)
predictedTT=EQ$t+temp$TT
## cors is the station corrections
weights=wts
## plot the weighting for each stations
dev.new(); plot(STAxy[,'x'],STAxy[,'y'],
		cex=(wts-min(wts))/(max(wts)-min(wts))*3+1,ann=F)
title(main='Example of station weighting',
      xlab='E-W',
      ylab='N-S')

## calculate the residual
res=observed-predictedTT
## calculate the (weighted) root mean square
rms=sqrt(mean(res^2))
wrms=sqrt(mean((weights*res)^2))

## Here starts the locating error estimation (in km)
dels=rep(NA,4)
S=svd(G)
LAM=diag(S$d/(S$d^2+lambdareg^2))
Gdagger=S$v %*% LAM %*% t(S$u)
covD=diag(err)
covB=Gdagger %*% covD %*% t(Gdagger)

uncer=1.96*sqrt(diag(covB)) ## 95% interval
print(paste('error in E-W direction:',format(uncer[1],digits=4),'km',
	    '\nerror in N-S direction:',format(uncer[2],digits=4),'km',
	    '\nerror in z direction:',format(uncer[3],,digits=4),'km'))
