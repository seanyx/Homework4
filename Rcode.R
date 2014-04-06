setwd('/Users/yangxiao/Dropbox/class2014spring/InverseTheory/hw4')

library(RSEIS)
library(Rquake)

## arrival time
vv=getpfile('03092106124p')
data.frame(vv$STAS)

vv$LOC

## all stations locations
sta=setstas('wash2.sta.now')