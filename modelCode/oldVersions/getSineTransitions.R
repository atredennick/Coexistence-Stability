x = c(1:101)
y = 0.15*sin(0.1*pi*(x-1))+0.5
y2 = 0.15*cos(0.1*pi*(x+4))+0.5

# plot(x, y, type="l", ylim=c(0,1))
# lines(x, y2, lty=2)

gVecTmp <- cbind(y, y2)
rm(y); rm(y2); rm(x)