rm(list=ls())
#################################################################
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

#library(latex2exp)

require(ggplot2)
library(cowplot)


meanLOS <- 6
del <-   (1/300) * meanLOS
omega <- 0.07 * meanLOS
lamb <- 0.285

mu <- 0.05 * meanLOS
sig <- 0.388

fC <- 0.075


fA <- 0.24
omega0 <- 0.181

Y <- 0.16                  #prevalence target
X <- 14.5 / 10000 * meanLOS  #incidence target 
targetR <- 3.6             #risk ratio target

losMGF <- function(x) 1/(1-x)

K <- function(r) -(1-losMGF(r))/r
K <- Vectorize(K)

states <- c('SN', 'CN', 'IN', 'SA', 'CA', 'IA', 'SV', 'CV', 'IV')

Tabx <- function(bet,p,ma,mp,omega,mu,lamb){
	rbind(c(-bet-omega    , del           , 0             , (1-sig)*mu, 0         , 0        , 0      , 0        , 0),
		c(bet           , -omega-p-del  , 0             , 0         , (1-sig)*mu, 0        , 0      , 0        , 0),
		c(0             , p             , -omega        , 0         , 0         ,(1-sig)*mu, 0      , 0        , 0),
		c((1-lamb)*omega, 0             , 0             , -mu-bet   , del       , 0        , 0      , 0        , 0),
		c(0             , (1-lamb)*omega, 0             , bet       , -p-del-mu , 0        , 0      , 0        , 0),
		c(0             , 0             , (1-lamb)*omega, 0         , p         , -mu      , 0      , 0        , 0),
		c(lamb*omega    , 0             , 0             , sig*mu    , 0         , 0        , -ma*bet, del      , 0),
		c(0             , lamb*omega    , 0             , 0         , sig*mu    , 0        , ma*bet , -mp*p-del, 0),
		c(0             , 0             , lamb*omega    , 0         , 0         , sig*mu   , 0      , mp*p     , 0))
}

getInitAbx <- function(omega0,lamb){
	c((1-fC)*(1-fA-omega0),
	  fC*(1-fA-omega0),
	  0,
	  (1-fC)*(1-lamb)*(fA+omega0),
	  fC*(1-lamb)*(fA+omega0),
	  0,
	  (1-fC)*lamb*(fA+omega0),
	  fC*lamb*(fA+omega0),
	  0)
}

initAbx <- getInitAbx(omega0,lamb)

names(initAbx) <- states

getEq <- function(Tfull,init){
	T <- Tfull[-1,-1] - Tfull[-1,1]
	rhs <- -Tfull[-1,1]

	Tinf <- solve(T,rhs)

	eig <- eigen(T)
	r <- eig$values
	V <- eig$vectors
	
	tmp <- solve(V,as.vector((init[-1]-Tinf)))
	eq <- V %*% (tmp*K(r)) + Tinf
	c(1-sum(eq),eq)
}

getStatesAtT <- function(Tfull,init,tm){
	T <- Tfull[-1,-1] - Tfull[-1,1]
	rhs <- -Tfull[-1,1]

	Tinf <- solve(T,rhs)

	eig <- eigen(T)
	r <- eig$values
	V <- eig$vectors

	tmp <- solve(V, as.vector((init[-1]-Tinf)))

	x <- V %*% (tmp*exp(r*tm)) + Tinf
	c(1-sum(x),x)
}

getXY <- function(bet,p,ma,mp){
	eq <- getEq(Tabx(bet,p,ma,mp,omega,mu,lamb),initAbx)
	names(eq) <- states
	Y <- 1-sum(eq[c('SN','SA','SV')])
	X <- p*(eq['CN']+eq['CA']+mp*eq['CV'])
	c(X,Y)
}

#"Bar" states:
#Recent or current abx (A)
#Susceptible, no recent or current abx (SNbar)
#Colonized, no recent or current abx (CNbar)
#Infected, no recent or current abx (INbar)

barStates <- c('A','SNbar','CNbar','INbar')

TabxBar <- function(bet,p){
	rbind(c(0, omega         , omega         , omega),
		c(0, -bet-omega    , del           , 0),
		c(0, bet           , -omega-p-del  , 0),
		c(0, 0             , p             , -omega))
}

initBar <- c(fA+omega0, (1-fC)*(1-fA-omega0), fC*(1-fA-omega0), 0)

getBP <- function(ma,mp) optim(fn = function(x) sum((getXY(x[1],x[2],ma,mp) - c(X,Y))^2), c(beta=.05,p=.03), control=list(abstol=1e-15))$par

getR <- function(ma,mp){

	bp <- getBP(ma,mp)

	bet <- bp['beta']
	p <- bp['p']

	xt <- function(tm){
		out <- getStatesAtT(Tabx(bet,p,ma,mp,omega,mu,lamb), initAbx, tm)
		names(out) <- states
		out
	}

	xtbar <- function(tm){
		out <- getStatesAtT(TabxBar(bet,p), initBar, tm)
		names(out) <- barStates
		out
	}

	R <- function(tm){
		st <- xt(tm)
		stBar <- xtbar(tm)
	
		mp * st['CV']/(st['SV']+st['CV']) / (stBar['CNbar']/(stBar['SNbar']+stBar['CNbar']))
	}

	denomFn <- function(tm){
		st <- xt(tm)
		stBar <- xtbar(tm)

		(st['SV']+st['CV']+stBar['SNbar']+stBar['CNbar'])*exp(-tm)
	}

	numerFn <- function(tm){
		st <- xt(tm)
		stBar <- xtbar(tm)

		R <- mp * st['CV']/(st['SV']+st['CV']) / (stBar['CNbar']/(stBar['SNbar']+stBar['CNbar']))
		R*(st['SV']+st['CV']+stBar['SNbar']+stBar['CNbar'])*exp(-tm)
	}
	
	integrate(Vectorize(numerFn),0,Inf)$value / integrate(Vectorize(denomFn),0,Inf)$value
}

numPlot <- 20

mp_max <- optimize(function(x) (getR(1,x)-targetR)^2, c(1,5))$minimum

mps <- seq(1,mp_max,len=numPlot)
mas <- rep(0,numPlot)
bets <- rep(0,numPlot)
ps <- rep(0,numPlot)

for(i in 1:(numPlot-1)){
	mas[i] <- optimize(function(x) (getR(x,mps[i])-targetR)^2, c(1,25))$minimum
	bp <- getBP(mas[i],mps[i])
	bets[i] <- bp[1]
	ps[i] <- bp[2]
}
mas[numPlot] <- 1
bpEnd <- getBP(mas[numPlot],mps[numPlot])
bets[numPlot] <- bpEnd[1]
ps[numPlot] <- bpEnd[2]


 

#x11()
#plot(mps,mas,type='l',xlab = 'Progression effect ($m_p$)',ylab='Acquisition effect ($m_a$)',cex.axis=1.25,cex.lab=1.25,lwd = 2)
 
#points(mps[3],mas[3],type='p',col = 'black',pch=15,cex = 1.5)
 
#points(mps[16],mas[16],type='p',col = 'black',pch=16,cex = 1.5)

#n = 2
#cols = gg_color_hue(n)

#dev.new(width = 4, height = 4)
#plot(1:n, pch = 16, cex = 2, col = cols)

fig1dat<-data.frame(mps,mas)

# Fig1<- ggplot(fig1dat,aes(x=mps*mu, y=mas*mu))+
# 	geom_line(size=1.2) +
# 	theme_bw(base_size=15) +
# 	labs(x=bquote("Progression Effect"~m[p]), y=bquote("Acquisition Effect"~m[a])) +
# 	geom_label(label='1', x=mps[3]*mu, 	y=mas[3]*mu, 
# 	           label.padding = unit(0.25, "lines"), # Rectangle size around label
# 	           label.size = 0.35,	label.r = unit(0,'lines'),	color = "#C77CFF",
# 	           fill="#C77CFF") +
# 	geom_label(label='1',	x=mps[16]*mu,	y=mas[16]*mu,
# 	           label.padding = unit(0.25, "lines"), # Rectangle size around label
# 	           label.size = 0.35,	label.r = unit(.5, "lines"),	color = "#7CAE00",
# 	           fill="#7CAE00") + 
# 	geom_segment(aes(x=mps[3]*mu,y=0, xend=mps[3]*mu, yend=mas[3]*mu),linetype='dotted',
# 	             size=1, color = "#C77CFF") +
# 	geom_segment(aes(x=0,y=mas[3]*mu, xend=mps[3]*mu, yend=mas[3]*mu),linetype='dotted',
# 	             size=1, color = "#C77CFF") + 
# 	geom_segment(aes(x=mps[16]*mu,y=0, xend=mps[16]*mu, yend=mas[16]*mu),linetype='dotted',
# 	             size=1, color = "#7CAE00") +
# 	geom_segment(aes(x=0,y=mas[16]*mu, xend=mps[16]*mu, yend=mas[16]*mu),linetype='dotted',
# 	             size=1, color = "#7CAE00") 

Fig1<- ggplot(fig1dat,aes(x=mps, y=mas))+
  geom_line(size=1.2) +
  theme_bw(base_size=15) +
  labs(x=bquote("Progression Effect"~m[p]), y=bquote("Acquisition Effect"~m[a])) +
  geom_label(label='1', x=mps[3], 	y=mas[3], 
             label.padding = unit(0.25, "lines"), # Rectangle size around label
             label.size = 0.35,	label.r = unit(0,'lines'),	color = "#C77CFF",
             fill="#C77CFF") +
  geom_label(label='1',	x=mps[16],	y=mas[16],
             label.padding = unit(0.25, "lines"), # Rectangle size around label
             label.size = 0.35,	label.r = unit(.5, "lines"),	color = "#7CAE00",
             fill="#7CAE00") + 
  geom_segment(aes(x=mps[3],y=0, xend=mps[3], yend=mas[3]),linetype='dotted',
               size=1, color = "#C77CFF") +
  geom_segment(aes(x=0,y=mas[3], xend=mps[3], yend=mas[3]),linetype='dotted',
               size=1, color = "#C77CFF") + 
  geom_segment(aes(x=mps[16],y=0, xend=mps[16], yend=mas[16]),linetype='dotted',
               size=1, color = "#7CAE00") +
  geom_segment(aes(x=0,y=mas[16], xend=mps[16], yend=mas[16]),linetype='dotted',
               size=1, color = "#7CAE00") 


Fig1

############################################################################
getBetaHat <- function(bet,p,ma,mp,omega,mu,lamb,initAbx,ms,mi){
	eq <- getEq(Tabx(bet,p,ma,mp,omega,mu,lamb),initAbx)
	names(eq) <- states
	betaHat <- as.numeric(bet/(eq['CN'] + eq['CA'] + ms*eq['CV'] + mi*(eq['IN'] + eq['IA'] + ms*eq['IV'])))
	betaHat
}

getInterventionResults <- function(bet,p,ma,mp,ms,mi,startReduct){
	betaHatFix <- getBetaHat(bet,p,ma,mp,omega,mu,lamb,initAbx,ms,mi)

	omegaNew <- omega*(1-startReduct)
	initAbxNew <- getInitAbx(omega0*(1-startReduct),lamb)

	betNew <- optimize(function(bet) (getBetaHat(bet,p,ma,mp,omegaNew,mu,lamb,initAbxNew,ms,mi) - betaHatFix)^2,c(0,10),tol=1e-15)$minimum

	TabxNew <- Tabx(betNew,p,ma,mp,omegaNew,mu,lamb)
	eqNew <- getEq(TabxNew,initAbxNew)
	names(eqNew) <- states
	
	acqNew <- betNew*(eqNew['SN'] + eqNew['SA'] + ma*eqNew['SV'])
	infNew <- p*(eqNew['CN'] + eqNew['CA'] + mp*eqNew['CV'])
	
	c(acq=acqNew,inf=infNew)
}

getInterventionEffect <- function(ind,ms,mi,startReduct){
	outBase <- getInterventionResults(bets[ind],ps[ind],mas[ind],mps[ind],ms,mi,0)
	
	startEffect <- matrix(0,length(startReduct),2)
	if(!identical(startReduct,0)){		
		for(i in 1:length(startReduct)){
			outNew <- getInterventionResults(bets[ind],ps[ind],mas[ind],mps[ind],ms,mi,startReduct[i])
 		#startEffect[i,] <- (outNew-outBase)/outBase
		startEffect[i,] <- outNew
		}
	}

	startEffect
}

sr <- seq(-0.5,0.5,len=20)

ie3 <- getInterventionEffect(3,1,10,sr)
ie16 <- getInterventionEffect(16,1,10,sr)


## figure 2

fig2dat <- data.frame(omega*(1-sr), ie3, ie16)
colnames(fig2dat)<-c('sr','ie3acq','ie3inf','ie16acq','ie16inf')

Fig2a <- ggplot(fig2dat) + 
	theme_bw(base_size=15)+
	geom_line(aes(sr,ie3acq),size=2, color= "#C77CFF") +
	geom_line(aes(sr,ie16acq),linetype='dashed',size=2,color = "#7CAE00") +
	xlab(bquote("Daily antibiotic prescribing rate,"~omega~"")) +
	ylab('Facilty onset acquisition rate')+
	theme(
		plot.title = element_text(color="black", size=14,
			 face="bold"))
Fig2a

Fig2b <- ggplot(fig2dat) + 
	theme_bw(base_size=14)+
	geom_line(aes(sr,ie3inf),size=2, color= "#C77CFF") +
	geom_line(aes(sr,ie16inf),linetype='dashed',size=2,color = "#7CAE00") +
	xlab(bquote("Daily antibiotic prescribing rate,"~omega~"")) +
	ylab('Facility onset infection rate')+
	theme(
		plot.title = element_text(color="black", size=14,
			 face="bold"))
Fig2b

plot_grid(Fig2a, Fig2b, labels = "AUTO")

#x11()
#plot(1-sr,ie3[,2],type='l',xlab=expression(paste('Decrease in antibiotic prescribing rate (',phi,')')),ylab='Change in facility-onset infections',cex.axis=1.5,cex.lab=1.5)
#lines(1-sr,ie16[,2],lty = 2)


#x11()
#plot(1-sr,ie3[,1],type='l',xlab=expression(paste('Decrease in antibiotic prescribing rate (',phi,')')),ylab='Change in facility acquisitions',cex.axis=1.5,cex.lab=1.5)
#lines(1-sr,ie16[,1],lty = 2)


############################################################################################################

getInterventionResultsMu <- function(bet,p,Ma,Mp,ms,mi,startReduct){
  ma <- Ma*mu
  mp <- Mp*mu
  betaHatFix <- getBetaHat(bet,p,Ma,Mp,omega,mu,lamb,initAbx,ms,mi)

  muNew <- mu*(1+startReduct)
  MaNew <- ma/muNew
  MpNew <- mp/muNew

  betNew <- optimize(function(bet) (getBetaHat(bet,p,MaNew,MpNew,omega,muNew,lamb,initAbx,ms,mi) - betaHatFix)^2,c(0,10),tol=1e-15)$minimum

  TabxNew <- Tabx(betNew,p,MaNew,MpNew,omega,muNew,lamb)
  eqNew <- getEq(TabxNew,initAbx)
  names(eqNew) <- states

  acqNew <- betNew*(eqNew['SN'] + eqNew['SA'] + MaNew*eqNew['SV'])
  infNew <- p*(eqNew['CN'] + eqNew['CA'] + MpNew*eqNew['CV'])

  c(acq=acqNew,inf=infNew)
}



getInterventionEffect <- function(ind,ms,mi,startReduct){
  outBase <- getInterventionResultsMu(bets[ind],ps[ind],mas[ind],mps[ind],ms,mi,0)

  startEffect <- matrix(0,length(startReduct),2)
  if(!identical(startReduct,0)){
    for(i in 1:length(startReduct)){
      outNew <- getInterventionResultsMu(bets[ind],ps[ind],mas[ind],mps[ind],ms,mi,startReduct[i])
    #startEffect[i,] <- (outNew-outBase)/outBase
	  startEffect[i,] <- outNew
    }
  }

  startEffect
}

sr <- seq(-0.2,0.5,len=20)
ie3 <- getInterventionEffect(3,1,10,sr)
ie16 <- getInterventionEffect(16,1,10,sr)

fig3dat <- data.frame(1/(mu*(1+sr)), ie3, ie16)
colnames(fig3dat)<-c('sr','ie3acq','ie3inf','ie16acq','ie16inf')

Fig3a <- ggplot(fig3dat) + 
	theme_bw(base_size=14)+
	geom_line(aes(sr,ie3acq),size=2, color= "#C77CFF") +
	geom_line(aes(sr,ie16acq),linetype='dashed',size=2,color = "#7CAE00") +
	xlab(bquote("Duration of antibiotic use,"~mu^-1~"")) +
	ylab('Facility onset acquisition rate')+
	theme(
		plot.title = element_text(color="black", size=14,
			 face="bold"))
Fig3a

Fig3b <- ggplot(fig3dat) + 
	theme_bw(base_size=14)+
	geom_line(aes(sr,ie3inf),size=2, color= "#C77CFF") +
	geom_line(aes(sr,ie16inf),linetype='dashed',size=2,color = "#7CAE00") +
	xlab(bquote("Duration of antibiotic use,"~mu^-1~"")) +
	ylab('Facility Onset Infection rate')+
	theme(
		plot.title = element_text(color="black", size=14,
			 face="bold"))
Fig3b


plot_grid(Fig3a, Fig3b, labels = "AUTO")

#x11()
#plot(1+sr,ie3[,2],type='l',xlab=expression(paste('Increase in antibiotic stopping rate (',phi,')')),ylab='Change in facility-onset infections',cex.axis=1.5,cex.lab=1.5)
#lines(1+sr,ie16[,2],lty = 2)

#x11()
#plot(1+sr,ie3[,1],type='l',xlab=expression(paste('Increase in antibiotic stopping rate (',phi,')')),ylab='Change in facility acquisitions',cex.axis=1.5,cex.lab=1.5)
#lines(1+sr,ie16[,1],lty = 2)


#######################################################################################

getInterventionResultsLamb <- function(bet,p,ma,mp,ms,mi,lamb,startReduct){
  betaHatFix <- getBetaHat(bet,p,ma,mp,omega,mu,lamb,initAbx,ms,mi)
  
  lambNew <- lamb*(1+startReduct)
  initAbxNew <- getInitAbx(omega0,lambNew)
  
  betNew <- optimize(function(bet) (getBetaHat(bet,p,ma,mp,omega,mu,lambNew,initAbxNew,ms,mi) - betaHatFix)^2, c(0,10),tol=1e-15)$minimum
  
  TabxNew <- Tabx(betNew,p,ma,mp,omega,mu,lambNew)
  eqNew <- getEq(TabxNew,initAbxNew)
  names(eqNew) <- states
  
  acqNew <- betNew*(eqNew['SN'] + eqNew['SA'] + ma*eqNew['SV'])
  infNew <- p*(eqNew['CN'] + eqNew['CA'] + mp*eqNew['CV'])
  
  c(acq=acqNew,inf=infNew)
}

getInterventionEffect <- function(ind,ms,mi,startReduct){
  outBase <- getInterventionResultsLamb(bets[ind],ps[ind],mas[ind],mps[ind],ms,mi,lamb,0)
  
  startEffect <- matrix(0,length(startReduct),2)
  if(!identical(startReduct,0)){		
    for(i in 1:length(startReduct)){
      outNew <- getInterventionResultsLamb(bets[ind],ps[ind],mas[ind],mps[ind],ms,mi,lamb,startReduct[i])
    #startEffect[i,] <- (outNew-outBase)/outBase
	  startEffect[i,] <- outNew
    }
  }
  
  startEffect
}

sr <- seq(-0.5,0.5,len=20)  
ie3 <- getInterventionEffect(3,1,10,sr)
ie16 <- getInterventionEffect(16,1,10,sr)

fig4dat <- data.frame(sr, ie3, ie16)
colnames(fig4dat)<-c('sr','ie3acq','ie3inf','ie16acq','ie16inf')

Fig4a <- ggplot(fig4dat) + 
	theme_bw(base_size=12)+
	geom_line(aes(lamb*(1+sr),ie3acq),size=2, color= "#C77CFF") +
	geom_line(aes(lamb*(1+sr),ie16acq),linetype='dashed',size=2,color = "#7CAE00") +
  xlab(bquote("Daily rate of prescribing high-risk antibiotics")) +
  ylab('Facility onset acquisition rate')+
	theme(
	 plot.title = element_text(color="black", size=14,
			 face="bold"))
Fig4a

Fig4b <- ggplot(fig4dat) + 
	theme_bw(base_size=12)+
	geom_line(aes(lamb*(1+sr),ie3inf),size=2, color= "#C77CFF") +
	geom_line(aes(lamb*(1+sr),ie16inf),linetype='dashed',size=2,color = "#7CAE00") +
	xlab(bquote("Daily rate of prescribing high-risk antibiotics")) +
	ylab('Facility onset infection rate')+
	theme(
	 plot.title = element_text(color="black", size=14,
			 face="bold"))
Fig4b
plot_grid(Fig4a, Fig4b, labels = "AUTO")


#x11()
#plot(1+sr,ie3[,2],type='l',xlab=expression(paste('Increase in high-risk antibiotic rate (',phi,')')),ylab='Change in facility-onset infections',cex.axis=1.5,cex.lab=1.5)
#lines(1+sr,ie16[,2],lty = 2)


#x11()
#plot(1+sr,ie3[,1],type='l',xlab=expression(paste('Increase in high-risk antibiotic rate (',phi,')')),ylab='Change in facility acquisitions',cex.axis=1.5,cex.lab=1.5)
#lines(1+sr,ie16[,1],lty = 2)
# 
# 
# ############################################################################################################
# 



sr <- 0.25
print("higher acq: acq vs inf")
print(getInterventionEffect(3,1,1,sr)) #Case 1 Reduced
print(getInterventionEffect(3,1,10,sr))
print(getInterventionEffect(3,3,1,sr))
print(getInterventionEffect(3,3,10,sr))
print("higher prog: acq vs inf")
print(getInterventionEffect(16,1,1,sr))
print(getInterventionEffect(16,1,10,sr))
print(getInterventionEffect(16,3,1,sr))
print(getInterventionEffect(16,3,10,sr))
