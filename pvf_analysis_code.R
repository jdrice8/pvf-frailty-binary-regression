library(BHH2)
library(lme4)

# Laplace transform for a gamma random variable
laplace.gamma <- function(x,sigma2) (1+sigma2*x)^(-1/sigma2)

# Laplace transform for the PVF random variable
laplace.pvf <- function(x,sigma2,xi) {
	if(xi!=0) {
		rho <- (xi+1)/(xi*sigma2)
		#if(xi>0) print(exp(-rho)) else print(rho)
		nu <- rho*xi
	
		if(sign(rho)==sign(xi)) L <- exp(-rho*(1-(nu/(nu+x))^xi)) else L <- rep(1,length(x)) #stop('rho and xi must have the same sign')
	} else L <- laplace.gamma(x,sigma2=sigma2)
	L
}

# number of intervals
K <- 6

# code to generate matrices needed for Conaway's method for calculating the exact marginal likelihood
d <- iA <- Amatch <- list()
for(k in 1:K) {
	d[[k]] <- ifelse(ffDesMatrix(k)==1,1,0)
	dk <- data.frame(d[[k]])
	if(k==1) colnames(dk) <- 'X1'
	Ak <- model.matrix(as.formula(paste('~',paste(paste('X',1:k,sep=''),collapse='*'),sep='')),data=dk)[rev(1:2^k),]
	iA[[k]] <- solve(Ak)
	Amatch[[k]] <- unlist(lapply(strsplit(gsub('[^0-9]','',colnames(Ak)),''),function(x) paste(x,collapse=''))) 
	
}

# code to get the exact marginal likelihood using Conaway's method (one subject)
get.exact0 <- function(y,eta,sigma2,xi=NULL) {
	if(is.null(xi)) xi <- 0
	k <- length(y)

	Ainv <- iA[[k]] 
	pos.index <- Amatch[[k]] 

	(Ainv %*% laplace.pvf(d[[k]] %*% eta,sigma2=sigma2,xi=xi))[which(pos.index== paste(which(y==1),collapse=''))] 	
}

# function to get the marginal likelihood under the PVF frailty model (all subjects)
obs.logL <- function(theta,z,y,id) {
	x <- z
	if(dim(z)[2]==length(theta)-2) {
		s <- exp(tail(theta,2))[1]
		xi <- exp(tail(theta,2))[2]-1
		b <- head(theta,-2)
	} else {
		s <- exp(tail(theta,1))
		xi <- 0
		b <- head(theta,-1)
	}

	sum(log(as.numeric(by(data.frame(y,x),id,function(u) get.exact0(y=u[!is.na(u[,1]),1],eta=exp(as.matrix(u[!is.na(u[,1]),2:dim(u)[2]]) %*% b),sigma2=s,xi=xi)))))

}

# conditional likelihood (given the random effect) for discrete mixture model
condtl.bern <- function(v,y,eta,tau,phi) {
	(plogis(-phi)*all(y==0) + plogis(phi)*prod(dbinom(y,1,1-exp(-eta*v))))*dlnorm(v,sdlog=tau)
}

# numerical integration to find marginal likelihood under discrete mixture model
exact.logL <- function(theta,y,x,id) {
	tau <- exp(tail(theta,2)[1])
	phi <- tail(theta,1)
	prodL <- function(u) {
		out <- try(integrate(function(v) sapply(v,condtl.bern,y=u[,1],eta=exp(as.matrix(u[,2:dim(u)[2]]) %*% head(theta,-2)),tau=tau,phi=phi),0,Inf),silent=TRUE)
		if(class(out)!='try-error') L <- out$val else L <- 1 #; print(out)}
		L

	}

	L <- as.numeric(by(data.frame(y,x),id, prodL))

	-sum(log(L))
}

# read in the data
example.data <- read.csv(file='example_data.csv')

# fit the model assuming no frailty
naive.fit <- glm(y ~ treat+as.factor(timepoint)+race+edu+cent.age+treat:as.factor(timepoint), data=example.data,
	family=binomial(link=cloglog), x=TRUE,y=TRUE)

# fit the Gaussian intercept model
gaussint.fit <- glmer(y ~ treat+as.factor(timepoint)+race+edu+cent.age+treat:as.factor(timepoint)+(1|id), data=example.data,
	family=binomial(link=cloglog))

# fit the Gaussian intercept with discrete mixture model
discmix.fit <- optim(c(summary(gaussint.fit)$coef[,1],0,0),exact.logL,hessian=TRUE,
	x=naive.fit$x,y=naive.fit$y,id=example.data[,1],
	method='BFGS',control=list(trace=1,REPORT=1,fnscale=1))

# fit the PVF frailty model
pvf.fit <- optim(c(summary(gaussint.fit)$coef[,1],0,0),obs.logL,hessian=TRUE,
	z=naive.fit$x,y=naive.fit$y,id=example.data[,1],
	method='BFGS',control=list(trace=1,REPORT=1,fnscale=-1))

# examine the results: left column is estimate, right column is standard error
summary(naive.fit)$coef[,1:2]
summary(gaussint.fit)$coef[,1:2]
cbind(pvf.fit$par,sqrt(diag(solve(-pvf.fit$hess))))
cbind(discmix.fit$par,sqrt(diag(solve(discmix.fit$hess))))
