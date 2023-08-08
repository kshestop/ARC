
#########################################################################################
### Functions to fit the individual curves and obtain unobserved species estimates

#Estimate the tail value given a function for a large number of reads
#Used if the integration returns an error
est.fn.tail <- function(fn.e, max.rds, max.obs.rds, blk.size, par1, par2) {
	n.blks <- ceiling((max.rds-max.obs.rds)/blk.size)
	out.blk <- rep(0, n.blks)
	for (i in 1:n.blks) {
		out.blk[i] <- sum(fn.e(((i-1)*blk.size+max.obs.rds+1):(i*blk.size+max.obs.rds), par1, par2))
	}
	tt4 <- sum(out.blk)
	return(tt4)
}

# Integrates the tail of a log-normal from max.obs.rds to max.rds
est.fn.tail.int.lnorm <- function(fn.e, max.rds, max.obs.rds, par1, par2) {
	output <- tryCatch(integrate(fn.e, mu.e=par1, sig.e=par2, lower=max.obs.rds, upper=max.rds)$value, warning=function(w){-1}, error=function(e){-1})
	if (output == -1) {
		output <- est.fn.tail(fn.e, 6e+6, max.obs.rds, 1e+6, par1, par2)
	}
	return( output )
}

# Integrates the tail of a log-logistic from max.obs.rds to max.rds
est.fn.tail.int.llog <- function(fn.e, max.rds, max.obs.rds, par1, par2) {
	output <- tryCatch(integrate(fn.e, alp=par1, gam=par2, lower=max.obs.rds, upper=max.rds)$value, warning=function(w){-1}, error=function(e){-1})
	if (output == -1) {
		output <- est.fn.tail(fn.e, 6e+6, max.obs.rds, 1e+6, par1, par2)
	}
	return( output )
}

# Integrates the tail of a weibull from max.obs.rds to max.rds
est.fn.tail.int.wb <- function(fn.e, max.rds, max.obs.rds, par1, par2) {
	output <- tryCatch(integrate(fn.e, lam=par1, gam=par2, lower=max.obs.rds, upper=max.rds)$value, warning=function(w){-1}, error=function(e){-1})
	if (output == -1) {
		output <- est.fn.tail(fn.e, 6e+6, max.obs.rds, 1e+6, par1, par2)
	}
	return( output )
}

#Fit and estimate the tail values using the log-normal distribution
est.lnorm.upd <- function(out.res.rds, sp.per.rd, wts.fit=NULL, g.st.inits=seq(.01, 20, by=0.1)) {
	#Data Frame:
	ds <- data.frame(x=out.res.rds, y=sp.per.rd); dim(ds)
	if (is.null(wts.fit)) {wts.fit <- rep(1, length(out.res.rds))}

	#Log-Normal Fit:
	rhs <- function(x, mu.e, sig.e) {
		return(1-pnorm((log(x)-mu.e)/sig.e))
	}

	alp.init <- log(mean(ds$x))
	for (med.init in alp.init) {
		for (g.st in g.st.inits) {
			x1 <- tryCatch(m.2 <- nls(y ~ rhs(x, mu.e, sig.e), data=ds, weights=wts.fit, start=list(mu.e=med.init,sig.e=g.st), trace=F, control=nls.control(maxiter=300, tol=1e-7)), warning=function(w){'warning'}, error=function(e){'error'})
			if (class(x1) != "character") {
				break
			}
		}
		if (class(x1) != "character") {
			break
		}
	}

	#Stop if there was no convergence
	if (class(x1) == 'character') {return(-1)}

	#Tail Behaviour:
	mu.est.ci <- summary(m.2)$coefficients[1,1]; sig.est.ci <- summary(m.2)$coefficients[2,1]
	out.prm <- matrix(c(mu.est.ci, sig.est.ci), nrow=1, ncol=2)

	fn.ests <- rhs(out.res.rds, mu.est.ci, sig.est.ci)
	obs.ests <- sp.per.rd; rds.ests <- out.res.rds

	#Estimate the tail value
	max.obs.rds <- max(out.res.rds)
	tail.sp.est <- est.fn.tail.int.lnorm(rhs, Inf, max.obs.rds, mu.est.ci, sig.est.ci)
	return(list(vals=tail.sp.est, prm=out.prm, tail.fn=fn.ests, tail.obs=obs.ests, tail.rds=rds.ests, fit=m.2))
}

#Estimate the tail values using the log-logistic
est.llog.upd <- function(out.res.rds, sp.per.rd, wts.fit=NULL, g.st.inits=seq(-10, 10, by=0.5)) {
	#Data Frame
	ds <- data.frame(x=out.res.rds, y=sp.per.rd); dim(ds)
	if (is.null(wts.fit)) {wts.fit <- rep(1, length(out.res.rds))}
	
	#Log-Logistic Fit:
	rhs <- function(x, alp, gam) {
		return(1/(1+alp*x^gam))
	}

	alp.init <- 1/median(ds$x)
	#for (med.init in seq(0.1, 3.5, by=0.1)*alp.init) {
	for (med.init in alp.init) {
		for (g.st in g.st.inits) {
			x1 <- tryCatch((m.2 <- nls(y ~ rhs(x, alp, gam), data=ds, weights=wts.fit, start=list(alp=med.init,gam=g.st), trace=F, control=nls.control(maxiter=300, tol=1e-7))), warning=function(w){'warning'}, error=function(e){'error'})
			#print(x1)
			if (class(x1) != "character") {
				break
			}
		}
		if (class(x1) != "character") {
			break
		}
	}
	#Stop if there was no convergence
	if (class(x1) == 'character') {return(-1)}

	#Tail Behaviour:
	mu.est.ci <- summary(m.2)$coefficients[1,1]; sig.est.ci <- summary(m.2)$coefficients[2,1]
	out.prm <- matrix(c(mu.est.ci, sig.est.ci), nrow=1, ncol=2)
	
	fn.ests <- rhs(out.res.rds, mu.est.ci, sig.est.ci)
	obs.ests <- sp.per.rd; rds.ests <- out.res.rds

	#Estimate the tail value
	max.obs.rds <- max(out.res.rds)
	tail.sp.est <- est.fn.tail.int.llog(rhs, Inf, max.obs.rds, mu.est.ci, sig.est.ci)
	return(list(vals=tail.sp.est, prm=out.prm, tail.fn=fn.ests, tail.obs=obs.ests, tail.rds=rds.ests, fit=m.2))
}

#Estimate the tail using the Weibull
est.weibull.upd <- function(out.res.rds, sp.per.rd, wts.fit=NULL, g.st.inits=seq(.01, 5, by=0.1)) {
	#Data Frame:
	ds <- data.frame(x=out.res.rds, y=sp.per.rd); dim(ds)
	if (is.null(wts.fit)) {wts.fit <- rep(1, length(out.res.rds))}
	
	#Weibull Fit:
	rhs <- function(x, lam, gam) {
		return(-lam*(x^gam))
	}
	rhs.e <- function(x, lam, gam) {
		return( exp(-lam*(x^gam)) )
	}

	alp.init <- -log(0.5)/median(ds$x)
	for (med.init in alp.init) {
		for (g.st in g.st.inits) {
			x1 <- tryCatch((m.2 <- nls(log(y) ~ rhs(x, lam, gam), data=ds, weights=wts.fit, start=list(lam=med.init,gam=g.st), trace=F, control=nls.control(maxiter=300, tol=1e-7))), warning=function(w){'warning'}, error=function(e){'error'})
			if (class(x1) != "character") {
				break
			}
		}
		if (class(x1) != "character") {
			break
		}
	}
	#Stop if there was no convergence
	if (class(x1) == 'character') {return(-1)}

	#Tail Behaviour:
	mu.est.ci <- summary(m.2)$coefficients[1,1]; sig.est.ci <- summary(m.2)$coefficients[2,1]
	out.prm <- matrix(c(mu.est.ci, sig.est.ci), nrow=1, ncol=2)
	
	fn.ests <- rhs.e(out.res.rds, mu.est.ci, sig.est.ci)
	obs.ests <- sp.per.rd; rds.ests <- out.res.rds

	max.obs.rds <- max(out.res.rds)
	tail.sp.est <- est.fn.tail.int.wb(rhs.e, Inf, max.obs.rds, mu.est.ci, sig.est.ci)
	return(list(vals=tail.sp.est, prm=out.prm, tail.fn=fn.ests, tail.obs=obs.ests, tail.rds=rds.ests, fit=m.2))
}




#########################################################################################
### END: Functions to fit the individual curves and obtain unobserved species estimates

#### Evalutation ##### (true rates are known)
## calculate the true covariance and count estimates based on the underlying rates
gen.true.vals <- function(orig.rates, ct.cut) {
	pois.prob <- sapply(orig.rates, function(x){c(dpois(0:ct.cut, lambda=x), ppois(ct.cut, lambda=x, lower=FALSE))})
	exp.cts <- apply(pois.prob, 1, sum) # expected counts from poisson probabilities
	exp.cts.tr <- exp.cts[-1]
	exp.var <- apply(pois.prob, 1, function(x){sum(x*(1-x))}) # poisson var

	## variance calculated from poisson probabilities
	out.cov.array <- array(0, dim=c(dim(pois.prob)[1], dim(pois.prob)[1], dim(pois.prob)[2]))
	for (i in 1:dim(pois.prob)[2]) {
		out.cov.array[,,i] <- est.mv.cov(pois.prob[,i])
	}
	out.cov.array <- apply(out.cov.array, c(1, 2), sum)
	out.cov.array <- out.cov.array + t(out.cov.array)
	diag(out.cov.array) <- diag(out.cov.array)/2
	return(list(cts=exp.cts.tr, ctsall=exp.cts, cov=out.cov.array))
}



#### Regularization
#Regularization - linear interpolation for p
scale.spc.est <- function(est.vec.orig, est.vec.reg, perc.var) {
	return(est.vec.orig*(1-perc.var) + est.vec.reg*perc.var)
}

#Regularization - The percentage change in the estimates after regularizing by p
calc.perc.change <- function(est.vec.orig, est.vec.reg, perc.var, obs.spc) {
	scaled.est <- scale.spc.est(est.vec.orig, est.vec.reg, perc.var) - obs.spc
	unreg.est <- est.vec.orig - obs.spc
	est.change <- (scaled.est - unreg.est)/unreg.est
	return(est.change)
}

#Regularization sensitivity calculation
calc.perc.cutoff <- function(est.vec.orig, est.vec.reg, obs.spc, est.cutoff, perc.cutoff, step.size=0.005) {
	perc.var <- perc.var.start <- 0
	m.diff <- 1
	while (1) {
		perc.var <- perc.var + step.size
		est.change <- calc.perc.change(est.vec.orig, est.vec.reg, perc.var, obs.spc)
		m.diff <- mean(abs(est.change) <= est.cutoff)
		#print(m.diff); print(perc.var); print(mean(abs(est.change)))
		if (m.diff <= perc.cutoff) {
			break
		} else {
			perc.var.start <- perc.var
		}
	}
	return(perc.var.start)
}



## Apply data regularization ### (part of esimation)
data.cts.smooth <- function(orig.dat, ct.cut, perc.var, smooth.type='log', use.wts=TRUE) {
	### uniform smoothing
	if (smooth.type == 'linear') {
		idx.use.len <- min(ct.cut, max(which(orig.dat > 0)))
		tau.reg <- get.yk.rds(orig.dat)*perc.var/get.yk.rds(rep(1, idx.use.len))
		scale.adj <- orig.dat*perc.var
		obs.dat.proc <- orig.dat - scale.adj
		obs.dat.proc[1:idx.use.len] <- obs.dat.proc[1:idx.use.len] + tau.reg
		obs.dat.proc <- ifelse(obs.dat.proc < 0, 0, obs.dat.proc)
		obs.dat.proc <- obs.dat.proc*(get.yk.rds(orig.dat)/get.yk.rds(obs.dat.proc))
		output <- obs.dat.proc
	} else if (smooth.type == 'log') {
	### log smoothing
		pos.idx <- 1:min(ct.cut, (min(which(orig.dat == 0))-1))
		orig.dat.pos <- orig.dat[pos.idx]
		outv <- log(orig.dat.pos); predv <- 1:length(orig.dat.pos)
		if (use.wts) {
			est.ls <- lm(outv ~ predv, weights=orig.dat.pos)
		} else {
			est.ls <- lm(outv ~ predv)
		}
		obs.dat.proc <- orig.dat
		exp.scale <- exp(est.ls$fitted.values)*perc.var*get.yk.rds(obs.dat.proc[pos.idx])/get.yk.rds(exp(est.ls$fitted.values))
		obs.dat.proc[pos.idx] <- obs.dat.proc[pos.idx]*(1-perc.var) + exp.scale
		obs.dat.proc <- ifelse(obs.dat.proc < 0, 0, obs.dat.proc)
		obs.dat.proc <- obs.dat.proc*(get.yk.rds(orig.dat)/get.yk.rds(obs.dat.proc))
		output <- obs.dat.proc
	}
	if (round(get.yk.rds(orig.dat) - get.yk.rds(output), 2) != 0) {stop('Scaling Error.')}
	return(output)
}



##########################################
### Data processing, reads and scaling ###

#Given count data and a cutoff, generate the observed data frequencies
est.strz.ry <- function(dat.out, count.cut, est.s=FALSE) {
	t1 <- data.frame(x=seq(0, count.cut+1)); t2 <- as.data.frame(table(dat.out)) #only keep counts below the cutoff
	t3 <- merge(t1, t2, by.x='x', by.y='dat.out', all.x=T); t3 <- t3[order(t3$x),]; t3[which(is.na(t3[,2])),2] <- 0
	r.y <- t3$Freq #create the output dataframe
	if (est.s) {
		return(r.y[2:length(r.y)])
	} else {
		return(r.y)
	}
}


#Aggregate data above the cutoff
trunc.data <- function(dat.out, count.cut) {
	high.counts <- which(dat.out > count.cut)
	dat.out[high.counts] <- count.cut + 1 #place holder for the truncation
	return(dat.out)
}

## Scale x_i data; convert to y_k data and scale
data.scale.conv <- function(obs.dat, ct.cut, dat.scale) {
	obs.dat.proc <- round(est.strz.ry(trunc.data(obs.dat, ct.cut), ct.cut, est.s=T)*dat.scale)
	obs.dat.scale <- rep(1:length(obs.dat.proc), times=obs.dat.proc)
	return(obs.dat.scale)
}

## Scale the y_k data to have to same number of reads as orig.dat
scale.curve <- function(yk.dat, orig.dat) {
	return( yk.dat*get.yk.rds(orig.dat)/get.yk.rds(yk.dat) )
}

## Get the number of reads from y_k data
get.yk.rds <- function(orig.dat) {
	return( sum(orig.dat*(1:length(orig.dat))) )
}

## Data processing
proc.input.data <- function(obs.dat, quant.data, ct.cut, perc.var=0, reg.type='linear') {
	orig.dat <- est.strz.ry(trunc.data(obs.dat, ct.cut), ct.cut, est.s=T) ## convert to y_t data

	exp.dat <- rep(1:length(obs.dat), times=obs.dat)
	exp.dat <- exp.dat[sample(1:length(exp.dat))]

	exp.dat.cut <- exp.dat[1:as.integer((length(exp.dat)*quant.data))]
	orig.dat.cut <- est.strz.ry(trunc.data(table(exp.dat.cut), ct.cut), ct.cut, est.s=T) ## the yk counts
	if (perc.var != 0) {
		orig.dat.cut <- data.cts.smooth(orig.dat.cut, ct.cut, perc.var, reg.type)
	}
	orig.dat.trc <- orig.dat.cut; orig.dat.trc[length(orig.dat.trc)] <- 0 ## the yk counts with yk > ct.cut removed
	
	n.spc <- round(get.yk.rds(orig.dat.trc)) ## the number of reads available for the truncated data
	return(list(orig=orig.dat, origcut=orig.dat.cut, origtrc=orig.dat.trc, rds=exp.dat, rdscut=exp.dat.cut, spc=n.spc))
}

proc.yk.data <- function(orig.dat, ct.cut, perc.var=0, reg.type='linear', use.wts=TRUE) {
	orig.dat.cut <- orig.dat
	if (perc.var != 0) {
		orig.dat.cut <- data.cts.smooth(orig.dat.cut, ct.cut, perc.var, reg.type, use.wts)
	}
	orig.dat.trc <- orig.dat.cut; orig.dat.trc[length(orig.dat.trc)] <- 0 ## the yk counts with yk > ct.cut removed
	
	n.spc <- round(get.yk.rds(orig.dat.trc)) ## the number of reads available for the truncated data
	return(list(orig=orig.dat, origcut=orig.dat.cut, origtrc=orig.dat.trc, spc=n.spc))
}




#######################
#### Bootstrapping ####

## Generate boostratp iterates from y_k format orig.dat
gen.boot.iterates <- function(orig.dat, n.split, run.type, boot.iter.fn=NA, save.data=FALSE) {
	n.split.run <- as.integer(n.split)
	if (is.na(boot.iter.fn) || !file.exists(boot.iter.fn)) {
		if (run.type == 'boot') {
			data.split.array <- array(NA, dim=c(length(orig.dat), 2, n.split.run))
			for (split.idx in 1:n.split.run) {
				exp.dat <- rep(1:length(orig.dat), times=orig.dat)
				idx.keep <- sample(1:length(exp.dat), length(exp.dat), replace=T)
				data.split.array[,1,split.idx] <- est.strz.ry(exp.dat[idx.keep], ct.cut, est.s=TRUE)
				data.split.array[,1,split.idx] <- scale.curve(data.split.array[,1,split.idx], orig.dat)
				data.split.array[,2,split.idx] <- est.strz.ry(exp.dat[setdiff(1:length(exp.dat), unique(idx.keep))], ct.cut, est.s=TRUE)
				data.split.array[,2,split.idx] <- scale.curve(data.split.array[,2,split.idx], orig.dat)
			}
		}
		if (run.type == 'bootv2') {
			data.split.array <- array(NA, dim=c(length(orig.dat), 2, n.split))
			for (split.idx in 1:n.split.run) {
				exp.dat <- rep(1:length(orig.dat), times=orig.dat)
				idx.keep <- sample(1:length(exp.dat), length(exp.dat), replace=T)
				data.split.array[,1,split.idx] <- est.strz.ry(exp.dat[idx.keep], ct.cut, est.s=TRUE)
				data.split.array[,2,split.idx] <- est.strz.ry(exp.dat[setdiff(1:length(exp.dat), unique(idx.keep))], ct.cut, est.s=TRUE)*length(exp.dat)/(length(exp.dat) - length(unique(idx.keep)))
			}
		}
		if (run.type == 'multi') {
			data.split.array <- array(NA, dim=c(length(orig.dat), 2, n.split.run))
			for (split.idx in 1:n.split.run) {
				data.split.array[,1,split.idx] <- c(rmultinom(1, size=sum(orig.dat), prob=orig.dat/sum(orig.dat)))
				data.split.array[,2,split.idx] <- c(rmultinom(1, size=sum(orig.dat), prob=orig.dat/sum(orig.dat)))
			}
		}
		if (save.data) {
			save(list=c('data.split.array'), file=boot.iter.fn)
		}
	} else {
		load(boot.iter.fn)
	}
	return(list(dat=data.split.array, nrun=n.split.run))
}

## Calculate L-2 distances between two sets
calc.split.dist <- function(data.split.array, n.split.run, orig.dat, run.type) {
	if (run.type %in% c('boot', 'bootv2', 'multi')) {
		dist.split <- sapply(1:n.split.run, function(x) {mean((data.split.array[,1,x] - data.split.array[,2,x])^2)})
	} else {
		dist.split <- sapply(1:n.split.run, function(x) {mean((data.split.array[,1,x] - orig.dat)^2)})
	}
	return(dist.split)
}



####################
#### Estimation ####

## Generate the species difference curve ##
gen.spc.diff <- function(orig.dat.trc, n.spc=NA) {
	if (is.na(n.spc)) {
		n.spc <- round(get.yk.rds(orig.dat.trc))
	}
	out.rd.edat <- 2:n.spc
	sp.rd.edat <- s.det(out.rd.edat, orig.dat.trc) - s.det(out.rd.edat-1, orig.dat.trc)
	curve.est.dat <- cbind(out.rd.edat, sp.rd.edat)
	return(curve.est.dat)
}

est.curves.spc <- function(obs.dat.proc, quant.curve.list, max.obs=Inf) {
	orig.dat <- obs.dat.proc$orig
	orig.dat.cut <- obs.dat.proc$origcut
	orig.dat.trc <- obs.dat.proc$origtrc
	n.spc <- obs.dat.proc$spc

	out.rd.edat <- 2:n.spc
	sp.rd.edat <- s.det(out.rd.edat, orig.dat.trc) - s.det(out.rd.edat-1, orig.dat.trc)

	if (length(out.rd.edat) > max.obs) {
		cut.idx <- unique(round(seq(1, length(out.rd.edat), length=max.obs)))
		out.rd.edat <- out.rd.edat[cut.idx]
		sp.rd.edat <- sp.rd.edat[cut.idx]
	}
	curve.est.dat <- cbind(out.rd.edat, sp.rd.edat)

	## store the fit for each index
	curve.fit.list <- list()
	dat.est.names <- c('obs', 'wtd.est', 'lnorm', 'llog', 'weib')
	dat.store <- array(NA, dim=c(length(quant.curve.list), length(dat.est.names)), dimnames=list(quant.curve.list, dat.est.names))

	prm.est.names <- c('ln1', 'ln2', 'll1', 'll2', 'wb1', 'wb2')
	prm.store <- array(NA, dim=c(length(quant.curve.list), length(prm.est.names)), dimnames=list(quant.curve.list, prm.est.names))

	wt.nm <- c('wt1.ln', 'wt2.ll', 'wt3.wb')
	wts.store <- array(NA, dim=c(length(wt.nm), length(quant.curve.list)), dimnames=list(wt.nm, quant.curve.list))

	for (quant.curve.idx in 1:length(quant.curve.list)) {
		quant.curve <- quant.curve.list[quant.curve.idx]
		
		## set the prop.of the curve to use (orig. function remove first obs)
		fit.idx <- which(out.rd.edat >= quantile(out.rd.edat, quant.curve)) 
		rds.fit <- out.rd.edat[fit.idx]; spc.fit <- sp.rd.edat[fit.idx]
		wts.fit <- (1:length(rds.fit))^2; wts.fit <- wts.fit/sum(wts.fit)*length(wts.fit)

		fit.est1 <- est.lnorm.upd(rds.fit, spc.fit, wts.fit, g.st.inits=seq(.01, 20, by=0.05))
		fit.est2 <- est.llog.upd(rds.fit, spc.fit, wts.fit, g.st.inits=seq(-5, 5, by=0.1))
		fit.est3 <- est.weibull.upd(rds.fit, spc.fit, wts.fit, g.st.inits=seq(.01, 5, by=0.1))

		param.vec <- c(fit.est1$prm, fit.est2$prm, fit.est3$prm) #order: lnorm, llog, weibull
		est.vals <- c(fit.est1$vals, fit.est2$vals, fit.est3$vals)
		est.vals.total <- c(fit.est1$vals, fit.est2$vals, fit.est3$vals)+sum(orig.dat)

		### Perform the weighted optimization ###
		x0 <- c(0, 0, 0)
		lbs <- rep(0, length(x0)); ubs <- rep(1, length(x0))
		local_opts <- list( "algorithm" = "NLOPT_LD_LBFGS", "xtol_rel"  = 1.0e-8, "maxeval" = 10000 )
		opts <- list( "algorithm" = "NLOPT_LD_AUGLAG", "xtol_rel" = 1.0e-8, "maxeval" = 10000, "local_opts"  = local_opts, "print_level" = 0 )

		obs.spc <- sum(orig.dat)
		weib.rest.true <- fit.est3$vals/fit.est1$vals*(fit.est1$vals/(obs.spc + fit.est1$vals))
		llog.rest.true <- fit.est1$vals/fit.est2$vals*(obs.spc/(obs.spc + fit.est1$vals))

		llog.rest.min <- min(weib.rest.true, llog.rest.true) #mdl2
		#llog.rest.avg <- min(llog.rest.true, mean(c(weib.rest.true, llog.rest.true)))

		weib.rat <- (fit.est1$vals/(obs.spc + fit.est1$vals))
		llog.rat <- (obs.spc/(obs.spc + fit.est1$vals))
		
		ln.obs.adj <- ifelse(weib.rat > llog.rat, mean(c(fit.est1$vals, fit.est2$vals)), mean(c(fit.est1$vals, fit.est3$vals)))
		weib.rat2 <- (ln.obs.adj/(obs.spc + ln.obs.adj))
		llog.rat2 <- (obs.spc/(obs.spc + ln.obs.adj))
		weib.rat <- mean(c(weib.rat, weib.rat2))
		llog.rat <- mean(c(llog.rat, llog.rat2))

		ubs.mod3.sel <- c(max(c(weib.rat, llog.rat)), llog.rest.min, weib.rest.true)
		ubs.mod3.sel[1] <- max(c(ubs.mod3.sel[1], 1 - sum(ubs.mod3.sel[c(2, 3)]) + 1e-6))
		ubs.mod3.sel <- ifelse(ubs.mod3.sel > 1, 1, ubs.mod3.sel)
		ubs.mod3.sel <- ifelse(ubs.mod3.sel < 0, 0, ubs.mod3.sel)
		res4mod3 <- nloptr( x0=x0, eval_f=eval_f0_mix_pen1, eval_grad_f=eval_f0_mix_pen1_grad, 
			lb = lbs, ub = ubs.mod3.sel, opts = opts, 
			rds.fit=rds.fit, spc.fit=spc.fit, wts.fit=wts.fit, param.vec=param.vec)
			wts.mod3 <- res4mod3$solution; unobs.mod3 <- wts.mod3%*%est.vals

		out.res.vec <- c(sum(orig.dat),	unobs.mod3 + sum(orig.dat),	est.vals.total)
		curve.fit.list[[quant.curve.idx]] <- list(fit.est1, fit.est2, fit.est3)
		dat.store[quant.curve.idx,] <- out.res.vec
		prm.store[quant.curve.idx,] <- param.vec
		wts.store[,quant.curve.idx] <- rbind(wts.mod3)
	}
	return(list(res=dat.store, prm=prm.store, wts=wts.store, curve=curve.est.dat, fit=curve.fit.list))
}

### Point Estimate Secondary Functions ###
s.det <- function(x.val, orig.dat, n.spc=NA) {
	if (is.na(n.spc)) {
		n.spc <- round(sum(orig.dat*(1:length(orig.dat)))) #total reads/observations
	}
	n.spc.obs <- 1:length(orig.dat) #yk = no. obs k times k = 1, ... 
	output <- sapply(x.val, function(y) {sum(( (1 - ( (1 - y/n.spc)^(n.spc.obs) ))*orig.dat ))*(as.numeric(y <= n.spc))})
	return(output)
}

s.p.det.cont <- function(x.val) {
	n.spc.obs <- 1:length(orig.dat)
	output <- sum( n.spc.obs*orig.dat/n.spc*( (1 - x.val/n.spc)^(n.spc.obs-1) ) )
	return(output)	
}

s.p.det.disc <- function(x.val) {
	n.spc.obs <- 1:length(orig.dat)
	output <- sum( -orig.dat*sapply(n.spc.obs, function(n) {sum( choose(n, 0:(n-1))*(1 - x.val/n.spc)^(0:(n-1))*(-1/n.spc)^(n-(0:(n-1))) )}) )
	return(output)
}

lnorm <- function(x, mu.e, sig.e) {
	return(1-pnorm((log(x)-mu.e)/sig.e))
}

llog <- function(x, alp, gam) {
	return(1/(1+alp*x^gam))
}

weib <- function(x, lam, gam) {
	return( exp(-lam*(x^gam)) )
}

## functions to find the optimal wights between the three curves
eval_f0_mix_pen1 <- function(x, rds.fit, spc.fit, wts.fit, param.vec) {
	y.lnorm <- 1-pnorm((log(rds.fit)-param.vec[1])/param.vec[2])
	y.llog <- 1/(1+param.vec[3]*rds.fit^param.vec[4])
	y.weib <- exp(-param.vec[5]*(rds.fit^param.vec[6]))
	return(sum( (((spc.fit - (x[1]*y.lnorm + x[2]*y.llog + x[3]*y.weib))^2) + (1/spc.fit)*x[2]*x[3])*wts.fit ))
}

eval_f0_mix_pen1_grad <- function(x, rds.fit, spc.fit, wts.fit, param.vec) {
	y.lnorm <- 1-pnorm((log(rds.fit)-param.vec[1])/param.vec[2])
	y.llog <- 1/(1+param.vec[3]*rds.fit^param.vec[4])
	y.weib <- exp(-param.vec[5]*(rds.fit^param.vec[6]))
	diff.sum <- (-2)*(spc.fit - (x[1]*y.lnorm + x[2]*y.llog + x[3]*y.weib))*wts.fit
	return( c(sum(diff.sum*y.lnorm), sum(diff.sum*y.llog + (1/spc.fit)*x[3]*wts.fit), sum(diff.sum*y.weib + (1/spc.fit)*x[2]*wts.fit)) )
}


fit.fn <- function(rds.fit, param.vec, wt.est) {
	y.lnorm <- 1-pnorm((log(rds.fit)-param.vec[1])/param.vec[2])
	y.llog <- 1/(1+param.vec[3]*rds.fit^param.vec[4])
	y.weib <- exp(-param.vec[5]*(rds.fit^param.vec[6]))
	return(cbind(y.lnorm, y.llog, y.weib)%*%wt.est)
}



