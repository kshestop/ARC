
#Pass in the parameters + weights to sample the rates from predefined distributions
#fix.wts = set the weight proportions exactly if true, sample from a multinomial if false
gen.rand.rates <- function(n.obs, param.vec, wts.vec, sim.type='gamma', fix.wts=TRUE) {
	if (sim.type == 'gamma') {
		dat.out <- list()
		mix.cts <- rmultinom(1, size=n.obs, prob=wts.vec) #sample from the multinomial
		for (i in 1:length(wts.vec)) {
			if (fix.wts) {
				n.block <- ceiling(n.obs*wts.vec[i])
			} else {
				n.block <- mix.cts[i,1] #randomized mixtures
			}
			if (param.vec[i,2] == 0) {
				dat.out[[i]] <- rep(0, n.block) #if zeros are included
			} else {
				dat.out[[i]] <- rgamma(n.block, shape=param.vec[i,1], rate=param.vec[i,2]) #generate gamma rates (no scaling by rds)
			}
		}
		exc.val <- length(unlist(dat.out))-n.obs
		if (exc.val > 0) {
			dat.out <- trim.list(dat.out, exc.val)
		}
	}
	if (length(unlist(dat.out)) != n.obs) {
		stop('mismatch in observations')
	}
	return(list(wts.vec=wts.vec, param.vec=param.vec, sim.type=sim.type, dat=dat.out))
}

#Given the rates, generate the random counts (rates are given as qi*T_bar)
#norm.sc = normalize the reads to be mean=1, otherwise use input
gen.rand.counts <- function(rand.rate, ind.rds, norm.sc=TRUE) {
	if (norm.sc) {
		rds.avg <- ind.rds/mean(ind.rds) #avg. reads scale
	} else {
		rds.avg <- ind.rds #do not normalize
	}
	rand.all <- unlist(rand.rate$dat) #extract the rates
	if (length(rds.avg == 1)) {
		rds.avg <- rep(rds.avg, length(rand.all)) #if reads are all the same, replicate
	}
	counts.out <- rep(0, length(rand.all))
	for (i in 1:length(rand.all)) {
		counts.out[i] <- rpois(1, rand.all[i]*rds.avg[i]) #generate counts
	}
	return(counts.out)
}

#Generate a set of random reads and standardize to have mu = 1
gen.rand.reads <- function(n.obs, min.rds, max.rds) {
	output <- as.integer(runif(n.obs, min.rds, max.rds))
	return( output/mean(output) )
}

#Helper Function: Remove excess obsevations from a list of vectors (1 each iteratively) if present
trim.list <- function(dat.out, exc.val) {
	counter <- exc.val
	while (counter != 0) {
		for (i in 1:length(dat.out)) {
			if (length(dat.out[[i]] > 0)) {
				dat.out[[i]] <- dat.out[[i]][1:(length(dat.out[[i]])-1)]
				counter <- counter-1
			}
			if (counter == 0) {
				break
			}
		}
	}
	return(dat.out)
}


