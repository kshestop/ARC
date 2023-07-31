
#given a vector computed the paired probabilities for the covariance matrix
est.mv.cov <- function(x) {
	cov.mat <- array(0, dim=c(length(x), length(x)))
	for (i in 1:dim(cov.mat)[1]) {
		for (j in i:dim(cov.mat)[1]) {
			if (i == j) {
				cov.mat[i,j] <- x[i]*(1-x[i])
			} else {
				cov.mat[i,j] <- -x[i]*x[j]
			}
		}
	}
	return(cov.mat)
}

eval_f_z <- function(x, orig.dat, ct.cut, pois.tr.mat, s.pen) {
	trp.trc.mat <- t(pois.tr.mat)[2:dim(pois.tr.mat)[2],] #transposed + (zero) truncated transition matrix (1, ..., c, c+) x (l1, ..., ln)
	bt.cts.est <- (trp.trc.mat%*%x)[,1]
	### penalize each by the sq. deviation from adjacent values
	x.idx <- 2:(length(x)-1)
	### l^2 penalty
	obj.pen <- s.pen * sum(x^2)
	grad.pen <- s.pen * 2*x

	objc <- sum((orig.dat - bt.cts.est)^2) + obj.pen
	gradc <- sapply(1:dim(trp.trc.mat)[2], function(x){ (-2)*sum((orig.dat - bt.cts.est)*trp.trc.mat[,x]) }) + grad.pen
	return( list("objective"=objc, "gradient"=gradc) )
}

internal.cov.loop <- function(s.pen.idx, s.pen.list, n.split, data.split.array, run.type, orig.dat, zvals, ct.cut, pois.tr.mat, eval_f_z) {
	out.serr <- array(NA, dim=c(n.split))
	out.est.store <- array(NA, dim=c(n.split, length(zvals)))

	s.pen <- s.pen.list[s.pen.idx]; #print(s.pen)
	x0 <- rep(0.5, length(zvals)); 

	opt.st.t <- Sys.time()
	for (split.idx in 1:n.split) {
		data.input <- data.split.array[,1,split.idx]
		lbs <- rep(0, length(zvals)); ubs <- rep(Inf, length(zvals))
		res5 <- nloptr( x0=x0, eval_f=eval_f_z, lb = lbs, ub = ubs,
			opts = list("algorithm"="NLOPT_LD_LBFGS", "xtol_rel"=1.0e-6, "maxeval"=10000),
			orig.dat = data.input, ct.cut = ct.cut, pois.tr.mat = pois.tr.mat, s.pen = s.pen ); #print(res4)
		#print(Sys.time() - opt.st.t)
		x0 <- res5$solution

		orig.dat.z <- res5$solution
		cts.est.bt <- t(pois.tr.mat)[2:(ct.cut+2),]%*%res5$solution; cts.est.bt

		if (run.type %in% c('boot', 'bootv2', 'multi')) {
			out.serr[split.idx] <- mean((cts.est.bt - data.split.array[,2,split.idx])^2)
		} else {
			out.serr[split.idx] <- mean((cts.est.bt - orig.dat)^2)
		}
		out.est.store[split.idx,] <- orig.dat.z
	}

	####################################
	## Compared to the resampled data ##
	opt.st.t <- Sys.time()
	lbs <- rep(0, length(zvals)); ubs <- rep(Inf, length(zvals))
	res5 <- nloptr( x0=x0, eval_f=eval_f_z, lb = lbs, ub = ubs,
		opts = list("algorithm"="NLOPT_LD_LBFGS", "xtol_rel"=1.0e-6, "maxeval"=10000),
		orig.dat = orig.dat, ct.cut = ct.cut, pois.tr.mat = pois.tr.mat, s.pen = s.pen ); #print(res4)
	#print(Sys.time() - opt.st.t)

	cts.est.orig <- t(pois.tr.mat)[2:(ct.cut+2),]%*%res5$solution; cts.est.orig
	return(list(outstat=out.serr, outest=out.est.store, zsoln=res5$solution))
}

## remove the first (zero) column and row when necc.
rem.zero.mat <- function(x.mat) {
	return(x.mat[2:dim(x.mat)[1],2:dim(x.mat)[2]])
}

#### the estimated - based on the z-values
est.cov.mat.z <- function(pois.tr.mat, orig.dat.z) {
	resamp.prob <- pois.tr.mat ### (0, 1, ..., ct.cut, ct.cut+)
	est.cov.mat <- array(0, dim=c(dim(resamp.prob)[2], dim(resamp.prob)[2]))
	for (i in 1:dim(est.cov.mat)[1]) {
		for (j in i:dim(est.cov.mat)[1]) {
			if (i == j) {
				est.cov.mat[i,j] <- orig.dat.z%*%(resamp.prob[,i]*(1-resamp.prob[,i]))
			} else {
				est.cov.mat[i,j] <- orig.dat.z%*%(-resamp.prob[,i]*resamp.prob[,j])
			}
		}
	}
	est.cov.mat <- est.cov.mat + t(est.cov.mat)
	diag(est.cov.mat) <- diag(est.cov.mat)/2
	est.cov.mat <- est.cov.mat[2:dim(est.cov.mat)[1],2:dim(est.cov.mat)[2]]
	return(est.cov.mat)
}

opt.res.calc <- function(opt.pen.idx, res.z.store, pois.tr.mat, out.cov.array, ct.cut, exp.cts.tr) {
	res.out.penalty <- c('no.reg', 'nz.pen.2', 'z.orig.min', 'deciles', 'diff.z.incr', 'decile.avg')
	res.out.meas <- c('mse.cts', 'mse.var', 'mse.mat')
	res.out <- array(NA, dim=c(length(opt.pen.idx), length(res.out.meas)))
	cov.est.store <- cov.true.store <- array(NA, dim=c(ct.cut+1, ct.cut+1, length(res.out.penalty)))
	cts.est.store <- cts.true.store <- array(NA, dim=c(ct.cut+1, length(res.out.penalty)))
	rownames(res.out) <- res.out.penalty; colnames(res.out) <- res.out.meas

	for (opt.idx.sel in 1:length(opt.pen.idx)) {
		pen.idx.use <- opt.pen.idx[opt.idx.sel]
		orig.dat.z <- res.z.store[,pen.idx.use]

		est.cov.mat <- est.cov.mat.z(pois.tr.mat, orig.dat.z)
		out.cov.array.cut <- rem.zero.mat(out.cov.array)

		#### estimated
		iter.cts.est <- t(pois.tr.mat)[2:(ct.cut+2),]%*%orig.dat.z
		iter.var.est <- diag(est.cov.mat)

		#### the true (by rate)
		iter.cts.true <- exp.cts.tr
		iter.var.true <- diag(out.cov.array.cut)

		cov.est.store[,,opt.idx.sel] <- est.cov.mat
		cov.true.store[,,opt.idx.sel] <- out.cov.array.cut
		cts.est.store[,opt.idx.sel] <- iter.cts.est
		cts.true.store[,opt.idx.sel] <- iter.cts.true

		res.out[opt.idx.sel,1] <- mean((iter.cts.est - iter.cts.true)^2)
		res.out[opt.idx.sel,2] <- mean((iter.var.est - iter.var.true)^2)
		res.out[opt.idx.sel,3] <- mean((est.cov.mat - out.cov.array.cut)^2)
	}
	return(list(cove=cov.est.store, covt=cov.true.store, ctse=cts.est.store, ctst=cts.true.store, res=res.out))
}

############### Select the optimal penalty index from s.pen.list
get.opt.pen.idx <- function(s.pen.list, output.ests, res.z.store, n.split, n.quant, output.stats) {
	n.obs.quant <- as.integer(n.split/n.quant)
	out.diff.init <- rep(NA, length(s.pen.list))
	for (j in 1:length(s.pen.list)) {
		out.diff.init[j] <- mean(apply(output.ests[,,j], 1, function(x){(x - res.z.store[,1])^2})) #compared to the original
	}
	out.avg.quant.err <- apply(output.stats, 2, function(x) {sapply(seq(1, dim(output.stats)[1], by=n.obs.quant), function(z){mean(x[z:min(length(x),(z+n.obs.quant))])})})

	#### five options for bootstrap penalization
	# 0. Baseline case (no regularization)
	# 1. take the smallest non-zero value
	opt.pen.idx1 <- 2

	# 2. Take the value closest to the original estimate, but stabilized
	opt.pen.idx2 <- which.min(out.diff.init)
	if (opt.pen.idx2 == 1) {
		opt.pen.idx2 <- 2
	}
	# 3. Take the first non-one value based on the deciles.
	min.by.quant <- apply(out.avg.quant.err, 1, function(x){which.min(x)})
	min.nz <- which(min.by.quant > 1)
	if (length(min.nz) == 0) {
		opt.pen.idx3 <- 2
	} else {
		opt.pen.idx3 <- min.by.quant[min.nz[1]]
	}
	# 4. Take the first increase
	opt.pen.idx4 <- which(diff(out.diff.init) > 0)[1]
	if (is.na(opt.pen.idx4) | length(opt.pen.idx4) == 0 | opt.pen.idx4 == 1) {
		opt.pen.idx4 <- 2
	}
	# 5. Average the bootstrap quantiles
	min.by.quant <- apply(out.avg.quant.err, 1, function(x){which.min(x)})
	min.by.quant[which(min.by.quant == 1)] <- 2
	opt.pen.idx5 <- round(mean(min.by.quant, na.rm=T))

	opt.pen.idx <- c(1, opt.pen.idx1, opt.pen.idx2, opt.pen.idx3, opt.pen.idx4, opt.pen.idx5)
	return(opt.pen.idx)
}



