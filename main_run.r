## Code to generate a point estimate using the ARC method ##

## Main parameters
scn.idx <- 1 ## the simulation scenario to load
idx.tag <- 1 ## the simulation iterated to use
run.type <- 'bootv2' ## for the CI estimate

#args = commandArgs(trailingOnly=TRUE)
#scn.idx <- as.numeric(args[1])
#run.type <- as.character(args[2])
#idx.tag <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
print(paste('Scenario:', scn.idx, 'Iterate:', idx.tag))

library(nloptr)
library(MASS)
library(mvtnorm)
library(tools)
library(gtools)
options(width=120)

####################
output.dir <- './'
data.dir <- './'

source('./data_gen.r') #data generation
source('./var_est.r') ### variance est
source('./point_est.r') ### point est

####################
## Point Estimate ##

## Main parameters
ver.run <- 'test'

## Fitting parameters
ct.cut.list <- c(20, 30, 50) ## the list of truncations to use
quant.data.list <- c(0.01, 0.05, 0.10, 0.20, 0.35, 0.75, 1.0) ## the list of proportion of reads to use
quant.curve.list <- c(0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30) ## the list of proportion of the acc. curve to sue
max.obs <- 1000
reg.prm <- 0.00 #(1 percent)
boot.point.type <- 'multi' ## bootstrap for the point estimate (if used)

## ci est and bootstrap prms
n.split <- 1500 # number of bootstrap iterates
n.quant <- 10 #number of quantiles to use for splitting the bootstrap

cov.use.idx <- 4 #the covariance estimate to use (bootstrap) - alternatively #6 (avg. over bts. quants)
n.ci.sim <- 2000 #no. of ci simulations

## for the robustness checks
perc.var.list <- c(0.0, 0.5, 1.0)
reg.type <- 'log'; use.wts <- TRUE

## placeholder parameters
quant.data.idx <- length(quant.data.list)
quant.curve.idx <- 1
ct.cut.idx <- 1

quant.data <- quant.data.list[quant.data.idx]
quant.curve <- quant.curve.list[quant.curve.idx]
ct.cut <- ct.cut.list[ct.cut.idx]

# Step 1: Load the Data
if (scn.idx %in% c(1:43)) {
	sim.dat.fn <- paste0(data.dir, '/simulated/simulated_species_data_scenario_upd_data_v2_', scn.idx, '_500.RData')
	load(sim.dat.fn)
	obs.dat <- rand.dat.storage[,idx.tag]; orig.rates <- unlist(rand.rate.storage[[idx.tag]]$dat)
	obs.dat.proc <- proc.input.data(obs.dat, quant.data, ct.cut)
} else if (scn.idx == 1011) {
	sim.dat.fn <- paste0(data.dir, '/steam_catalogue_cleaned_example_2.RData')
	load(sim.dat.fn)
	obs.dat <- vect.list[[idx.tag]]; orig.rates <- NA
	obs.dat.proc <- proc.input.data(obs.dat, quant.data, ct.cut)
} else if (scn.idx == 1005) {
	sim.dat.fn <- paste0(data.dir, '/wagner_tree_microbiome_cleaned_example_1.RData')
	load(sim.dat.fn)
	obs.dat <- ldat.ex[idx.tag,]; orig.rates <- NA
	obs.dat.proc <- proc.input.data(obs.dat, quant.data, ct.cut)
}
# Step 2: Specify a particular dataset to use
#t.st <- Sys.time()
curve.pt.est <- est.curves.spc(obs.dat.proc, quant.curve.list, max.obs=max.obs)
curve.point.est <- apply(curve.pt.est$res, 2, mean)
#Sys.time() - t.st

###################################
## Estimate the confidence interval
# parameters for the ci
orig.dat <- obs.dat.proc$origcut

iter.true.vals <- gen.true.vals(orig.rates, ct.cut)
exp.cts.tr <- iter.true.vals$cts
out.cov.array <- iter.true.vals$cov

### Specify the parameters of the CI fit
ct.zval <- min(ct.cut, max(which(orig.dat != 0))+2) ## the lambdas for the z-values (underlying dist'ns)
zvals <- exp(seq(-6, log(ct.zval), length=40))
pois.tr.mat <- t(sapply(zvals, function(x) {c(dpois(0:ct.cut, lambda=x), ppois(ct.cut, lambda=x, lower=FALSE))}))
s.pen.list <- c(0, exp(seq(-2.5, 4.5, length=30)))/sum(orig.dat)

#### Generate the bootstrap iterates
boot.iter.fn <- paste0(output.dir, '/boot_iter_df_sample_scn', scn.idx, '_iter', idx.tag, '_run', run.type, '_nsplit', n.split, '_v', ver.run, '.RData')

data.split.array <- gen.boot.iterates(orig.dat, n.split, run.type, boot.iter.fn)$dat

## sort the bootstrap iterates
dist.split <- calc.split.dist(data.split.array, n.split, orig.dat, run.type)
data.split.array <- data.split.array[,,order(dist.split)]

pen.res.fn <- paste0(output.dir, '/pen_results_df_sample_run_scn_', scn.idx, '_iter', idx.tag, '_run', run.type, '_ct', ct.cut, '_nsplit', n.split, '_v', ver.run, '.RData')

### Storage array for the optimal penalty estimation ###
output.ests <- array(NA, dim=c(n.split, length(zvals), length(s.pen.list)))
res.z.store <- array(NA, dim=c(length(zvals), length(s.pen.list)))
output.stats <- array(NA, dim=c(n.split, length(s.pen.list)), dimnames=list(1:n.split, 1:length(s.pen.list)))

if (file.exists(pen.res.fn)) {
	load(pen.res.fn)
} else {
	t.st <- Sys.time()
	r.pen <- list()
	for (s.pen.idx in 1:length(s.pen.list)) {
		r.pen[[s.pen.idx]] <- internal.cov.loop(s.pen.idx, s.pen.list, n.split, data.split.array, run.type, orig.dat, zvals, ct.cut, pois.tr.mat, eval_f_z)
	}
	print(Sys.time() - t.st)
	for (s.pen.idx in 1:length(s.pen.list)) {
		output.ests[,,s.pen.idx] <- r.pen[[s.pen.idx]]$outest ## z-weight for each resample iterate
		res.z.store[,s.pen.idx] <- r.pen[[s.pen.idx]]$zsoln ## the z-weight estimates for the observed data
		output.stats[,s.pen.idx] <- r.pen[[s.pen.idx]]$outstat ## error measures by iteration
	}
	#save(list=c('output.ests', 'res.z.store', 'output.stats', 'ct.cut', 's.pen.list', 'pois.tr.mat', 'orig.dat'), file=pen.res.fn)
}
opt.pen.idx <- get.opt.pen.idx(s.pen.list, output.ests, res.z.store, n.split, n.quant, output.stats)
### evaluate each penalty measure
opt.iter.eval <- opt.res.calc(opt.pen.idx, res.z.store, pois.tr.mat, out.cov.array, ct.cut, exp.cts.tr)
cov.est.store <- opt.iter.eval$cove; cov.true.store <- opt.iter.eval$covt
cts.est.store <- opt.iter.eval$ctse; cts.true.store <- opt.iter.eval$ctst
res.out <- opt.iter.eval$res
print(res.out)

sample.filt.fn <- paste0(output.dir, '/sample_var_est_results_run_', scn.idx, '_iter', idx.tag, '_run', run.type, '_nquant', n.quant, '_ct', ct.cut, '_nsplit', n.split, '_v', ver.run, '.RData')
#save(list=c('res.out', 'cov.est.store', 'cov.true.store', 'cts.est.store', 'cts.true.store', 'opt.pen.idx'), file=sample.filt.fn)

#stop('')
est.cov.mat <- cov.est.store[,,cov.use.idx]
est.cts.mat <- cts.est.store[,cov.use.idx]



###############################################################
### the ci estimate
ci.dat.sim <- mvrnorm(n.ci.sim, mu=orig.dat, Sigma=est.cov.mat)
ci.pos.sim <- round(ifelse(ci.dat.sim<0, 0, ci.dat.sim))

int.boot.loop <- function(boot.idx, boot.dat.train, boot.dat.test, perc.var.list, quant.data, quant.curve.list, ct.cut, reg.type, use.wts, max.obs) {
	if (boot.idx %% 100 == 0) {print(boot.idx)}
	res.out.arr <- array(NA, dim=c(length(perc.var.list), 4)); rownames(res.out.arr) <- perc.var.list
	for (dat.scale.idx in 1:length(perc.var.list)) {
		perc.var <- perc.var.list[dat.scale.idx]
		obs.dat.boot <- proc.yk.data(boot.dat.train, ct.cut, perc.var, reg.type, use.wts)
		curve.pt.est <- est.curves.spc(obs.dat.boot, quant.curve.list, max.obs=max.obs)
		res.out.arr[dat.scale.idx,] <- apply(curve.pt.est$res, 2, mean)[c(2:5)]
	}
	return(list(estres=res.out.arr))
}

ci.res.fn <- paste0(output.dir, '/ci_results_df_sample_run_scn_', scn.idx, '_iter', idx.tag, '_run', run.type, '_ct', ct.cut, '_nsplit', n.ci.sim, '_reg', reg.prm, '_cov', cov.use.idx, '_v', ver.run, '.RData')

if (file.exists(ci.res.fn)) {
	load(ci.res.fn)
} else {
	#Generate the results
	t.st <- Sys.time()
	r.boot <- list()
	for (boot.idx in 1:n.ci.sim) {
		r.boot[[boot.idx]] <- int.boot.loop(boot.idx, ci.pos.sim[boot.idx,], ci.pos.sim[boot.idx,], perc.var.list, quant.data, quant.curve.list, ct.cut, reg.type, use.wts, max.obs)
	}
	print(Sys.time() - t.st)

	#Store the results
	est.st <- est.av.st <- array(NA, dim=c(dim(r.boot[[1]]$estres)[1], dim(r.boot[[1]]$estres)[2], n.ci.sim), dimnames=list(perc.var.list, 1:4, 1:n.ci.sim))
	for (boot.idx in 1:n.ci.sim) {
		est.st[,,boot.idx] <- r.boot[[boot.idx]]$estres
		res.out.arr <- r.boot[[boot.idx]]$estres
		est.av.st[,,boot.idx] <- sapply(res.out.arr[which(perc.var.list == 0),], function(x){x*(1-perc.var.list)}) + sapply(res.out.arr[which(perc.var.list == 1),], function(x){x*perc.var.list})
	}
	#save(list=c('est.st', 'est.av.st', 'ci.pos.sim', 'quant.data', 'ct.cut', 'quant.curve.list', 'curve.pt.est', 'curve.point.est', 'reg.prm', 'reg.type', 'use.wts', 'perc.var.list'), file=ci.res.fn)
}







