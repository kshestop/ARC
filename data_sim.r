


res.save.dir <- './'
source('./data_gen.r')

#Generate the list of scenarios:
out.scen <- expand.grid(seq(0, 1, by=0.20), seq(0, 1, by=0.20), seq(0, 1, by=0.20), seq(0, 1, by=0.20))
out.scen <- out.scen[apply(out.scen, 1, sum) == 1,]
out.scen <- out.scen[out.scen[,1] > 0 | out.scen[,2] > 0,]
out.scen <- out.scen[out.scen[,4] <= 0.4,]
out.scen <- out.scen[order(out.scen[,1], out.scen[,2], out.scen[,3], decreasing=T),]
rownames(out.scen) <- 1:dim(out.scen)[1]
param.vec <- cbind(c(1, 1, 2, 3), c(2, 1, 1, 1))

zero.pr.obs <- c(2/3, 1/2, 1/4, 1/8)
k1 <- as.matrix(out.scen)%*%zero.pr.obs
scen.zer.ord <- order(k1, decreasing=T)

#Parameters
n.obs <- 1000 #number of ### species ###
n.rds <- 1 #number of reads
n.gen <- 500 #number of times to generate data for a scenario
res.type <- 'data_v2'

ii.run.list <- 1:dim(out.scen)[1]
jj.run.list <- 1:n.gen

##########################################
### Generate data for the iterate here ###
for (ii in ii.run.list) {
	data.filename <- paste(res.save.dir, 'simulated_species_data_scenario_upd_', res.type, '_', ii, '_', n.gen, '.RData', sep='')
	if (file.exists(data.filename)) {
		load(data.filename)
	} else {
		rand.dat.storage <- array(0, dim=c(n.obs, n.gen))
		rand.rate.storage <- list() #simulated rates
		for (jj in 1:n.gen) {
			#Generate the Data:
			rand.rate <- gen.rand.rates(n.obs, param.vec, out.scen[ii,], sim.type='gamma', fix.wts=FALSE)
			dat.out <- gen.rand.counts(rand.rate, 1, norm.sc=FALSE) #dat.out = obs. counts, ind.rds = total reads
			rand.rate.storage[[jj]] <- rand.rate
			rand.dat.storage[,jj] <- dat.out #store the simulated data for the iteration
		}
		save(list=c('rand.rate.storage', 'rand.dat.storage'), file=data.filename)
	}
}







