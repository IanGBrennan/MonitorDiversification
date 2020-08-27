#### 1. Loading the functions and packages ####

# chose your working directory : put the path of the downloaded Rcode folder below
setwd("~/Google.Drive/R.Analyses/ClaDS/")
library(TESS)

# loading the needed functions
source("likelihood_ClaDS0.R")
source("sim_ClaDS.R")
source("utils.R")
source("proposal.R")
source("run_ClaDS0.R")
source("ClaDS1_likelihood_functions.R")
source("ClaDS2_likelihood_functions.R")
source("fit_ClaDS.R")

# this one needs to be compiled ; for this you need to open a terminal, go to your 
# working directory, and enter
# R CMD SHLIB diversif_lognormal.c
dyn.load("diversif_lognormal.so")



#### 2. Simulating trees from the model #####

set.seed(1)

obj= sim_ClaDS( lamb_par=0.1,      # initial speciation rate
                mu_par=0.5,          # turnover rate (extinction/speciation)
                sigma=0.3,           # standard deviation of the new rates law
                lamb_shift=0.9,      # trend parameter alpha
                condition="taxa",    # the stoppping condition to use (can also be time)
                taxa.stop = 20,      # the number of tips to simulate (if condition=="taxa")
                prune.extinct = T)   # should extincti taxa be pruned from the result? (default to T)


# the function returns a list with the tree and the associated rates. Here is how you get to them :
tree = obj$tree
speciation_rates = obj$lamb[obj$rates]
extinction_rates = obj$mu[obj$rates]

# and the plotting function
plot.with.rate(tree,speciation_rates,
               log=T,        # should the rates be plotted on a log scale?
               lwd=3)        # width of the branches




#### 3. Infering the model parameters on a tree, for the pure birth version (ClaDS0) ####

# this function runs the mcmc until the gelman stopping criterion (that takes several chains and 
# measures the difference between the variance within and among chains) is reached. It should be below 
# 1.05, the functions prints the gelman factor every "iteration" iterations. The function reach 3 chains 
# at the same time, that can be ran in parallel (see nCPU)

# here I run it on the tree simulated above, which is very small

sampler.varanus = run_ClaDS0(tree=vtree,        # the data
              name="Monitor_ClaDS0.Rdata",          # the name the results will be saved on (if NULL it won't be)
              nCPU=3,             # the number of CPUs to use (3 chains are run so it only makes sense to make it 1 or 3)               
              pamhLocalName = "local",   # the function is writing in a text file to make the execution quicker, this is the name of this file
              iteration=500000,   # number of iterations after which the gelman factor is computed and printed. The function stops if it is below 1.05
              thin=2000,           # number of iterations after which the chains state is recorded
              update=1000, adaptation=5000)  # options for the initial proposal adaptation phase



#### 4. Plotting the results (ClaDS0) ####

# if the exemple in 3. takes too long to run, it is saved in the attached file "example.Rdata"

load("example_ClaDS0.Rdata")
sampler_ClaDS0=sampler

ntips=tree$Nnode+1
nedges=2*tree$Nnode
npar=nedges+3


# prune the first iterations of the chains
chains=mcmc.list(lapply(1:3,function(j){mcmc(sampler[[j]]$chain[-(1:10),])}))
plot(mcmc.list(lapply(1:3,function(j){mcmc(sampler[[j]]$chain[-(1:10),c(1:3,nedges+4)])})))




# the mcmc above is run on the log of the relative rates (rate/(parentRate*alpha)), 
# below I compute the log of the absolute rates within the chains to be then able to get
# to the MAPs (posterior maximas) of the rates
for(k in 1:length(chains)){
  for(l in 1:nrow(chains[[k]])){
    chains[[k]][l,4:npar]=get_rates(tree,c(chains[[k]][l,3],
                                            chains[[k]][l,4:npar]+chains[[k]][l,2]))[-1]
  }
}

# extract the MAPs (argmax of the posterior marginals)
MAPS=sapply(1:npar, function(i){D=density(c(chains[[1]][,i],
                                            chains[[2]][,i],
                                            chains[[3]][,i]))
return(D$x[which.max(D$y)])})
MAPS[2:npar]=exp(MAPS[2:npar])

print(paste0("sigma = ",MAPS[1]," ; alpha = ",MAPS[2]," ; lambda_0 = ", MAPS[3]))

plot.with.rate(tree, MAPS[4:npar],log=T,lwd=3)

# with the true simulated rates on the left panel :
plot.with.rate(tree, speciation_rates, MAPS[4:npar],log=T,lwd=3, same.scale = T)

# on such a small trees the inference is not very good. It gets better on larger trees.



#### 5. Infering the model parameters (constant turnover, ClaDS2) ####

# prepare the sampler
sampler.v = prepare_ClaDS(tree=vtree,          # the phylogeny
                        sample_fraction=(76/83),  # the fraction of species of the clade present in the reconstructed phylogeny
                        Nchain = 3,         # number of MCMC chains
                        nlambda = 1000,     # number of lambdaspace steps in likelihood function
                        nt = 30,            # number of time steps in likelihood function
                        model_id="ClaDS2")  # use model_id = "ClaDS1" for constant extinction rate
                        #res_ClaDS0 = sampler_ClaDS0)  # the output of ClaDS0 to use as a startpoint. If NULL (the default)a random startpoint is used 
# and run it
sampled.v = fit_ClaDS(           
  mcmcSampler = sampler.v,       # an object as returned by prepare ClaDS or fit_ClaDS
  iterations = 10000,            # number of iterations
  thin = 2000,                   # number of iterations after which the chains state is recorded
  nCPU = 3)                    # number of cores to use (between 1 and the number of chains)


# you can continue to run a fit
sampler = fit_ClaDS(mcmcSampler = sampler,iterations = 100,thin = 50, nCPU = 3)

# and plot the results
plotChains_ClaDS(sampled.v)
MAPS=getMAPS_ClaDS(sampled.v)
plot.with.rate(vtree,speciation_rates,exp(MAPS[-(1:4)]),log=T,lwd=3)  # it should, of course, be run for many more iterations

plot.with.rate(vtree, exp(MAPS[-(1:4)]), log=T, lwd=4)

dtree <- extract.clade(vtree, 85)
plot(dtree)

plot.with.rate(dtree,exp(MAPS[-(1:4)]), log=T, lwd=4)


library(RPANDA)
c2.varanus <- fit_ClaDS(tree = vtree,
                           sample_fraction = (76/83),
                           iterations = 1000000,
                           thin = 1000,
                           #it_save = 10000,
                           file_name = "C2_Varanus",
                           model_id = "ClaDS2",
                           nCPU = 3)
plot_ClaDS_chains(c2.varanus)
maps <- getMAPS_ClaDS(c2.varanus, thin = 10)
print(paste0("sigma = ", maps[1], " ; alpha = ",
             maps[2], " ; epsilon = ", maps[3], " ; l_0 = ", maps[4] ))

# plot the inferred branch specific speciation rates
plot_ClaDS_phylo(vtree, maps[-(1:4)], lwd = 4, 
                 log=T, show.tip.label = T)

plot(vtree); edgelabels(bg = "white", fr="c", cex=0.4)

extract.ClaDS.rates <- function(rates, min.branch, max.branch){
  clade.rates <- rates[(min.branch + 4):(max.branch + 4)]
  return(clade.rates)
}

odatria.rates <- extract.ClaDS.rates(maps, 7, 49)
hapturo.rates <- extract.ClaDS.rates(maps, 69, 83)
euprep.rates  <- extract.ClaDS.rates(maps, 83, 101)
varanus.rates <- extract.ClaDS.rates(maps, 50, 66)
sotero.rates  <- extract.ClaDS.rates(maps, 119, 137)
empagus.rates <- extract.ClaDS.rates(maps, 112, 118)
african.rates <- extract.ClaDS.rates(maps, 139, 149)

all.rates <- data.frame(sp_rates = odatria.rates, group = "Odatria")
all.rates <- add_row(all.rates, sp_rates = hapturo.rates, group = "Hapturosaurus")
all.rates <- add_row(all.rates, sp_rates = euprep.rates , group = "Euprepiosaurus")
all.rates <- add_row(all.rates, sp_rates = varanus.rates, group = "Varanus_Papusaurus")
all.rates <- add_row(all.rates, sp_rates = sotero.rates , group = "Soterosaurus")
all.rates <- add_row(all.rates, sp_rates = empagus.rates, group = "Empagusia")
all.rates <- add_row(all.rates, sp_rates = african.rates, group = "Polydaedalus_Psammosaurus")

library(ggridges)
library(ggplot2)

ggplot(all.rates, aes(x=sp_rates, y=group, 
                       point_color=sp_rates)) +
  geom_density_ridges(jittered_points = T, scale = 2,
                      point_shape = "|", alpha = 0.75,
                      point_size = 2, size = 0.5,
                      position = position_points_jitter(height=0),
                      aes(fill = group)) +
  scale_fill_brewer(palette = "RdYlBu") +
  theme_bw()




par(mfrow=c(1,2)); densityplot(maps); densityplot(odatria.rates)

data("Caprimulgidae_ClaDS2")

Caprimulgidae_ClaDS2$sampler
# extract the Maxima A Posteriori for each parameter
maps = getMAPS_ClaDS(Caprimulgidae_ClaDS2$sampler, thin = 1)

print(paste0("sigma = ", maps[1], " ; alpha = ",
             maps[2], " ; epsilon = ", maps[3], " ; l_0 = ", maps[4] ))



plot_ClaDs2 <- function (phylo, rates, rates2 = NULL, same.scale = T, main = NULL, 
          lwd = 2, log = T, show.tip.label = F, ...) 
{
  Colors = colorRampPalette(c("steelblue2", "paleturquoise3", 
                              "palegreen2", "yellow2", "salmon1", "darkorange", "red", 
                              "red4"))(100)
  if (is.null(rates2)) {
    if (log) 
      rates = log(rates)
    if (isTRUE(all.equal(rep(as.numeric(rates[1]), length(rates)), 
                         as.numeric(rates)))) {
      col = rep(1, length(rates))
      plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
           main = main, edge.width = lwd, ...)
      if (log) {
        image.plot(z = c(exp(rates[1]), 2 * exp(rates[1])), 
                   col = Colors, horizontal = T, legend.only = T)
      }
      else {
        image.plot(z = c(rates[1], 2 * rates[1]), col = Colors, 
                   horizontal = T, legend.only = T)
      }
    }
    else {
      col = round((rates - min(rates))/diff(range(rates)) * 
                    99) + 1
      plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
           main = main, edge.width = lwd, ...)
      if (log) {
        min = min(rates)
        max = max(rates)
        m10 = floor(min/log(10))
        M10 = ceiling(max/log(10))
        if ((M10 - m10) < 4) {
          ticks = c(1, 2, 5)
        }
        else {
          ticks = 1
        }
        ticks = as.vector(sapply(m10:M10, function(k) {
          return(ticks * 10^k)
        }))
        lt = length(ticks[ticks > exp(min) & ticks < 
                            exp(max)])
        if (lt < 4) {
          ticks = c(round(exp(min), max(0, -1 * m10 + 
                                          (lt < 2))), ticks, round(exp(max), max(0, 
                                                                                 -1 * M10 + 1 + (lt < 2))))
        }
        image.plot(z = c(min, max), col = Colors, horizontal = T, 
                   legend.only = T, axis.args = list(at = log(ticks), 
                                                     labels = ticks))
      }
      else {
        image.plot(z = as.matrix(rates), col = Colors, 
                   horizontal = T, legend.only = T)
      }
    }
  }
  else {
    if (log) {
      rates = log(rates)
      rates2 = log(rates2)
    }
    if (same.scale) {
      min = min(min(rates), min(rates2))
      max = max(max(rates), max(rates2))
      par(mfrow = c(1, 2))
      col = round(((rates - min)/(max - min)) * 99) + 1
      plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
           edge.width = lwd, ...)
      col = round(((rates2 - min)/(max - min)) * 99) + 
        1
      plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
           edge.width = lwd, ...)
      par(mfrow = c(1, 1))
      if (log) {
        m10 = floor(min/log(10))
        M10 = ceiling(max/log(10))
        if ((M10 - m10) < 4) {
          ticks = c(1, 2, 5)
        }
        else {
          ticks = 1
        }
        ticks = as.vector(sapply(m10:M10, function(k) {
          return(ticks * 10^k)
        }))
        lt = length(ticks[ticks > exp(min) & ticks < 
                            exp(max)])
        if (lt < 4) {
          ticks = c(round(exp(min), max(0, -1 * m10 + 
                                          (lt < 2))), ticks, round(exp(max), max(0, 
                                                                                 -1 * M10 + 1 + (lt < 2))))
        }
        image.plot(z = c(min, max), col = Colors, horizontal = T, 
                   legend.only = T, axis.args = list(at = log(ticks), 
                                                     labels = ticks))
      }
      else {
        image.plot(z = c(min, max), col = Colors, horizontal = T, 
                   legend.only = T)
      }
    }
    else {
      par(mfrow = c(1, 2))
      if (isTRUE(all.equal(rep(rates[1], length(rates)), 
                           rates))) {
        col = rep(1, length(rates))
        plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
             edge.width = lwd, ...)
        if (log) {
          image.plot(z = c(exp(rates[1]), 2 * exp(rates[1])), 
                     col = Colors, horizontal = T, legend.only = T)
        }
        else {
          image.plot(z = c(rates[1], 2 * rates[1]), col = Colors, 
                     horizontal = T, legend.only = T)
        }
      }
      else {
        col = round(((rates - min(rates))/(max(rates) - 
                                             min(rates))) * 99) + 1
        plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
             edge.width = lwd, ...)
        if (log) {
          min = min(rates)
          max = max(rates)
          m10 = floor(min/log(10))
          M10 = ceiling(max/log(10))
          if ((M10 - m10) < 4) {
            ticks = c(1, 2, 5)
          }
          else {
            ticks = 1
          }
          ticks = as.vector(sapply(m10:M10, function(k) {
            return(ticks * 10^k)
          }))
          lt = length(ticks[ticks > exp(min) & ticks < 
                              exp(max)])
          if (lt < 4) {
            ticks = c(round(exp(min), max(0, -1 * m10 + 
                                            (lt < 2))), ticks, round(exp(max), max(0, 
                                                                                   -1 * M10 + 1 + (lt < 2))))
          }
          image.plot(z = c(min, max), col = Colors, horizontal = T, 
                     legend.only = T, axis.args = list(at = log(ticks), 
                                                       labels = ticks))
        }
        else {
          image.plot(z = as.matrix(rates), col = Colors, 
                     horizontal = T, legend.only = T)
        }
      }
      if (isTRUE(all.equal(rep(rates2[1], length(rates2)), 
                           rates2))) {
        col = rep(1, length(rates2))
        plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
             edge.width = lwd, ...)
        if (log) {
          image.plot(z = c(exp(rates2[1]), 2 * exp(rates2[1])), 
                     col = Colors, horizontal = T, legend.only = T)
        }
        else {
          image.plot(z = c(rates2[1], 2 * rates2[1]), 
                     col = Colors, horizontal = T, legend.only = T)
        }
      }
      else {
        col = round(((rates2 - min(rates2))/(max(rates2) - 
                                               min(rates2))) * 99) + 1
        plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
             edge.width = lwd, ...)
        if (log) {
          min = min(rates2)
          max = max(rates2)
          m10 = floor(min/log(10))
          M10 = ceiling(max/log(10))
          if ((M10 - m10) < 4) {
            ticks = c(1, 2, 5)
          }
          else {
            ticks = 1
          }
          ticks = as.vector(sapply(m10:M10, function(k) {
            return(ticks * 10^k)
          }))
          lt = length(ticks[ticks > exp(min) & ticks < 
                              exp(max)])
          if (lt < 4) {
            ticks = c(round(exp(min), max(0, -1 * m10 + 
                                            (lt < 2))), ticks, round(exp(max), max(0, 
                                                                                   -1 * M10 + 1 + (lt < 2))))
          }
          image.plot(z = c(min, max), col = Colors, horizontal = T, 
                     legend.only = T, axis.args = list(at = log(ticks), 
                                                       labels = ticks))
        }
        else {
          image.plot(z = as.matrix(rates2), col = Colors, 
                     horizontal = T, legend.only = T)
        }
      }
    }
    par(mfrow = c(1, 1))
  }
}





