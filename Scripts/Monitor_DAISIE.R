source("~/Google.Drive/R.Analyses/Convenient Scripts/Calculate_AICs.R")
library(DAISIE)

setwd("~/Documents/GitHub/MonitorDiversification/Data/")

vtree <- read.tree("~/Documents/GitHub/MonitorDiversification/Data/UsedTree.tre")

# extract branching times for clades of interest
ov.bt <- branching.times(extract.clade(vtree, getMRCA(vtree, c("komodoensis", "tristis"))))
names(ov.bt) <- NULL; sort(ov.bt, decreasing=T)

vk.bt <- max(nodeHeights(vtree)) - nodeheight(vtree, getParent(vtree, which(vtree$tip.label == "keithhornei")))
vd.bt <- max(nodeHeights(vtree)) - nodeheight(vtree, getParent(vtree, which(vtree$tip.label == "doreanus")))
vc.bt <- max(nodeHeights(vtree)) - nodeheight(vtree, getParent(vtree, which(vtree$tip.label == "chlorostigma")))

D.oz <- read.csv("DAISIE_Australia.csv", header=T)

oz.list    <- DAISIE_dataprep(datatable = D.oz,
                              island_age = 35, M = 30)
split.list <- DAISIE_dataprep(datatable = D.oz,
                              island_age = 35, M = 30,
                              number_clade_types = 2,
                              list_type2_clades = c("V_keithornei", "V_doreanus", "V_chlorostigma"),
                              prop_type2_pool = 0.1)

eqr <- DAISIE_ML(datalist = oz.list,
                 initparsopt = c(0.5, 0.5, 30, 0.01, 0.1), # lambdaC (clado), mu (ex), K (cap), gamma (imm), lamdaA (anagen)
                 idparsopt = 1:5,
                 ddmodel = 11,
                 parsfix = NULL, 
                 idparsfix = NULL) # started at 12:55 ended at 1:02, took ~7 min

eqr_noDD <- DAISIE_ML(datalist = oz.list,
                      initparsopt = c(0.5, 0.5, 0.01, 0.1), # lambdaC (clado), mu (ex), K (cap), gamma (imm), lamdaA (anagen)
                      idparsopt = c(1,2,4,5),
                      ddmodel = 11,
                      parsfix = Inf, 
                      idparsfix = 3) # started at 2:15 ended at 2:30, took ~15 min

two_lambda <- DAISIE_ML(datalist = split.list,
                        ddmodel = 11,
                        initparsopt = c(1, 0.1, 30, 0.005, 0.1, 
                                        1),
                        idparsopt = c(1,2,3,4,5,6),
                        parsfix = NULL,
                        idparsfix = NULL,
                        idparsnoshift = c(7,8,9,10)) # started at 1:50 ended at 1:56, took ~6 min

two_lambda_noDD <- DAISIE_ML(datalist = split.list,
                             ddmodel = 11,
                             initparsopt = c(1, 0.1, 0.005, 0.1, 
                                             1),
                             idparsopt = c(1,2,4,5,6),
                             parsfix = Inf,
                             idparsfix = 3,
                             idparsnoshift = c(7,8,9,10))


two_lambda_mu <- DAISIE_ML(datalist = split.list,
                           ddmodel = 11,
                           initparsopt = c(1, 0.1, 30, 0.005, 0.1, 
                                           1, 0.1),
                           idparsopt = c(1,2,3,4,5,6,7),
                           parsfix = NULL,
                           idparsfix = NULL,
                           idparsnoshift = c(8,9,10)) # started at 2:00 


all.models <- c("eqr", "eqr_noDD", "two_lambda", "two_lambda_mu")
DAISIE.AIC(all.models)


regions <- read.delim("geo_file_Regions.txt", sep="\t", row.names = 1)
dotTree(vtree, regions, labels=T,length=8)
phylo.heatmap(vtree, regions)




dsamp <- diamonds[sample(nrow(diamonds), 1000), ]
d <- ggplot(dsamp, aes(carat, price)) +
    geom_point(aes(colour = clarity))

d + scale_color_brewer()

head(morph.data)

ggplot(morph.data, aes(Mean_logSVL, PC1_LinearMeasurements)) +
  geom_point(aes(color = Oz_Group))

