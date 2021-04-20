library(dplyr)
library(ggplot2)

setwd("/Applications/PyRate/pyrate_mcmc_logs/Varanus_Set_2/Plates")
log.files <- dir(getwd(), "mcmc.log")


mcmc.res <- NULL
for(k in 1:length(log.files)){
  # read in the log file
  curr.log <- read.table(log.files[k], header=T)
  
  # remove a certain percentage as burn-in (here 25%)
  curr.log <- curr.log[(nrow(curr.log)*0.25):nrow(curr.log),]
  
  # extract mean values for all the parameters
  means <- apply(curr.log, 2, mean)
  
  # put the information into a dataframe
  mcmc.res <- rbind(mcmc.res, data.frame(Value = means, Parameter = names(means),
                         Type = "Plates", Rep = k))
}
rownames(mcmc.res) <- NULL


#mcmc.constant <- mcmc.res
mcmc.plates <- mcmc.res

all.mcmc <- rbind(mcmc.constant, mcmc.plates)

posterior <- dplyr::filter(all.mcmc, Parameter == "posterior")

like.diff <- NULL
for (j in 1:length(unique(posterior$Rep))){
  curr.rep <- filter(posterior, Rep == j)
  diff <- filter(curr.rep, Type == "Constant")$Value - filter(curr.rep, Type == "Plates")$Value
  if(diff > 0){diff.color = "#3d98d3"}
  if(diff < 0){diff.color = "#8cdb5e"}
  like.diff <- rbind(like.diff, data.frame(Value = diff, Color = diff.color, Rep = j))
}

ggplot(like.diff, aes(x=Rep, y=Value, label=round(Value, 1))) +
  geom_point(stat = "identity", size=10, color=like.diff$Color) +
  geom_text(color="white", size=3) +
  scale_y_continuous(name="likelihood Difference (Constant Rate - 2 Rates)", limits=c(-20,20), breaks=c(-20,-10,0,10,20)) +
  scale_x_continuous(name="PyRate Dataset Replicate", breaks=c(1:7)) +
  geom_hline(yintercept=0) +
  theme_classic() + coord_flip()



