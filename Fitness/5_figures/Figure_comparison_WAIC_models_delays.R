## Library
library(rstan)
library(RColorBrewer)
library(binom)

setwd('/Users/noemielefrancq/Documents/Project_fitness_Strep_Pneumo/SPneumoMobility/Fitness/fitness_clean_NL/')

data = read.csv('4_1_per_province_NVT_PCV7_PCV13_testswicth/results_waic.csv')
# data = data[c(1, 3:10),]
data$delay = c(seq(-4,4,1), 6, 8)
data$waic_true = data$waic
data$waic = data$waic - min(data$waic)

############################################################################################
setwd('/Users/noemielefrancq/Documents/Project_fitness_Strep_Pneumo/SPneumoMobility/Fitness/fitness_clean_NL/5_figures/')
pdf(width = 4.5, height = 4.5, file = "Figure_comparison_delays.pdf", onefile = T)
par(mai = c(1, 1, 0.5, 0.2))
############################################################################################
plot(data$delay, data$waic, pch = 16, col = 'black',
     main = 'Comparison models \n for different timings of fitness changes',
     ylab = 'Difference to the best WAIC',
     xlab = 'Year of the swicth in fitness \n (year since implementation of PCV7)')
abline(v = 0, lty = 2)
polygon(x = c(-20, 20, 20, -20), 
        y = c(-1, -1, 2, 2), 
        col = adjustcolor('grey60', alpha.f = 0.25),
        border = F)
polygon(x = c(-20, 20, 20, -20), 
        y = c(-1, -1, 7, 7), 
        col = adjustcolor('grey60', alpha.f = 0.1),
        border = F)
############################################################################################
dev.off()
############################################################################################

write.csv2(data, 'Results_WAIC.csv')
