# function to plot species density over time based on hierarchical distance sampling analysis

library(dplyr)
library(ggplot2)

plot_HDS_density <- function(mod = "output/bobo_2025.rds"){
  
  # get alpha code
  alpha <- substr(mod, 8,11)
    
  # read in JAGS model object
  out <- readRDS(mod)
  
  # format data for plotting
  plot_data <- data.frame(field = rep(c("S", "K"), 7),
                          year = c(2010, 2010, 2012, 2012, 2013, 2013, 
                                   2018, 2018, 2019, 2019, 2024, 2024, 
                                   2025, 2025),
                          Nmed = c(median(out$sims.list$lam.S.Y10), median(out$sims.list$lam.K.Y10), 
                                   median(out$sims.list$lam.S.Y12), median(out$sims.list$lam.K.Y12), 
                                   median(out$sims.list$lam.S.Y13), median(out$sims.list$lam.K.Y13),
                                   median(out$sims.list$lam.S.Y18), median(out$sims.list$lam.K.Y18), 
                                   median(out$sims.list$lam.S.Y19), median(out$sims.list$lam.K.Y19), 
                                   median(out$sims.list$lam.S.Y24), median(out$sims.list$lam.K.Y24),
                                   median(out$sims.list$lam.S.Y25), median(out$sims.list$lam.K.Y25)),
                          Nq2.5 = c(quantile(out$sims.list$lam.S.Y10, 0.025), 
                                    quantile(out$sims.list$lam.K.Y10, 0.025), 
                                    quantile(out$sims.list$lam.S.Y12, 0.025), 
                                    quantile(out$sims.list$lam.K.Y12, 0.025), 
                                    quantile(out$sims.list$lam.S.Y13, 0.025), 
                                    quantile(out$sims.list$lam.K.Y13, 0.025),
                                    quantile(out$sims.list$lam.S.Y18, 0.025), 
                                    quantile(out$sims.list$lam.K.Y18, 0.025), 
                                    quantile(out$sims.list$lam.S.Y19, 0.025), 
                                    quantile(out$sims.list$lam.K.Y19, 0.025), 
                                    quantile(out$sims.list$lam.S.Y24, 0.025), 
                                    quantile(out$sims.list$lam.K.Y24, 0.025),
                                    quantile(out$sims.list$lam.S.Y25, 0.025),
                                    quantile(out$sims.list$lam.K.Y25, 0.025)),
                          Nq97.5 = c(quantile(out$sims.list$lam.S.Y10, 0.975), 
                                     quantile(out$sims.list$lam.K.Y10, 0.975), 
                                     quantile(out$sims.list$lam.S.Y12, 0.975), 
                                     quantile(out$sims.list$lam.K.Y12, 0.975), 
                                     quantile(out$sims.list$lam.S.Y13, 0.975), 
                                     quantile(out$sims.list$lam.K.Y13, 0.975),
                                     quantile(out$sims.list$lam.S.Y18, 0.975), 
                                     quantile(out$sims.list$lam.K.Y18, 0.975), 
                                     quantile(out$sims.list$lam.S.Y19, 0.975), 
                                     quantile(out$sims.list$lam.K.Y19, 0.975), 
                                     quantile(out$sims.list$lam.S.Y24, 0.975), 
                                     quantile(out$sims.list$lam.K.Y24, 0.975),           
                                     quantile(out$sims.list$lam.S.Y25, 0.975), 
                                     quantile(out$sims.list$lam.K.Y25, 0.975)),           
                          Nq10 = c(quantile(out$sims.list$lam.S.Y10, 0.1), 
                                   quantile(out$sims.list$lam.K.Y10, 0.1), 
                                   quantile(out$sims.list$lam.S.Y12, 0.1), 
                                   quantile(out$sims.list$lam.K.Y12, 0.1), 
                                   quantile(out$sims.list$lam.S.Y13, 0.1), 
                                   quantile(out$sims.list$lam.K.Y13, 0.1),
                                   quantile(out$sims.list$lam.S.Y18, 0.1), 
                                   quantile(out$sims.list$lam.K.Y18, 0.1), 
                                   quantile(out$sims.list$lam.S.Y19, 0.1), 
                                   quantile(out$sims.list$lam.K.Y19, 0.1), 
                                   quantile(out$sims.list$lam.S.Y24, 0.1), 
                                   quantile(out$sims.list$lam.K.Y24, 0.1),
                                   quantile(out$sims.list$lam.S.Y25, 0.1), 
                                   quantile(out$sims.list$lam.K.Y25, 0.1)),
                          Nq90 = c(quantile(out$sims.list$lam.S.Y10, 0.9), 
                                   quantile(out$sims.list$lam.K.Y10, 0.9), 
                                   quantile(out$sims.list$lam.S.Y12, 0.9), 
                                   quantile(out$sims.list$lam.K.Y12, 0.9), 
                                   quantile(out$sims.list$lam.S.Y13, 0.9), 
                                   quantile(out$sims.list$lam.K.Y13, 0.9),
                                   quantile(out$sims.list$lam.S.Y18, 0.9), 
                                   quantile(out$sims.list$lam.K.Y18, 0.9), 
                                   quantile(out$sims.list$lam.S.Y19, 0.9), 
                                   quantile(out$sims.list$lam.K.Y19, 0.9), 
                                   quantile(out$sims.list$lam.S.Y24, 0.9), 
                                   quantile(out$sims.list$lam.K.Y24, 0.9),
                                   quantile(out$sims.list$lam.S.Y25, 0.9), 
                                   quantile(out$sims.list$lam.K.Y25, 0.9))) %>%
    mutate(Dmed = Nmed/pi,
           Dq2.5 = Nq2.5/pi,
           Dq97.5 = Nq97.5/pi,
           Dq10 = Nq10/pi,
           Dq90 = Nq90/pi)
  
  # define upper y-axis plotting limit
  plot_y_upper <- max(plot_data$Dq97.5) + max(plot_data$Dq97.5)*0.15

  # define species name y-axis labels
  if(alpha %in% "bobo"){lab_sp <- "Bobolink"}
  if(alpha %in% "eame"){lab_sp <- "Eastern Meadowlark"}
  if(alpha %in% "grsp"){lab_sp <- "Grasshopper Sparrow"}
  
  # create the time series plot
  plot <- ggplot(plot_data) +
    geom_point(aes(x = year, y = Dmed, color = field), size = 4,
               position = position_dodge(0.5)) +
    geom_errorbar(aes(x = year, ymin = Dq2.5, ymax = Dq97.5,
                      color = field),
                  width = 0,
                  position = position_dodge(0.5)) +
    geom_errorbar(aes(x = year, ymin = Dq10, ymax = Dq90,
                      color = field),
                  linewidth = 1,
                  width = 0,
                  position = position_dodge(0.5)) +
    geom_line(aes(x = year, y = Dmed, color = field),
              position = position_dodge(0.5)) +
    scale_color_manual(values = c("steelblue", "darkred")) +
    scale_x_continuous(breaks = seq(2012,2025, by = 4)) +
    ylim(0,plot_y_upper) +
    labs(y = paste0(lab_sp, " density (no./ha)"), x = "", 
         color = "Field") +
    theme_bw() +
    theme(text = element_text(size = 13))
  
  # save the plot
  ggsave(file.path("output", paste0(alpha, "_density.png")), 
         height = 4, width = 6, dpi = 600)
  
  return(plot)
  
}
