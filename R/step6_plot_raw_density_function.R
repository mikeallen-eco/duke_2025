# function to plot average raw counts over time +/- 2 SE

plot_raw_density <- function(data = d, alpha = "GRSP"){
  
  if(alpha == "GRSP") spname <- "Grasshopper Sparrow"
  if(alpha == "EAME") spname <- "Eastern Meadowlark"
  if(alpha == "RWBL") spname <- "Red-winged Blackbird"
  if(alpha == "BOBO") spname <- "Bobolink"
  
  # get complete list of points to fill in zero samples
  allpts <- data %>% 
    select(field, year, pt) %>% 
    distinct()
  
  # format raw count data
  raw_plot_data <- data %>% 
    filter(species %in% alpha) %>%
    filter(dist < 101 | is.na(dist)) %>% 
    group_by(field, year, pt) %>% 
    tally() %>%
    right_join(allpts %>% select(year, field, pt) %>% distinct(), 
               by = join_by(pt, year, field)) %>%
    replace_na(list(n = 0)) %>%
    arrange(pt) %>%
    mutate(dens = n/(4*3.1415)) %>% 
    group_by(field, year) %>% 
    summarize(mean.num = mean(dens),
              sd = sd(dens),
              .groups = "drop") %>%
    mutate(npts = c(rep(10,7), rep(9,7)),
           yr = c(2010, 2012, 2013, 2018, 2019, 2024, 2025,
                  2010, 2012, 2013, 2018, 2019, 2024, 2025),
           se = sd/sqrt(npts))
  
  if(alpha == "RWBL") raw_plot_data <- raw_plot_data %>% filter(!yr %in% 2010:2013)
  
  # set upper y-axis limit based on data
  y_upper <- max(raw_plot_data$mean.num+2*raw_plot_data$se) * 1.1
  
  # plot
  raw_plot <- ggplot(raw_plot_data) +
    geom_point(aes(x = yr, y = mean.num, color = field), size = 4,
               position = position_dodge(0.5)) +
    geom_errorbar(aes(x = yr, ymin = mean.num-(2*se), ymax = mean.num+(2*se),
                      color = field),
                  width = 0,
                  position = position_dodge(0.5)) +
    geom_line(aes(x = yr, y = mean.num, color = field),
              position = position_dodge(0.5)) +
    scale_color_manual(values = c("steelblue", "darkred")) +
    scale_x_continuous(breaks = seq(2010,2024, by = 4)) +
    ylim(0,y_upper) +
    labs(y = paste0(spname, " density (count / ha)"), x = "", 
         color = "Field") +
    theme_bw() +
    theme(text = element_text(size = 13))
  
  
  ggsave(file.path("output", paste0(tolower(alpha), "_raw.png")), height = 4, width = 6, dpi = 600)

  raw_plot
}