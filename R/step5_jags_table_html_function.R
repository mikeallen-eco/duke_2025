library(knitr)
library(kableExtra)

jags_table_html <- function(jags_model,
                            digits = 3,
                            full_table = TRUE,
                            caption = "JAGS Model Summary") {
  
  if (!"summary" %in% names(jags_model)) {
    stop("Input object does not look like a jagsUI model (no $summary element).")
  }
  
  # Extract summary table
  df <- as.data.frame(jags_model$summary)
  
  # Clean names to match your example
  df$Parameter <- rownames(df)
  rownames(df) <- NULL
  
  df <- df[, c("Parameter",
               "mean", "sd",
               "2.5%", "50%", "97.5%",
               "overlap0", "f", "Rhat", "n.eff")]
  
  # Round numeric columns
  num_cols <- sapply(df, is.numeric)
  df[num_cols] <- lapply(df[num_cols], round, digits = digits)
  
  # Format with kable + kableExtra
  
  df %>%
    kable(format = "html",
          caption = caption,
          escape = FALSE,
          align = "lrrrrrcccc") %>%
    kable_classic(full_width = full_table, html_font = "Arial") %>%
    scroll_box(width = "100%")
}
