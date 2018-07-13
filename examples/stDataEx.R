\dontrun{
  
library(ggplot2)

data("coprecip")

# view response for 1982
plot(coprecip, t = 1982)

# view local covariate for 1982
plot(coprecip, t = 1982, type = 'covariate')

# view sea surface temperatures for 1982
plot(coprecip, t = 1982, type = 'remote')

# view sea surface EOFs (a.k.a. principal components)
# visually, the 1982 sea surface temperatures appear to have a strong 
# El Nino pattern
plot(coprecip, type = 'eof', pattern = 1)

# we can look at the scores of the EOFs over time to confirm this
plot(coprecip, type = 'eof_scores', pattern = 1) + 
  geom_vline(xintercept = 1982, lty=3)

# the first EOF explains <40% of variability in sea surface temperatures
plot(coprecip, type = 'eof_scree')

}
