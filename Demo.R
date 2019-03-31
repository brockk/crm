
library(dfcrm)

skeleton = c(0.1, 0.15, 0.25, 0.35, 0.5)
num_doses <- length(skeleton)
dose_indices <- 1:num_doses

library(ggplot2)
dose_data <- data.frame(
  Index = dose_indices,
  PriorProbDLT = skeleton
)


# Plot the prior curve.

## Using "base R" graphics
plot(dose_indices, skeleton, type = 'b', ylim = c(0,1))
# Meaningful but perhaps a bit ugly

## Using ggplot2
ggplot(dose_data, aes(x = Index, y = PriorProbDLT)) + 
  geom_point() + 
  geom_line() + 
  ylim(0, 1) + 
  labs(x = 'Dose-level', y = 'Prob(DLT)', title = 'Prior beliefs')

