# Plot Description
## Gaussian Fits
- In PBR submission, every individual participant is fit to two models: a two Gaussian model and a three Gaussian model
- Each participant is plotted in this folder
  - Participants are separated by group (auto for automatic and man for manual) and are ordered by the size of the mean of the third Gaussian representing estimated mean below-threshold adjustment size
  - AIC difference between the two and three Gaussian model is also provided

## Performance over task
- Every participant's first and third block are plotted in the pdf Every-PPT-Adjustment-Over-Task
- On the x axis is trial (333 total per block)
- On the y axis is the Bernoulli parameter
- The dotted line represents the true Bernoulli parameter used to generate samples of red/blue dots
- The coloured line represents the participant's estimate of the Bernoulli parameter
- The colour of the line represents whether the adjustment to the estimate was sample-consistent (green) sample inconsistent (red) or no change (yellow)
- Participants are ordered by number of adjustments made over the course of te task (999 trials) which is given in the heading of each plot