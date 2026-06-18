# Load necessary libraries
library(mgcv)
library(pso)
source("data_simulation/functions.R")
source("module/ok_gam.R")

#start by simulating some univariate data

set.seed(101)
n <- 1000
x1 <- runif(n,0, 1)


y <-
  f2(x1) +  rnorm(n, sd = 5)
my_data <- data.frame(x1 = x1,
                      y = y)

# Define the initial gam model
initial_gam_model <-
  gam(y ~ s(x1, bs = "cr", k = 30),
      data = my_data, )

x_seq <- seq(min(x1), max(x1), length.out = 100)
y_seq <- f2(x_seq) - mean(f2(x_seq))

plot(
  initial_gam_model,
  pages = 1,
  residuals = TRUE,
  shade = TRUE,
  lwd = 3
)

# Add grid lines
abline(h = pretty(range(c(y, y_seq))),
       col = "lightgray",
       lty = "dotted")
abline(v = pretty(range(x1)),
       col = "lightgray",
       lty = "dotted")

# Add the data generating function line
lines(x_seq, y_seq, col = "red", lwd = 2, lty=2)

# Call the fit_gam_pso function ---> NO knot slection irl
result <- fit_gam_optim(
  gam_model = initial_gam_model,
  data = my_data,
  n_knots = 30,
  alpha = 1e-07,
  smoothing_method = "REML",
  max_iterations = 10000
)

plot(
  result$model,
  pages = 1,
  residuals = TRUE,
  shade = TRUE,
  lwd = 3
)
# Add grid lines
abline(h = pretty(range(c(y, y_seq))),
       col = "lightgray",
       lty = "dotted")
abline(v = pretty(range(x1)),
       col = "lightgray",
       lty = "dotted")
lines(x_seq,
      y_seq,
      col = "red",
      lwd = 2,
      lty = 2)


# My algorithm
# Required Packages
library(caret)
library(corpcor)
library(corrplot)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(gridGraphics)
library(MASS)
library(Matrix)
library(patchwork)
library(pROC)
library(pracma)
library(splines)
library(splines2)
library(tidyverse)
library(aspline)
library(xtable)
library(cowplot)

# Our functions
source("AKSSAM/Functions.R")

m = 1

model1 <- GAM.asplines3.wood2(X = as.matrix(x1), y = y, family = 'gaussian', 
                              lambda.init = c(100), ndx = c(120), bdeg = c(3), 
                              maxiter1 = 100, maxiter2 = 300, maxiter3 = 50, 
                              tol1 = 1e-5, tol2 = 1e-5, tol3 = 1e-5, 
                              epsilon = 1e-5)

# List with selected knots and parameter vector
K_sel <- model1$K_sel
alpha <- as.vector(model1$alpha.new)

Grid_list = vector("list", m)
# Vector of covariates' basis sizes
grid_basis_length = rep(0,m)

grid <- as.matrix(x_seq)

for (t in 1:m){                 # Covariates
  Grid_list[[t]] = my.bbase4(grid[,t], K_sel[[t]], 3)
  grid_basis_length[t] = dim(Grid_list[[t]])[2]
}

grid1 <- grid # Initialize matrix to store grid values of each covariate

for (t in 1:m){
  B = Grid_list[[t]] # Obtain the design matrix
  
  # Calculate the indexes of coeffs.
  if (t == 1){
    index_0 = 2
    index_1 = 1 + grid_basis_length[1]
  }else{
    index_0 = 2 + sum(grid_basis_length[1:(t-1)])
    index_1 = index_0 + grid_basis_length[t] - 1
  }
  
  # Covariate values in the grid
  grid1[ ,t] <- B %*% alpha[index_0:index_1]
}


plot(x1,y)
lines(x_seq, grid1[,1] + alpha[1], col = 'red', lwd = 2)
rug(model1$K_sel[[1]], lwd = 2)


#####################################

result <- fit_gam_pso(
  gam_model = initial_gam_model,
  data = my_data,
  n_knots = 30,
  alpha = 1e-07,
  smoothing_method = "REML",
  max_iterations = 100
)

plot(
  result$model,
  pages = 1,
  residuals = TRUE,
  shade = TRUE,
  lwd = 3,
  ylim = range(c(y + 3, y_seq + 2))
)
# Add grid lines
abline(h = pretty(range(c(y, y_seq))),
       col = "lightgray",
       lty = "dotted")
abline(v = pretty(range(x1)),
       col = "lightgray",
       lty = "dotted")
lines(x_seq,
      y_seq,
      col = "red",
      lwd = 2,
      lty = 2)


### OVERALL: They do not effectively do knot selection but teh smoothing is remarkable: we overfit a bit more
### 
### 
### OJO No hacemos selección de nodos !!!!!!!! A no ser que ndx sea suficientemente grande

