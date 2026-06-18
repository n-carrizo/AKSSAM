# AKSSAM: Automatic Knot Selection in Smooth Additive Models
Repository associated with the paper:

> Carrizosa, N., Durban, M., Guerrero, V. (2026) `Automatic Knot Selection in Smooth Additive Models`. Manuscript to submit for publication.

This repository contains the implementation in `R` of AKSSAM, which is an algorithm designed to perform automatic knot selection in generalized additive models (GAMs).
Additionally, the scripts and data required to reproduce the computational experiments presented in the article are included.

## The main function

```r
AKSSAM(X, y, family, lambda.init, ndx, bdeg, 
       maxiter1, maxiter2, maxiter3, tol1, tol2, tol3, 
       epsilon)
```

**Description:** Performs automatic knot selection in GAMs via an alternating implementation of the adaptive ridge and Fellner-Schall algorithm.

**Input:**
<ul>
  <li><code>X</code>: Matrix with observed variables for each covariate (columns)</li>
  <li><code>y</code>: Vector of observed values of the objective variable</li>
  <li><code>family</code>: String; either 'gaussian', 'poisson', or 'binomial', indicating the distribution of the objective variable</li>
  <li><code>lambda.init</code>: Vector of penalization parameters used as initialization </li>
  <li><code>ndx</code>: Vector of inner intervals for each covariate </li>
  <li><code>bdeg</code>: Vector containing the degree of the B-spline basis for each covariate  </li>
  <li><code>maxiter1</code>: Maximum iterations for weight optimization </li>
  <li><code>maxiter2</code>: Maximum iterations for penalization optimizatio*n </li>
  <li><code>maxiter3</code>: Maximum iterations for IRLS </li>
  <li><code>tol1</code>: Relative tolerance for convergence in weight optimization </li>
  <li><code>tol2</code>: Relative tolerance for penalization optimization </li>
  <li><code>tol3</code>: Relative tolerance for IRLS </li>
  <li><code>epsilon</code>: Epsilon term for the adaptive ridge procedure </li>
</ul>


**Output:** A list with the following elements:
<ul>
  <li><code>lambda</code>: Vector of penalizations maximizing the restricted log-likelihood </li>
  <li><code>K_sel</code>: Matrix with selected knots for each covariate (columns) for optimal penalization </li>
  <li><code>New_Design_list</code>: List of design matrices after knot selection for each covariate and intercept corresponding to the optimal penalization </li>
  <li><code>alpha.new</code>: Ordered vector (intercept first, then covariates) storing resulting parameter estimates for the optimal penalization </li>
</ul>

## Example

You may find below an use example of AKSSAM in a simulated Binomial GAM scenario and an illustration of the estimated model.
```r
## Function for data simulation
simulate.data.binom = function(m, n){

  set.seed(12345)                      # Set seed for reproducibility
  
  X = matrix(1, n, m)                  # Initialize the instance mat.
  FF = matrix(1, n, m)                 # Initialize the func. values mat.
  f_total = rep(0, n)                  # Linear estimator
  
  for (i in 1:m) {
    x <- runif(n, 0, 1)
    X[, i] <- scale(x, center = TRUE, scale = FALSE)
    
    if (i == 1) {
      f <- sin(5.5 * pi * x) + 3 * x - 5 * x^2 + cos(10 * pi * x)
    } else {
      f <- 2.5* sin(pi^2 * x)
    }
    
    f_centered <- f - mean(f)
    FF[, i] <- f_centered
    f_total <- f_total + f_centered
  }
  
  
  p <- plogis(f_total)                 # Binomial probabilities
  y <- rbinom(n, 1, p)                 # Response

  
  list(X = X,
       y = y,
       FF = FF)
}


## Simulate the data
m = 2                           # Nº of covariates
n = 1000                        # Nº of instances
ll = simulate.data.binom(m,n)   # Simulate the data
X = ll$X                        # Covariate samples
y = ll$y                        # Objective variable
FF = ll$FF                      # Real additive terms values 
rm(ll)                          # Erase the list


## AKSSAM 
model = AKSSAM(X = X, y = y, family = "binomial", 
               lambda.init = rep(10, m), ndx = rep(30, m), bdeg = rep(3L, m), 
               maxiter1 = 100, maxiter2 = 100, maxiter3 = 100, 
               tol1 = 1e-5, tol2 = 1e-5, tol3 = 1e-5, 
               epsilon = 1e-5)


## Obtain the estimated sets of knots and coefficients
K_sel = model$K_sel                   # Selected knots
Design_list = model$New_Design_list   # Design matrices
alpha = model$alpha.new               # Coefficients
B_star <- do.call(cbind, Design_list) # Joint design matrix B*
predictions <- B_star %*% alpha       # Predictions in X
 
 
## Make predictions and estimated effects
# Initialize terms
Design_list = vector("list", m + 1)
Design_list[[1]] = matrix(1, 100, 1) 
basis_length = rep(0,m)
covariate_predictions = vector("list", m)
bdeg = rep(3, m)
# Iterate through covariates
for (t in 1:m){
  # Build the design matrix with selected knots
  B = my.bbase4(seq(-0.5,0.5, length.out = 100), K_sel[[t]], bdeg[t])
  Design_list[[t+1]] = B
  basis_length[t] = dim(B)[2]
  if (t == 1){index_0 = 2; index_1 = 1 + basis_length[1]
  }else{index_0 = 2 + sum(basis_length[1:(t-1)]); index_1 = index_0 + basis_length[t] - 1}
  # Effect of the t-th covariate
  covariate_predictions[[t]] = B %*% alpha[index_0:index_1]
}
# Entire prediction
B_star <- do.call(cbind, Design_list) # Joint design matrix B*
predictions <- B_star %*% alpha     
```

<div align="center">
<img src="Images/frontpage.png" width="1200"/>
</div>

## Repository contents

The structure of the repository is the following:
```
/AKSSAM
├── /Computational Experiments 
│   ├── /Required_Scripts              # Folder containing the neccesary functions and datasets
│   │
│   ├── RealData1ElectricLoad.qmd      # Script for the ElectricLoad (Gaussian) real data scenario
│   ├── real_data1.RData               # Results from executing RealData1ElectricLoad.qmd
│   ├── RealData2PimaIndians.qmd       # Script for the PimaIndians (Binomial) real data scenario
│   ├── real_data2.RData               # Results from executing RealData2PimaIndians.qmd
│   ├── Simulation1Gaussian.qmd        # Script for the Gaussian GAM simulation scenario
│   ├── simulation1.RData              # Results from executing Simulation1Gaussian.qmd
│   ├── Simulation2Poisson.qmd         # Script for the Poisson GAM simulation scenario
│   └── simulation2.RData              # Results from executing Simulation2Poisson.qmd
│
├── AKSSAM.R                           # R implementation of the AKSSAM method
│
├── EXAMPLE.qmd                        # Illustrative examples of AKSSAM
|
├── /Images                            # Miscellaneous folder
│
├── LICENSE.R                          # License information
│
└── README.md                          # Project documentation
```

## Contact

This project was developed by *Nicolás Carrizosa* (https://github.com/n-carrizo) as part of a research project within the Universidad Carlos III de Madrid.

It benefited from the support of the grant PRTR-CNS2023, funded by MICIU/AEI /10.13039/501100011033 and by European Union NextGenerationEU/PRTR and is part of the project *Modelos Aditivos con Restricciones para la Optimización Global*.

## Disclaimer

This is a preliminary version of the algorithm. It is unstable and may exhibit issues.

## License

This repository is licensed under the [MIT License](LICENSE).