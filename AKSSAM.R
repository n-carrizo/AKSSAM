# ============================================================
# Script Name : AKSSAM.R
# Project     : AKSSAM Algorithm Implementation
# Author      : Nicolas Carrizosa Arias
# Date        : 2026-06-15
# Description : This script implements the AKSSAM algorithm.
#               It includes a function implementing the 
#               proposed algorithm, all the lower-level 
#               functions required and an extension to the 
#               additive level of the A-Splines algorithm.
# ============================================================

### ----------------------------------------------------------------------------
### Utility Functions 
### ----------------------------------------------------------------------------

## Basic Functions -------------------------------------------------------------

# Function    : my.knots
# Description : Calculates adequate equally-spaced knots
#
# Input:
#   - x    : Vector of covariate values
#   - xl   : Left boundary of the desired knot interval
#   - xr   : Right boundary of the desired knot interval
#   - ndx  : Number of inner intervals
#   - bdeg : Degree of the B-spline basis (for extra knots)
#
# Output:
#   - Vector of adequate equally-spaced knots
# 
my.knots = function(x, xl = min(x), xr = max(x), ndx = 10, bdeg = 3){
  
  dx <- (xr - xl) / ndx
  
  return(seq(xl - bdeg * dx, xr + bdeg * dx, by = dx))
}


# Function    : my.bbase4
# Description : Constructs the design B-spline matrix 
#               given the knot vector
#
# Input:
#   - x     : Vector of covariate values
#   - knots : Vector of knots
#   - bdeg  : Degree of the B-spline basis
#
# Output:
#   - Design matrix of the B-spline basis
# 
my.bbase4 = function (x, knots, bdeg = 3){
  
  # outer.ok = TRUE to avoid numerical issues when the knot vector is large
  # (ndx = 60 onwards)
  
  return(splineDesign(knots, x, bdeg + 1, outer.ok = TRUE))
}

## Penalizations ---------------------------------------------------------------

# Function    : construct.penalizations
# Description : Creates the penalization matrices (Adaptive ridge + 
#               Identifiability) for the main loop
#
# Input:
#   - Design_list : List with the design matrices for each covariate
#   - order_diffs : Vector with the difference order for each covariate
#   - lambda      : Vector with the penalization for each covariate
#   - w           : List containing vectors of penalization weights 
#                   for each covariate
#
# Output:
#   - A block-wise constructed matrix representing the overall
#     penalization structure
# 
construct.penalizations = function(Design_list, order_diffs, lambda, w){
  
  # Nº of instances
  n = dim(Design_list[[1]])[1]
  ones = matrix(1, n, 1)
  
  # Size (including intercept) of the list
  m = length(Design_list)
  
  # Initialize penalization list
  Pen_list = vector("list", m)
  
  # Intercept penalization
  Pen_list[[1]] = matrix(0, 1, 1) 
  
  # Covariates' penalizations
  Pen_list[-1] <- lapply(2:m, function(i) {
    B = Design_list[[i]]                                    
    D = diff(diag(ncol(B)), differences = order_diffs[i-1]) # Differences matrix
    P = lambda[i-1] * crossprod(D * sqrt(w[[i-1]]))         # Adaptive ridge penalization
    PI = tcrossprod(crossprod(B, ones))                     # Identifiability penalization
    return(P + PI)                                          
  })
  
  # Unify all the penalizations block-wise and return
  return(as.matrix(bdiag(Pen_list)))
}


# Function    : construct.penalizations2
# Description : Constructs the identifiability penalty matrices
#
# Input:
#   - Design_list : List with the design matrices for each covariate
#
# Output:
#   - A block-wise constructed matrix representing the overall
#     identifiability penalization
# 
construct.penalizations2 = function(Design_list){
  
  # Nº of instances
  n = dim(Design_list[[1]])[1]
  ones = matrix(1, n, 1)
  
  # Size (including intercept) of the list
  m = length(Design_list)
  
  # Initialize penalization list
  Pen_list = vector("list", m)
  
  # Intercept penalization
  Pen_list[[1]] = matrix(0, 1, 1) 
  
  # Covariates' penalizations
  
  Pen_list[-1] <- lapply(2:m, function(i) {
    # Identifiability penalization
    B = Design_list[[i]]                           
    return(tcrossprod(crossprod(B, ones)))
  })
  
  # Unify all the penalizations block-wise and return
  return(as.matrix(bdiag(Pen_list)))
}


# Function    : construct.penalizations.deriv
# Description : Constructs the matrix Sj used in Fellner-Schall algorithm
#
# Input:
#   - Design_list : List with the design matrices for each covariate
#   - order_diffs : Vector with the difference order for each covariate
#   - w           : List containing vectors of penalization weights for each 
#                   covariate
#   - j           : Integer indicating the covariate of interest
#
# Output:
#   - A block-wise constructed matrix Sj, where Dj Wj Dj is placed 
#     in the j-th covariate's position and zeros elsewhere
# 
construct.penalizations.deriv = function(Design_list, order_diffs, w, j){
  
  # Nº of instances
  n = dim(Design_list[[1]])[1]
  ones = matrix(1, n, 1)
  
  # Size (including intercept) of the list
  m = length(Design_list)
  
  # Index j accounts for covariate, j+1 to index accounting intercept
  j = j + 1
  
  # Initialize penalization list
  Pen_list = vector("list", m)
  
  # Intercept penalization
  Pen_list[[1]] = matrix(0, 1, 1) 
  
  # Covariates' penalization (without identifiability)
  Pen_list[-1] <- lapply(2:m, function(i) {
    B = Design_list[[i]]
    if (i == j){ 
      # Only calculate block if i == j
      # Differences matrix
      D = diff(diag(ncol(B)), differences = order_diffs[i-1]) 
      
      # Weighted Smoothing Penalization without lambda_j
      t(D) %*% diag(w[[i-1]]) %*% D           
    } else{ 
      # If i != j then make the block 0
      matrix(0, dim(t(B)%*%B)[1], dim(t(B)%*%B)[2])
    }
  })
  
  # Unify all the penalizations block-wise and return
  return(as.matrix(bdiag(Pen_list)))
}


## Robust Matrix Inversion -----------------------------------------------------

# Function    : robust_solve
# Description : Solves systems of the form Ax = b, either via QR decomposition
#               or, if the former raises en error, employs SVD. If none of both
#               work, an error is raised. Alternatively, if no b is provided,
#               the inverse of A is computed.
#
# Input:
#   - A  : Matrix of real numbers
#   - b  : Vector of real numbers.
#
#
# Output:
#   - A vector with the solution/approximation to the equation Ax = b if b is 
#     provided, or the inverse of A if not. 
#
robust_solve <- function(A, b = NULL){
  result <- withCallingHandlers(
    tryCatch({
      if (is.null(b)) qr.solve(A, tol = 1e-18) else qr.solve(A, b, tol = 1e-18)
    }, error = function(e) {
      
      svd_decomp <- svd(A)
      keep <- which(abs(svd_decomp$d) > 1e-10)
      
      if (length(keep) == 0) stop("All singular values below tolerance. Matrix is effectively zero.")
      
      v_keep <- svd_decomp$v[, keep]
      u_keep <- svd_decomp$u[, keep]
      
      if (is.null(b)) {
        v_keep %*% (1 / svd_decomp$d[keep] * t(u_keep))
      } else {
        v_keep %*% (crossprod(u_keep, b) / svd_decomp$d[keep])
      }
    }),
    warning = function(w) {
      warning(w)  
      invokeRestart("muffleWarning")
    }
  )
  
  if (!is.null(b)) result <- result %>% as.vector else result <- result %>% as.matrix
  
  return(result)
}


## IRLS Algorithm --------------------------------------------------------------

# We took the following stability measures:
#   
# - Initialize via an oversmoothed and unpenalized IRLS instance to obtain an 
# initial feasible seed for the algorithm
# - Added extra +1e-6 term in potentially unstable denominators.


# Function    : IRLS.init
# Description : Fits unpenalized B-splines in a GLM context 
#               recursively using Newton-Raphson for ndx = 5.
#               Used as starting seed for IRLS.
#
# Input:
#   - X       : Matrix with observed variables for each covariate (columns)
#   - y       : Vector of observed values of the objective variable
#   - family  : String; either 'gaussian', 'poisson', or 'binomial', 
#               indicating the distribution of the objective variable
#   - maxiter : Maximum number of iterations for the IRLS algorithm
#   - tol     : Relative tolerance for convergence of the IRLS algorithm
#
# Output:
#   - Linear predictor (eta) at convergence
#                
IRLS.init = function(X, y, family, bdeg, maxiter = 50, tol = 1e-5){
  
  ## Correctly store the GLM family 
  family <- match.arg(family, choices = c("poisson", "binomial"))
  
  # Distribution family match
  if (family == "poisson") {
    g_inv <- function(x) exp(x)
    dg <- function(x) 1 / (x + 1e-6)
    V <- function(x) x
  } else if (family == "binomial") {
    g_inv <- function(x) 1 / (1 + exp(-x))
    dg <- function(x) 1 / (x * (1 - x))
    V <- function(x) x * (1 - x)
  }
  
  ## Construct the smaller basis 
  # Compute the equally-spaced knot vectors for ndx = 5
  K = vector('list', m)
  for (i in 1:m){
    K[[i]] = my.knots(X[,i],min(X[,i]), max(X[,i]), 5, bdeg[i])
  }
  
  # List with design matrices
  Design_list = vector("list", m + 1)
  # Intercept term
  Design_list[[1]] = matrix(1, n, 1)
  # Design matrices for each covariate:
  Design_list[-1] <- lapply(1:m, function(i) my.bbase4(X[,i], K[[i]], bdeg[i]))
  # Construct the added design matrix in order B^* = [1:B1:...:Bm]
  B <- do.call(cbind, Design_list)
  
  # Initialize a list containing the parameters
  par_list = vector('list', m) 
  # Adequate initialization by family
  if (family == "binomial"){
    par_list <- lapply(1:m, function(i) rep(0, ncol(Design_list[[i+1]])))
  } else{
    par_list <- lapply(1:m, function(i) rep(1, ncol(Design_list[[i+1]])))
  }
  # Adequate intercept
  if (family == 'binomial') intercept = log(mean(y)/(1 - mean(y))) else intercept = log(mean(y))
  # Parameter initialization
  par <- c(intercept, unlist(par_list))
  
  # GLM Scenario
  
  # Initialize the linear predictor
  eta_old <- B %*% par
  
  ## IRLS Algorithm with Max iteration tolerance 
  for (iter in 1:maxiter){
    
    # Calculate mu and Omega
    mu <- g_inv(as.vector(eta_old))
    Omega_vec <- V(mu)
    
    # Solve the IRLS algorithm instance
    z <- eta_old +  (y - mu) * dg(mu)
    B_weighted <- B * Omega_vec 
    par <- qr.solve(crossprod(B_weighted, B), crossprod(B_weighted, z), tol = 1e-18)
    
    # Newly estimated linear predictor
    eta_new <- B %*% par
    
    # Convergence check
    convergence <- sum((eta_old - eta_new)^2) / sum((eta_new)^2)
    if (convergence < tol){
      break 
    }
    
    # Update the linear estimate
    eta_old <- eta_new 
  }
  
  ## Output: Estimated linear predictor
  return(eta_new)
}


# Function    : Fit.BSplines.Penalized
# Description : Fits penalized B-splines in a GLM context using IRLS.
#
# Input:
#   - B        : Design matrix
#   - y        : Vector of observed values of the objective variable
#   - P        : Penalization matrix for the overall problem
#   - family   : String; either 'gaussian', 'poisson', or 'binomial', 
#                indicating the distribution of the objective variable
#   - maxiter  : Maximum number of iterations for the IRLS algorithm
#   - tol      : Relative tolerance for convergence of the IRLS algorithm
#   - eta_init : Initialization value for the linear predictor
#
# Output:
#   - par        : Ordered vector (intercept first, then covariates) 
#                  storing resulting IRLS parameter estimates
#   - Omega_vec  : Vector of weights of the IRLS-associated weight matrix.
#                  NA if family == 'gaussian'
#   - eta        : Linear predictor at convergence
#                  NA if family == 'gaussian'
#                
Fit.BSplines.Penalized = function(B, y, P, family, maxiter = 50, tol = 1e-5, eta_init = NULL){
  
  family <- match.arg(family, choices = c("gaussian", "poisson", "binomial"))
  
  if (family == "poisson") {
    g_inv <- function(x) exp(x)
    dg <- function(x) 1 / (x + 1e-6)
    V <- function(x) x
  } else if (family == "binomial") {
    g_inv <- function(x) 1 / (1 + exp(-x))
    dg <- function(x) 1 / (x * (1 - x))
    V <- function(x) x * (1 - x)
  } else if (family == "gaussian") {
    par <- robust_solve(crossprod(B) + P, crossprod(B, y))
    return(list(par = par))
  }
  
  eta_old <- eta_init
  
  for (iter in 1:maxiter){
    
    mu <- g_inv(as.vector(eta_old))
    Omega_vec <- V(mu)
    
    z <- eta_old + (y - mu) * dg(mu)
    B_weighted <- B * Omega_vec 
    par <- robust_solve(crossprod(B_weighted, B) + P, crossprod(B_weighted, z))
    
    eta_new <- B %*% par
    
    # Handle situations in which y = mu
    if (anyNA(eta_new)){
      z[which(is.na(z))] = eta_old[which(is.na(z))]
      par <- robust_solve(crossprod(B_weighted, B), crossprod(B_weighted, z))
      eta_new <- B %*% par
    }
    
    numerator <- sum((eta_old - eta_new)^2)
    denominator <- sum((eta_new)^2)
    convergence <- if (denominator < 1e-8 || numerator == Inf) Inf else numerator / denominator
    if (convergence < tol) break
    
    eta_old <- eta_new 
  }
  
  return(list(
    par = par,
    Omega_vec = Omega_vec,
    eta = eta_new
  ))
}


## Adaptive Ridge Algorithm ----------------------------------------------------

# Function    : adridge
# Description : Performs the adaptive ridge (WPSS) algorithm
#
# Input:
#   - Design_list  : List containing the design matrices for each covariate and the intercept
#   - basis_length : Vector specifying the length of the B-spline basis for each covariate
#   - order_diffs  : Vector containing the difference order for each covariate
#   - family       : String; either 'gaussian', 'poisson', or 'binomial', 
#                    indicating the distribution of the objective variable
#   - lambda       : Vector containing the penalizations for each covariate
#   - w            : List of initial weights for each covariate
#   - old_sel      : List of initial indicators of selected inner knots for each covariate
#   - par          : List of initial parameter estimates for each covariate
#   - B.new        : Full B-spline basis matrix
#   - y            : Vector of observed values of the objective variable
#   - maxiter      : Maximum iterations for the WPSS algorithm
#   - tol          : Relative tolerance for convergence of the WPSS algorithm
#
# Output:
#   - A list with the following elements:
#       * sel       : List containing indicators of selected inner knots for each covariate
#       * w         : List of updated weights for each covariate
#       * par.new   : Vector of parameter estimates for the WPSS solution
#       * converge  : Boolean indicating whether the algorithm converged within the specified tolerance
#       * eta       : Linear predictor at convergence.
#                     NA if family == 'gaussian'
# 
adridge <- function(Design_list, basis_length, family, lambda, w, old_sel,
                    B, y, bdeg, epsilon, maxiter, tol, eta_init = NULL){
  
  # Initialize the actual selected knots indexes 
  # (Empty because it is filled within the first iteration)
  m = length(Design_list) - 1
  sel <- vector('list', m) 
  
  # Indexes for each covariate's parameters
  index_start <- cumsum(c(2, head(basis_length, -1)))
  index_end <- index_start + basis_length - 1
  
  # Penalization orders for each covariate
  pen_order <- bdeg + 1
  
  # Initialize the parameter list
  par_list = vector("list", m)
  
  # Initialize the linear estimator
  eta = eta_init
  
  # IRLS Max iterations tolerance check
  for (iter in 1:maxiter){
    
    # Construction of the penalization matrix
    P = construct.penalizations(Design_list, pen_order, lambda, w)
    
    # Estimation via P-Splines of the new coefficients
    ll = Fit.BSplines.Penalized(B, y, P, family, maxiter, tol, eta)
    par.new = ll$par
    if (family != 'gaussian') eta = ll$eta
    
    # Assign the coefficients in a list and
    # update the weights and selected knots indexes
    for (i in 1:m){
      # Coefficients
      par_list[[i]] = par.new[index_start[i]:index_end[i]]
      # Differences of given order
      D <- diff(par_list[[i]], differences = pen_order[i])
      # Weights
      w[[i]] = 1 / (D ^ 2 + epsilon ^ 2)
      # Selected indexes
      sel[[i]] = w[[i]] * D ^ 2
    }
    
    # Convergence criterion 
    crit_list <- sapply(1:m, function(i) max(abs(old_sel[[i]] - sel[[i]])) < tol)
    
    # Product to ensure the criterion is met for each covariate
    converge <- all(crit_list)
    
    # Convergence chack
    if (converge) break
    
    # Update the selection index
    old_sel <- sel
  }
  
  return(list(
    sel = sel,           # Indicator of knot selection
    w = w,               # List of weights
    par.new = par.new,   # Parameters at convergence
    converge = converge, # Boolean indicating if the algorithm converged
    eta = eta            # Linear predictor at convergence
  ))
}


### ----------------------------------------------------------------------------
### Main Functions 
### ----------------------------------------------------------------------------

## AKSSAM.0 Algorithm ------------------------------------------------------------

# Function    : AKSSAM.0
# Description : Performs automatic knot selection in GAMs via an 
#               alternating implementation of the adaptive ridge 
#               and Fellner-Schall algorithm
#
# Input:
#   - X           : Matrix with observed variables for each covariate (columns)
#   - y           : Vector of observed values of the objective variable
#   - family      : String; either 'gaussian', 'poisson', or 'binomial',
#                   indicating the distribution of the objective variable
#   - lambda.init : Vector of penalization parameters used as initialization
#   - ndx         : Vector of inner intervals for each covariate
#   - bdeg        : Vector containing the degree of the B-spline basis for each covariate
#   - maxiter1    : Maximum iterations for weight optimization
#   - maxiter2    : Maximum iterations for penalization optimization
#   - maxiter3    : Maximum iterations for IRLS
#   - tol1        : Absolute tolerance for convergence in weight optimization
#   - tol2        : Relative tolerance for penalization optimization
#   - tol3        : Relative tolerance for IRLS
#   - epsilon     : Epsilon term for the adaptive ridge procedure
#
# Output:
#   - A list with the following elements:
#       * lambda         : Vector of penalizations maximizing the restricted log-likelihood
#       * K_sel          : Matrix with selected knots for each covariate (columns) for optimal penalization
#       * New_Design_list: List of design matrices after knot selection for each covariate and intercept
#                          corresponding to the optimal penalization
#       * alpha.new      : Ordered vector (intercept first, then covariates) storing resulting parameter 
#                          estimates for the optimal penalization
#
AKSSAM.0 = function(X, y, family, lambda.init, ndx, bdeg, 
                    maxiter1, maxiter2, maxiter3, tol1, tol2, tol3, 
                    epsilon){
  
  ## Initalize terms 
  # Identify the number of covariates and instances
  m = dim(X)[2]
  n = length(y)
  
  # Check for adequate ndx
  ndx_coerced <- pmin(ndx, floor(0.8 * sapply(1:m, function(i) length(unique(X[, i])))))
  # Coerce those inadequate ones
  if (any(ndx != ndx_coerced)){
    warning('Too many initial knots: the number of knots for certain covariates has been coerced', call. = FALSE)
    ndx = ndx_coerced
  }
  
  # Compute the equally-spaced knot vectors
  K <- lapply(1:m, function(i) my.knots(X[, i], min(X[, i]), max(X[, i]), ndx[i], bdeg[i]))
  
  # List with design matrices
  Design_list = vector("list", m + 1)
  
  # Intercept term
  Design_list[[1]] = matrix(1, n, 1)
  
  # Design matrices for each covariate:
  Design_list[-1] <- lapply(1:m, function(i) my.bbase4(X[,i], K[[i]], bdeg[i]))
  
  # Obtain the sizes of each basis
  basis_length = sapply(Design_list[-1], ncol)
  
  # Indexes for each covariate's parameters
  index_start <- cumsum(c(2, head(basis_length, -1)))
  index_end <- index_start + basis_length - 1
  
  # Construct the added design matrix in order B^* = [1:B1:...:Bm]
  B.new <- do.call(cbind, Design_list)
  # Obtain B.new' * B.new
  cross.B.New = crossprod(B.new)
  
  # Initialize the old selected knots list 
  old_sel <- lapply(1:m, function(i) rep(0, ncol(Design_list[[i+1]]) - bdeg[[i]] - 1))
  
  # Initialize the selected knots list
  sel <- vector('list', m) 
  
  # Penalization orders for each covariate
  pen_order <- bdeg + 1
  
  # Initialize a list containing the parameters
  par_list = vector('list', m)
  
  # Initialize the weights list
  w <- lapply(1:m, function(i) rep(1, ncol(Design_list[[i+1]]) - bdeg[[i]] - 1))
  
  # Initialize the penalizations
  lambda = lambda.init  # Old penalizations
  lambda.new = lambda   # New penalizations: Updated within first iteration
  
  # Initialize eta as the seed for the linear predictor
  if (family != "gaussian"){
    eta = IRLS.init(X, y, family, bdeg, maxiter3, tol3)
  } else {
    eta = NULL
  }
  
  
  ## Main Loop
  # Max iterations for weight optimization
  for (l in 1:maxiter1){
    
    # Weighted B-splines fit
    PP = construct.penalizations(Design_list, pen_order, lambda, w)
    
    # IRLS in GLM case
    ll = Fit.BSplines.Penalized(B.new, y, PP, family, maxiter3, tol3, eta)
    par.old = ll$par
    
    # In GLM case, retain eta and Omega weight matrix for penalization update
    if (family != "gaussian"){
      eta = ll$eta
      Omega_vec = ll$Omega_vec
      basis_and_pen_inv = robust_solve(crossprod(B.new * Omega_vec, B.new) + PP)
    } else{
      basis_and_pen_inv = robust_solve(cross.B.New + PP)
    }
    
    # Max iterations for lambda optimization for fixed weights
    for (k in 1:maxiter2){
      
      if(family == "gaussian"){
        # Obtain RSS 
        predicted = B.new %*% par.old   # Compute X*beta
        RSS = sum((y - predicted)^2)    # RSS
        
        # Obtain tr((X'X + S_lambda + P^I)^{-1}X'X)
        trace1 = sum(basis_and_pen_inv * t(cross.B.New))
        # Estimation of sigma^2
        sigma2 = RSS / (n - trace1)
      }
      
      # New penalization terms
      for (j in 1:m){
        
        Sj = construct.penalizations.deriv(Design_list, pen_order, w, j)
        
        # tr([(S_lambda + p^I)^{-1} - (X'X + S_lambda + p^I)^{-1}]Sj) --- Gaussian
        # tr([(S_lambda + p^I)^{-1} - (X' Omega X + S_lambda + p^I)^{-1}]Sj) --- GLM
        trace4 = sum((pseudoinverse(PP) - basis_and_pen_inv) * t(Sj))
        
        # Denominator truncation
        denom = as.numeric(crossprod(par.old, Sj %*% par.old))
        if (abs(denom) < 1e-6){
          denom = 1e-6
        }
        
        # Updated j-th element of lambda vector
        if (family == "gaussian"){ # Gaussian scenario: Omega = I, phi = sigma^2
          lambda.new[j] =  sigma2 * (trace4 / denom) * lambda[j] 
        } else{                    # Binomial or Poisson scenarios: Omega by IRLS, phi = 1
          lambda.new[j] = (trace4 / denom) * lambda[j] 
        } 
        
        # Negative update check
        if (lambda.new[j] < 0) lambda.new[j] = 0.1
        
      }
      
      # IRLS in GLM case
      PP = construct.penalizations(Design_list, pen_order, lambda.new, w)
      ll = Fit.BSplines.Penalized(B.new, y, PP, family, maxiter3, tol3, eta)
      par.new = ll$par
      
      if (family != "gaussian"){
        eta = ll$eta
        Omega_vec = ll$Omega_vec
        basis_and_pen_inv = robust_solve(crossprod(B.new * Omega_vec, B.new) + PP)
      } else{
        basis_and_pen_inv = robust_solve(cross.B.New + PP)
      }
      
      # Inner loop convergence: Relative convergence of the linear predictor
      numerator <- sum((B.new %*% (par.old - par.new))^2)
      denominator <-  sum((B.new %*% par.new)^2)
      
      converge1 = if (denominator < 1e-8 || numerator == Inf) Inf else numerator / denominator
      
      # Update the penalization and parameters
      lambda = lambda.new
      par.old = par.new
      
      # Convergence check
      if (converge1 < tol1) break
    }
    
    # Perform adaptive ridge once (weight update) after lambda optimization
    # Assign the coefficients and update the weights and selected vector
    for (i in 1:m){
      # Coefficients
      par_list[[i]] = par.new[index_start[i]:index_end[i]]
      # Differences of given order
      D <- diff(par_list[[i]], differences = pen_order[i])
      # Weights
      w[[i]] = 1 / (D ^ 2 + epsilon ^ 2)
      # Selected indexes
      sel[[i]] = w[[i]] * D ^ 2
    }
    
    # Outer loop convergence: Absolute convergence of selected vectors
    sel1 <- unlist(sel)
    old_sel1 <- unlist(old_sel)
    
    # Compute the relative squared error
    numerator <- sum((sel1 - old_sel1)^2)
    denominator <- sum(sel1^2)
    
    # Compute convergence score (numeric value)
    converge2 <- if (denominator < 1e-8 || numerator == Inf) Inf else numerator / denominator
    
    # Stop if convergence ratio is below threshold
    if (converge2 < tol2) break
    
    
    # Update the selection index
    old_sel <- sel
  }
  
  if (l == maxiter1 || k == maxiter2){
    # Iteration tolerance warning
    warning('Maximum number of iterations reached')
  }
  
  ## After-convergence knot selection
  # Obtain the selected knots in each covariate
  K_sel = lapply(1:m, function(i) {
    # Initial knots
    knots = K[[i]]
    # Identify the inner, selected knots (taking into account indexing, there are bdeg[i] + 1 outer knots)
    selected_index = which(sel[[i]] > 0.99) + bdeg[i] + 1
    # Get the external knots
    extra_index = c(1:(bdeg[i] + 1), (length(knots) - bdeg[i]):length(knots))
    # Join the selected and extra knots
    return(knots[sort(unique(c(selected_index, extra_index)))])
  })
  
  # Obtain the new design matrices
  New_Design_list = Design_list
  New_Design_list[-1] <- lapply(1:m, function(i) my.bbase4(X[,i], K_sel[[i]], bdeg[i]))
  
  # Obtain the final design matrix B^* = [1:B1:...:Bm]
  B.new <- do.call(cbind, New_Design_list)
  
  # Obtain the identifiability penalization
  PP = construct.penalizations2(New_Design_list)
  
  # Solve the B-Splines regression (with identifiability correction)
  ll =  Fit.BSplines.Penalized(B.new, y, PP, family, maxiter3, tol3, eta)
  alpha.new = ll$par
  
  ## Output 
  return(list(
    K_sel = K_sel,                       # Selected knots
    New_Design_list = New_Design_list,   # Design matrices
    alpha.new = alpha.new,               # Coefficients
    lambda_sel = lambda.new              # Selected penalization
  ))
}



## AKSSAM.0_fallback Algorithm ------------------------------------------------------------

# Function    : AKSSAM.0_fallback
# Description : Performs automatic knot selection in GAMs via an 
#               alternating implementation of the adaptive ridge 
#               and Fellner-Schall algorithm with an upper bound
#               in the updated parameters to avoid numerical instabilities.
#
# Input:
#   - X           : Matrix with observed variables for each covariate (columns)
#   - y           : Vector of observed values of the objective variable
#   - family      : String; either 'gaussian', 'poisson', or 'binomial',
#                   indicating the distribution of the objective variable
#   - lambda.init : Vector of penalization parameters used as initialization
#   - ndx         : Vector of inner intervals for each covariate
#   - bdeg        : Vector containing the degree of the B-spline basis for each covariate
#   - maxiter1    : Maximum iterations for weight optimization
#   - maxiter2    : Maximum iterations for penalization optimization
#   - maxiter3    : Maximum iterations for IRLS
#   - tol1        : Absolute tolerance for convergence in weight optimization
#   - tol2        : Relative tolerance for penalization optimization
#   - tol3        : Relative tolerance for IRLS
#   - epsilon     : Epsilon term for the adaptive ridge procedure
#
# Output:
#   - A list with the following elements:
#       * lambda         : Vector of penalizations maximizing the restricted log-likelihood
#       * K_sel          : Matrix with selected knots for each covariate (columns) for optimal penalization
#       * New_Design_list: List of design matrices after knot selection for each covariate and intercept
#                          corresponding to the optimal penalization
#       * alpha.new      : Ordered vector (intercept first, then covariates) storing resulting parameter 
#                          estimates for the optimal penalization
#

# If the original algorithm does not converge, we set an upper bound in lambda
AKSSAM.0_fallback = function(X, y, family, lambda.init, ndx, bdeg, 
                             maxiter1, maxiter2, maxiter3, tol1, tol2, tol3, 
                             epsilon){
  
  ## Initalize terms 
  # Identify the number of covariates and instances
  m = dim(X)[2]
  n = length(y)
  
  # Check for adequate ndx
  ndx_coerced <- pmin(ndx, floor(0.8 * sapply(1:m, function(i) length(unique(X[, i])))))
  # Coerce those inadequate ones
  if (any(ndx != ndx_coerced)){
    warning('Too many initial knots: the number of knots for certain covariates has been coerced', call. = FALSE)
    ndx = ndx_coerced
  }
  
  # Compute the equally-spaced knot vectors
  K <- lapply(1:m, function(i) my.knots(X[, i], min(X[, i]), max(X[, i]), ndx[i], bdeg[i]))
  
  # List with design matrices
  Design_list = vector("list", m + 1)
  
  # Intercept term
  Design_list[[1]] = matrix(1, n, 1)
  
  # Design matrices for each covariate:
  Design_list[-1] <- lapply(1:m, function(i) my.bbase4(X[,i], K[[i]], bdeg[i]))
  
  # Obtain the sizes of each basis
  basis_length = sapply(Design_list[-1], ncol)
  
  # Indexes for each covariate's parameters
  index_start <- cumsum(c(2, head(basis_length, -1)))
  index_end <- index_start + basis_length - 1
  
  # Construct the added design matrix in order B^* = [1:B1:...:Bm]
  B.new <- do.call(cbind, Design_list)
  # Obtain B.new' * B.new
  cross.B.New = crossprod(B.new)
  
  # Initialize the old selected knots list 
  old_sel <- lapply(1:m, function(i) rep(0, ncol(Design_list[[i+1]]) - bdeg[[i]] - 1))
  
  # Initialize the selected knots list
  sel <- vector('list', m) 
  
  # Penalization orders for each covariate
  pen_order <- bdeg + 1
  
  # Initialize a list containing the parameters
  par_list = vector('list', m)
  
  # Initialize the weights list
  w <- lapply(1:m, function(i) rep(1, ncol(Design_list[[i+1]]) - bdeg[[i]] - 1))
  
  # Initialize the penalizations
  lambda = lambda.init  # Old penalizations
  lambda.new = lambda   # New penalizations: Updated within first iteration
  
  # Initialize eta as the seed for the linear predictor
  if (family != "gaussian"){
    eta = IRLS.init(X, y, family, bdeg, maxiter3, tol3)
  } else {
    eta = NULL
  }
  
  
  ## Main Loop
  # Max iterations for weight optimization
  for (l in 1:maxiter1){
    
    # Weighted B-splines fit
    PP = construct.penalizations(Design_list, pen_order, lambda, w)
    
    # IRLS in GLM case
    ll = Fit.BSplines.Penalized(B.new, y, PP, family, maxiter3, tol3, eta)
    par.old = ll$par
    
    # In GLM case, retain eta and Omega weight matrix for penalization update
    if (family != "gaussian"){
      eta = ll$eta
      Omega_vec = ll$Omega_vec
      basis_and_pen_inv = robust_solve(crossprod(B.new * Omega_vec, B.new) + PP)
    } else{
      basis_and_pen_inv = robust_solve(cross.B.New + PP)
    }
    
    # Max iterations for lambda optimization for fixed weights
    for (k in 1:maxiter2){
      
      if(family == "gaussian"){
        # Obtain RSS 
        predicted = B.new %*% par.old   # Compute X*beta
        RSS = sum((y - predicted)^2)    # RSS
        
        # Obtain tr((X'X + S_lambda + P^I)^{-1}X'X)
        trace1 = sum(basis_and_pen_inv * t(cross.B.New))
        # Estimation of sigma^2
        sigma2 = RSS / (n - trace1)
      }
      
      # New penalization terms
      for (j in 1:m){
        
        Sj = construct.penalizations.deriv(Design_list, pen_order, w, j)
        
        # tr([(S_lambda + p^I)^{-1} - (X'X + S_lambda + p^I)^{-1}]Sj) --- Gaussian
        # tr([(S_lambda + p^I)^{-1} - (X' Omega X + S_lambda + p^I)^{-1}]Sj) --- GLM
        trace4 = sum((pseudoinverse(PP) - basis_and_pen_inv) * t(Sj))
        
        # Denominator truncation
        denom = as.numeric(crossprod(par.old, Sj %*% par.old))
        if (abs(denom) < 1e-6){
          denom = 1e-6
        }
        
        # Updated j-th element of lambda vector
        if (family == "gaussian"){ # Gaussian scenario: Omega = I, phi = sigma^2
          lambda.new[j] =  sigma2 * (trace4 / denom) * lambda[j] 
        } else{                    # Binomial or Poisson scenarios: Omega by IRLS, phi = 1
          lambda.new[j] = (trace4 / denom) * lambda[j] 
        } 
        
        # Negative update check
        if (lambda.new[j] < 0) lambda.new[j] = 0.1
        
        # Exaggerate update check if non-normal (Uncomment if needed)
        if (family != "gaussian"){
          if (lambda.new[j] > 1e+4) lambda.new[j] = 1e+4
        }
        
      }
      
      # IRLS in GLM case
      PP = construct.penalizations(Design_list, pen_order, lambda.new, w)
      ll = Fit.BSplines.Penalized(B.new, y, PP, family, maxiter3, tol3, eta)
      par.new = ll$par
      
      if (family != "gaussian"){
        eta = ll$eta
        Omega_vec = ll$Omega_vec
        basis_and_pen_inv = robust_solve(crossprod(B.new * Omega_vec, B.new) + PP)
      } else{
        basis_and_pen_inv = robust_solve(cross.B.New + PP)
      }
      
      # Inner loop convergence: Relative convergence of the linear predictor
      numerator <- sum((B.new %*% (par.old - par.new))^2)
      denominator <-  sum((B.new %*% par.new)^2)
      
      converge1 = if (denominator < 1e-8 || numerator == Inf) Inf else numerator / denominator
      
      # Update the penalization and parameters
      lambda = lambda.new
      par.old = par.new
      
      # Convergence check
      if (converge1 < tol1) break
    }
    
    # Perform adaptive ridge once (weight update) after lambda optimization
    # Assign the coefficients and update the weights and selected vector
    for (i in 1:m){
      # Coefficients
      par_list[[i]] = par.new[index_start[i]:index_end[i]]
      # Differences of given order
      D <- diff(par_list[[i]], differences = pen_order[i])
      # Weights
      w[[i]] = 1 / (D ^ 2 + epsilon ^ 2)
      # Selected indexes
      sel[[i]] = w[[i]] * D ^ 2
    }
    
    # Outer loop convergence: Absolute convergence of selected vectors
    sel1 <- unlist(sel)
    old_sel1 <- unlist(old_sel)
    
    # Compute the relative squared error
    numerator <- sum((sel1 - old_sel1)^2)
    denominator <- sum(sel1^2)
    
    # Compute convergence score (numeric value)
    converge2 <- if (denominator < 1e-8 || numerator == Inf) Inf else numerator / denominator
    
    # Stop if convergence ratio is below threshold
    if (converge2 < tol2) break
    
    
    # Update the selection index
    old_sel <- sel
  }
  
  if (l == maxiter1 || k == maxiter2){
    # Iteration tolerance warning
    warning('Maximum number of iterations reached')
  }
  
  ## After-convergence knot selection
  # Obtain the selected knots in each covariate
  K_sel = lapply(1:m, function(i) {
    # Initial knots
    knots = K[[i]]
    # Identify the inner, selected knots (taking into account indexing, there are bdeg[i] + 1 outer knots)
    selected_index = which(sel[[i]] > 0.99) + bdeg[i] + 1
    # Get the external knots
    extra_index = c(1:(bdeg[i] + 1), (length(knots) - bdeg[i]):length(knots))
    # Join the selected and extra knots
    return(knots[sort(unique(c(selected_index, extra_index)))])
  })
  
  # Obtain the new design matrices
  New_Design_list = Design_list
  New_Design_list[-1] <- lapply(1:m, function(i) my.bbase4(X[,i], K_sel[[i]], bdeg[i]))
  
  # Obtain the final design matrix B^* = [1:B1:...:Bm]
  B.new <- do.call(cbind, New_Design_list)
  
  # Obtain the identifiability penalization
  PP = construct.penalizations2(New_Design_list)
  
  # Solve the B-Splines regression (with identifiability correction)
  ll =  Fit.BSplines.Penalized(B.new, y, PP, family, maxiter3, tol3, eta)
  alpha.new = ll$par
  
  ## Output 
  return(list(
    K_sel = K_sel,                       # Selected knots
    New_Design_list = New_Design_list,   # Design matrices
    alpha.new = alpha.new,               # Coefficients
    lambda_sel = lambda.new              # Selected penalization
  ))
}



## AKSSAM Algorithm ------------------------------------------------------------

# Function    : AKSSAM
# Description : Performs automatic knot selection in GAMs via an 
#               alternating implementation of the adaptive ridge 
#               and Fellner-Schall algorithm and includes an upper
#               bound on the updated parameters if numerical instabilities
#               emerge
#
# Input:
#   - X           : Matrix with observed variables for each covariate (columns)
#   - y           : Vector of observed values of the objective variable
#   - family      : String; either 'gaussian', 'poisson', or 'binomial',
#                   indicating the distribution of the objective variable
#   - lambda.init : Vector of penalization parameters used as initialization
#   - ndx         : Vector of inner intervals for each covariate
#   - bdeg        : Vector containing the degree of the B-spline basis for each covariate
#   - maxiter1    : Maximum iterations for weight optimization
#   - maxiter2    : Maximum iterations for penalization optimization
#   - maxiter3    : Maximum iterations for IRLS
#   - tol1        : Relative tolerance for convergence in weight optimization
#   - tol2        : Relative tolerance for penalization optimization
#   - tol3        : Relative tolerance for IRLS
#   - epsilon     : Epsilon term for the adaptive ridge procedure
#
# Output:
#   - A list with the following elements:
#       * lambda         : Vector of penalizations maximizing the restricted log-likelihood
#       * K_sel          : Matrix with selected knots for each covariate (columns) for optimal penalization
#       * New_Design_list: List of design matrices after knot selection for each covariate and intercept
#                          corresponding to the optimal penalization
#       * alpha.new      : Ordered vector (intercept first, then covariates) storing resulting parameter 
#                          estimates for the optimal penalization
#
AKSSAM= function(X, y, family, lambda.init, ndx, bdeg, 
                 maxiter1, maxiter2, maxiter3, tol1, tol2, tol3, 
                 epsilon){
  
  result <- tryCatch({
    withCallingHandlers({
      # Attempt to perform non-robustified AKSSAM
      AKSSAM.0(X, y, family, lambda.init, ndx, bdeg,
                   maxiter1, maxiter2, maxiter3, tol1, tol2, tol3,
                   epsilon)
    }, warning = function(w) {
      warning("Robustification measures were needed", call. = FALSE)
      invokeRestart("muffleWarning")
    })
    
  }, error = function(e) {
    tryCatch({
      withCallingHandlers({
        # Attempt to perform robustified AKSSAM after numerical instabilities
        AKSSAM.0_fallback(X, y, family, lambda.init, ndx, bdeg,
                              maxiter1, maxiter2, maxiter3, tol1, tol2, tol3,
                              epsilon)
      }, warning = function(w) {
        warning("Robustification measures were needed", call. = FALSE)
        invokeRestart("muffleWarning")
      })
      
    }, error = function(e2) {
      # Raise an error
      stop("AKSSAM did not converge due to numerical instabilities.", call. = FALSE)
    })
  })
  
  return(list(
    K_sel = result$K_sel,                       # Selected knots
    New_Design_list = result$New_Design_list,   # Design matrices
    alpha.new = result$alpha.new,               # Coefficients
    lambda_sel = result$lambda_sel              # Selected penalization
  ))
}


### ----------------------------------------------------------------------------
### Extra Functions 
### ----------------------------------------------------------------------------

## GAM A-Splines ---------------------------------------------------------------

# Function    : GAM.asplines
# Description : Performs automatic knot selection in GAMs via an 
#               extension of the A-Splines algorithm
#
# Input:
#   - X           : Matrix with observed variables for each covariate (columns)
#   - y           : Vector of observed values of the objective variable
#   - K           : Matrix containing initial knots for each covariate (columns)
#   - lambda      : Vector containing the penalizations for each covariate
#   - bdeg        : Vector containing the degree of the B-spline basis for each covariate
#   - order_diffs : Vector containing the difference order for each covariate
#   - family      : String; either 'gaussian', 'poisson', or 'binomial',
#                    indicating the distribution of the objective variable
#   - maxiter     : Maximum iterations for the IRLS algorithm
#   - tol         : Relative tolerance for convergence of the IRLS algorithm
#   - epsilon     : Epsilon term for the adaptive ridge procedure
#   - eta_init    : Initialization of the linear predictor
#
# Output:
#   - A list with the following elements:
#       * K_sel          : Matrix with selected knots for each covariate (columns)
#       * New_Design_list: List of design matrices after knot selection for each covariate and intercept
#       * alpha.new      : Ordered vector (intercept first, then covariates) of resulting parameter estimates
#       * eta            : Estimated linear predictor
#
GAM.asplines = function(X, y, ndx, lambda, bdeg, family, maxiter, 
                        tol, epsilon, eta_init = NULL){
  
  ## Initalize terms 
  # Identify the number of covariates and instances
  m = dim(X)[2]
  n = length(y)
  
  # Check for adequate ndx
  ndx_coerced <- pmin(ndx, floor(0.8 * sapply(1:m, function(i) length(unique(X[, i])))))
  # Coerce those inadequate ones
  if (any(ndx != ndx_coerced)){
    warning('Too many initial knots: the number of knots for certain covariates has been coerced')
    ndx = ndx_coerced
  }
  
  # Compute the equally-spaced knot vectors
  K <- lapply(1:m, function(i) my.knots(X[, i], min(X[, i]), max(X[, i]), ndx[i], bdeg[i]))
  
  # List with design matrices
  Design_list = vector("list", m + 1)
  
  # Intercept term
  Design_list[[1]] = matrix(1,n,1)
  
  # Design matrices for each covariate:
  Design_list[-1] <- lapply(1:m, function(i) my.bbase4(X[,i], K[[i]], bdeg[i]))
  
  # Obtain the sizes of each basis
  basis_length <- sapply(Design_list[-1], ncol)
  
  # Construct the added design matrix in order B^* = [1:B1:...:Bm]
  B.new <- do.call(cbind, Design_list)
  
  # Initialize the old selected knots list 
  old_sel <- lapply(1:m, function(i) rep(0, ncol(Design_list[[i+1]]) - bdeg[[i]] - 1))
  
  # Initialize the weights list
  w <- lapply(1:m, function(i) rep(1, ncol(Design_list[[i+1]]) - bdeg[[i]] - 1))
  
  # Initialize eta 
  if (is.null(eta_init) && family != "gaussian"){
    eta = IRLS.init(X = X, y = y, family = family, maxiter = 50, tol = 1e-5)
  } else{
    eta = eta_init
  }
  
  ## Main Loop 
  # Adaptive Ridge
  ll = adridge(Design_list, basis_length, family, lambda, w, old_sel,
               B.new, y, bdeg, epsilon, maxiter, tol, eta)
  
  sel = ll$sel                           # Selected inner knot indexes
  w = ll$w                               # Weights after convergence
  par.new = ll$par.new                   # Vector of parameters
  if (family != 'gaussian') eta = ll$eta # Linear predictor
  
  # Boolean indicating convergence
  converge = ll$converge 
  
  
  # Final assignment of the algorithm
  if (converge){
    # Obtain the selected knots in each covariate
    K_sel = lapply(1:m, function(i) {
      # Initial knots
      knots = K[[i]]
      # Identify the inner, selected knots 
      selected_index = which(sel[[i]] > 0.99) + bdeg[[i]] + 1
      # Get the external knots
      extra_index = c(1:(bdeg[i] + 1), (length(knots) - bdeg[i]):length(knots))
      # Join the selected and extra knots
      return(knots[sort(unique(c(selected_index, extra_index)))])
    }) 
    
    # Obtain the new design matrices
    New_Design_list = Design_list
    New_Design_list[-1] <- lapply(1:m, function(i) my.bbase4(X[,i], K_sel[[i]], bdeg[i]))
    # Obtain the final design matrix B^* = [1:B1:...:Bm]
    B.new <- do.call(cbind, New_Design_list)
    
    # Obtain the identifiability penalization
    PP = construct.penalizations2(New_Design_list)
    
    # Solve the B-Splines regression (with identifiability correction)
    ll = Fit.BSplines.Penalized(B.new, y, PP, family, maxiter, tol, eta)
    alpha.new = ll$par
    if (family != 'gaussian') eta = ll$eta
  } else{
    # Error output
    warning('WPSS Algorithm didnt converge')
    return(NULL)
  }
  
  
  ## Output: List with the selected knots, design matrices and coefficients
  return(list(
    K_sel = K_sel,                       # Selected knots
    New_Design_list = New_Design_list,   # Design matrices
    alpha.new = alpha.new,               # Coefficients
    eta = eta                            # Linear predictor
  ))
}

