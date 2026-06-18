#### Functions for GAM Asplines and Wood's algorithm ####

## Basic Functions -------------------------------------------------------------

# my.knots
# Function which calculates the adequate equally-spaced knots
# Input:  - x: The vector of covariate values
#         - xl: Left-side of the desired knot interval
#         - xr: Right-side of the desired knot interval
#         - ndx: Nº of inner intervals
#         - bdeg: Degree of the B-Spline basis (For extra knots)
# Output: - knots: Vector of adequate equally-spaced knots
# 
my.knots = function(x, xl = min(x), xr = max(x), ndx = 10, bdeg = 3){
  
  dx <- (xr - xl) / ndx
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
  
  knots
}


# my.bbase4
# Function which constructs the design B-splines matrix given the knot vector
# Input:  - x: The vector of covariate values
#         - knots: Vector of knots
#         - bdeg: Degree of the  B-Spline basis
# Output: - B: Design matrix of the B-Spline basis
# 
my.bbase4 = function (x, knots, bdeg = 3) 
{
  B <- splineDesign(knots, x, bdeg + 1, outer.ok = TRUE)
  # outer.ok = TRUE to avoid issues when the knot vector is large
  # (ndx = 60 onwards)
  
  B
}

## Penalizations ---------------------------------------------------------------

# construct.penalizations
# Function which creates the penalizations (A-splines + Identifiability)
# for the main loop
# Input:  - Design_list: List with the design matrices for each covariate
#         - order_diffs: Vector with the difference order for each covariate
#         - lambda: Vector with the penalization for each covariate
#         - w: List containing vectors of penalization weights for each 
#              covariate
# Output: - PP: Matrix constructed block-wise to account for the overall 
#               penalization
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
    B = Design_list[[i]]                                    # Design matrix
    D = diff(diag(ncol(B)), differences = order_diffs[i-1]) # Difference matrix
    P = lambda[i-1] * crossprod(D * sqrt(w[[i-1]]))         # Weighted Smoothing Penalization
    PI = tcrossprod(crossprod(B, ones))                     # Identifiability penalization
    return(P + PI)                                          # Join both penalizations
  })
  
  # Unify all the penalizations block-wise and return
  PP = as.matrix(bdiag(Pen_list))
}


# construct.penalizations2
# Function which just adds the identifiability penalty 
# (After-convergence B-splines problem)
# Input:  - Design_list: List with the design matrices for each covariate
# Output: - PP: Matrix constructed block-wise to account for the overall 
#               penalization
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
    B = Design_list[[i]]                            # Design matrix
    return(tcrossprod(crossprod(B, ones)))          # Identifiability penalization
  })
  
  # Unify all the penalizations block-wise and return
  PP = as.matrix(bdiag(Pen_list))
}


# contruct.penalizations.deriv
# Function which constructs the matrix Sj in Woods algorithm.
# 
# Input:  - Design_list: List with the design matrices for each covariate
#         - order_diffs: Vector with the difference order for each covariate
#         - w: List containing vectors of penalization weights for each 
#              covariate
#         - j : Integer indicating the covariate.
# Output: - PP: Matrix Sj, which is a block-wise constructed matrix which has
#               Dj Wj Dj in the j-th covariate's spot an 0's otherwise. 
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
    B = Design_list[[i]] # Obtain the design mat
    if (i == j){ # Only calculate block if i == j
      
      # Difference matrix
      D = diff(diag(ncol(B)), differences = order_diffs[i-1]) 
      
      # Weighted Smoothing Penalization without lambda_j
      t(D) %*% diag(w[[i-1]]) %*% D           
    } else{ # If i != j then make the block 0
      matrix(0, dim(t(B)%*%B)[1], dim(t(B)%*%B)[2])
    }
  })
  
  # Unify all the penalizations block-wise and return
  PP = as.matrix(bdiag(Pen_list))
}


## IRLS ------------------------------------------------------------------------

# We took the following stability measures:
#   
# - Initialize via an oversmoothed and unpenalized IRLS instance to obtain an initial feasible seed for the algorithm
# - Added extra +1e-6 term in potentially unstable denominators.

### Initialitation for IRLS

# IRLS.init
# Function which fits unpenalized B-splines in a GLM context recursively using 
# Newton-Rhapson for ndx = 5. To be used as starting seed for IRLS.
# Input:  - X: Matrix which stores, column-wise the observed variables for each
#              covariate 
#         - y: Vector of the objective variable's observed values
#         - family: String, either 'gaussian', 'poisson' or 'binomial' which 
#                   indicates the distribution of the objective variable.
#         - maxiter: Maximum iterations for the IRLS algorithm tolerance
#         - tol: Relative tolerance for the convergence of the IRLS algorithm
# Output: - eta: Linear predictor at convegence
#                
IRLS.init = function(X, y, family, maxiter = 50, tol = 1e-5){
  
  ## Correctly store the GLM family --------------------------------------------
  family <- match.arg(family, choices = c("poisson", "binomial"))
  
  # Distribution family match
  if (family == "poisson") {
    g_inv <- function(x) exp(x)
    dg <- function(x) 1 / (x + 1e-6)
    V <- function(x) x
  } else if (family == "binomial") {
    g_inv <- function(x) 1 / (1 + exp(-x))
    dg <- function(x) 1 / (x * (1 - x) + 1e-6)
    V <- function(x) x * (1 - x)
  }
  
  ## Construct the smaller basis -----------------------------------------------
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
  
  ## IRLS Algorithm with Max iteration tolerance -------------------------------
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
  
  ## Output --------------------------------------------------------------------
  return(eta_new) # Estimated linear predictor
}

### Main Algorithm

# Fit.BSplines.Penalized
# Function which fits penalized B-splines in a GLM context recursively using 
# Newton-Rhapson
# Input:  - B: Design matrix 
#         - y: Vector of the objective variable's observed values
#         - P: Penalization matrix of the overall problem
#         - family: String, either 'gaussian', 'poisson' or 'binomial' which 
#                   indicated the distribution of the objective variable.
#         - maxiter: Maximum iterations for the IRLS algorithm tolerance
#         - tol: Relative tolerance for the convergence of the IRLS algorithm
#         - eta_init: Initialization value for the linear predictor
# Output: - par: Ordered vector by intercept and then by covariate order which
#                stores the resulting parameter estimations of the IRLS 
#                algorithm
#         - Omega_vec: Vector of weights of the IRLS-associated weight matrix.
#                
Fit.BSplines.Penalized = function(B, y, P, family, maxiter = 50, tol = 1e-5, eta_init = NULL){
  
  # Correctly store the GLM family
  family <- match.arg(family, choices = c("gaussian", "poisson", "binomial"))
  
  # Distribution family match
  if (family == "poisson") {
    g_inv <- function(x) exp(x)
    dg <- function(x) 1 / (x + 1e-6)
    V <- function(x) x
  } else if (family == "binomial") {
    g_inv <- function(x) 1 / (1 + exp(-x))
    dg <- function(x) 1 / (x * (1 - x) + 1e-6)
    V <- function(x) x * (1 - x)
  } else if (family == "gaussian") {
    par <- qr.solve(crossprod(B) + P, crossprod(B, y), tol = 1e-18)
    return(list(par = par))
  }
  
  # GLM Scenario
  
  # Initialize the linear predictor
  eta_old <- eta_init
  
  # Penalized IRLS Algorithm with Max iteration tolerance
  
  # Iteration Tolerance
  for (iter in 1:maxiter){
    
    # Calculate mu and Omega
    mu <- g_inv(as.vector(eta_old))
    Omega_vec <- V(mu)
    
    # Solve the IRLS algorithm instance
    z <- eta_old +  (y - mu) * dg(mu)
    B_weighted <- B * Omega_vec 
    par <- qr.solve(crossprod(B_weighted, B) + P, crossprod(B_weighted, z), tol = 1e-18)
    
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
  
  list(
    par = par,             # Estimated parameters at IRLS convergence
    Omega_vec = Omega_vec, # Estimated weights vector
    eta = eta_new          # Estimated linear predictor
  )
}



## Adaptive Ridge Algorithm ----------------------------------------------------

# adridge
# Auxiliary function to ease the algorithm formulation.
# Performs the adaptive ridge algorithm. 
# 
# Input:  - Design_list: A list containing the design matrices for each covariate 
#                        and the intercept.
#         - basis_length: A vector specifying the length of the B-spline basis 
#                         for each covariate.        
#         - order_diffs: Vector containing the difference order for each 
#                        covariate.
#         - family: String, either 'gaussian', 'poisson' or 'binomial' which 
#                   indicated the distribution of the objective variable.
#         - lambda: Vector containing the penalizations for each covariate.
#         - w: A list of initial weights for each covariate.
#         - old_sel: A list of initial indicator of selected inner knots for 
#                    each covariate. 
#         - par: A list of initial parameter estimates for each covariate.
#         - B.new: The full B-spline basis matrix.
#         - y: Vector which stores the observed values of the objective variable.
#         - maxiter: Maximum iterations for the WPSS algorithm tolerance.
#         - tol: Relative tolerance for the convergence of the WPSS algorithm.
# Output: - List with the following:
#              - sel: A list containing the indicator of selected inner knot for 
#                     each covariate.
#              - w: A list of updated weights for each covariate.
#              - par.new: Vector which contains the parameter estiamtes for the 
#                         WPSS solution.
#              - converge: A boolean indicating whether the algorithm converged 
#                          within the specified tolerance.
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



## GAM A-Splines ---------------------------------------------------------------

# GAM.asplines3
# Function which  performs Automatic knot selection for B-splines
# Input:  - X: Matrix which stores, column-wise the observed variables for each
#              covariate
#         - y: Vector which stores the observed values of the objective variable
#         - K: Matrix containing column-wise the initial knots for each 
#              covariate 
#         - lambda: Vector containing the penalizations for each covariate
#         - bdeg: Vector containing the degree of the B-spline basis for each 
#                 covariate
#         - order_diffs: Vector containing the difference order for each 
#                        covariate 
#         - family: String, either 'gaussian', 'poisson' or 'binomial' which 
#                   indicated the distribution of the objective variable.
#         - maxiter: Maximum iterations for the IRLS algorithm tolerance
#         - tol: Relative tolerance for the convergence of the IRLS algorithm
#         - epsilon: Epsilon term for the adaptive ridge procedure
#         - par.init: List which contains, by covariate, an initial approx
#                     to the WPSS solution. 
#                     If ommitted, it is automatically passed by as 1s.
# Output: - List with the following:
#              - K_sel: Matrix containing column-wise the selected knots for each 
#                       covariate
#              - New_Design_list: List containing the design matrices after knot 
#                                 selection for each covariate and intercept
#              - alpha.new: Ordered vector by intercept and then by covariate 
#                           order which stores the resulting parameter 
#                           estimations of the algorithm
#             - par.WPSS: List which contains, by covariate, the WPSS solution.
#
GAM.asplines3 = function(X, y, ndx, lambda, bdeg, family, maxiter, 
                         tol, epsilon, eta_init = NULL){
  ## Initalize terms -----------------------------------------------------------
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
  
  ## Main Loop -----------------------------------------------------------------
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
  
  
  ## Output --------------------------------------------------------------------
  # Output a list with the selected knots, design matrices and coefficients
  list(K_sel = K_sel,                       # Selected knots
       New_Design_list = New_Design_list,   # Design matrices
       alpha.new = alpha.new,               # Coefficients
       eta = eta                            # Linear predictor
  )
}




## Wood's Algorithm ------------------------------------------------------------

# Notice that we apply the following modifications to the original algorithm formulation: 
#   
#   - We apply the linearity of the trace operator: $\textrm{tr}(A) + \textrm{tr}(B) = \textrm{tr}(A + B)$. Therefore, in the update of the penalization term, we obtain $$
#   \text{tr}\left(\left(S_{\lambda} + p^I\right)^{-} S_j \right) - \text{tr}\left(\left(X'X + S_{\lambda} + p^I\right)^{-1} S_j\right) = \text{tr}\left( \left[\left(S_{\lambda} + p^I\right)^{-} - \left(X'X + S_{\lambda} + p^I\right)^{-1}\right] S_j \right).$$
#   
#   - In order to avoid numerical errors via the quotient  $$\lambda_j^* = \sigma^2 \frac{\text{tr}((S_{\lambda}+P^I)^{-} S_j) - \text{tr}((\mathbf{X}^\top \mathbf{X} + S_{\lambda}+P^I)^{-1} S_j)}{\hat{\boldsymbol{\beta}}_{\lambda}^\top S_j \hat{\boldsymbol{\beta}}_{\lambda}} \lambda_j,$$ we will truncate the term $\hat{\boldsymbol{\beta}}_{\lambda}^\top S_j \hat{\boldsymbol{\beta}}_{\lambda}$ by `1e-6`.
# 
# - Notice that in the GLM case, the penalization update is $$\lambda_j^* = \phi \frac{\text{tr}((S_{\lambda}+P^I)^{-} S_j) - \text{tr}((\mathbf{X}^\top \Omega \mathbf{X} + S_{\lambda}+P^I)^{-1} S_j)}{\hat{\boldsymbol{\beta}}_{\lambda}^\top S_j \hat{\boldsymbol{\beta}}_{\lambda}} \lambda_j, $$ where $\Omega$ is the diagonal weight matrix at convergence of the IRLS algorithm.
# 
# - We'll add as well `qr.solve` to invert the required matrices to avoid ill-conditioning errors, since those seem to vanish after a few iterations.
# 
# - We added a modification such that, if the updated lambda is negative (numerical issues), it is set back to 0.1.
# 
# - We made an extra initialization step where we fit an unpenalized and oversmoothed IRLS instance in the GLM part for which we obtain a feasible initial value for the linear predictor. 


# GAM.asplines3.wood2
# Function which performs iterative maximization of the log-restricted marginal 
# likelihood of a Gaussian GAM model via the penalization in the aridge algorithm.
# After convergence of this algorithm, the function performs aspline regression
# with the selected parameter
# 
# Input:  - X: Matrix which stores, column-wise the observed variables for each
#              covariate
#         - y: Vector which stores the observed values of the objective variable
#         - family: String, either 'gaussian', 'poisson' or 'binomial' which 
#                   indicated the distribution of the objective variable.
#         - lambda.init: Vector of penalization parameters to use as initialization
#         - ndx: Vector of inner intervals for each covariate 
#         - bdeg: Vector containing the degree of the B-spline basis for each 
#                 covariate
#         - maxiter1: Maximum iterations for weight optimization
#         - maxiter2: Maximum iterations for penalization optimization
#         - maxiter3: Maximum iterations for IRLS
#         - tol1: Absolute tolerance for the convergence for weight optimization
#         - tol2: Relative tolerance for penalization optimization
#         - tol3: Relative tolerance for IRLS
#         - epsilon: Epsilon term for the adaptive ridge procedure
# Output: - List with the following:
#              - lambda: Vector of penalizations rendering the maximum restricted 
#                        log-likelihood coefficients.
#              - K_sel: Matrix containing column-wise the selected knots for each 
#                       covariate for the optimal penalization.
#              - New_Design_list: List containing the design matrices after knot 
#                                 selection for each covariate and intercept
#                                 for the optimal penalization.
#              - alpha.new: Ordered vector by intercept and then by covariate 
#                           order which stores the resulting parameter 
#                           estimations of the algorithm for the optimal 
#                           penalization.
#
GAM.asplines3.wood2 = function(X, y, family, lambda.init, ndx, bdeg, 
                               maxiter1, maxiter2, maxiter3, tol1, tol2, tol3, 
                               epsilon){
  ## Initalize terms -----------------------------------------------------------
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
    eta = IRLS.init(X, y, family, maxiter3, tol3)
  } else {
    eta = NULL
  }
  
  
  ## Main Loop -----------------------------------------------------------------
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
      basis_and_pen_inv = qr.solve(crossprod(B.new * Omega_vec, B.new) + PP, tol = 1e-18)
    } else{
      basis_and_pen_inv = qr.solve(cross.B.New + PP, tol = 1e-18)
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
        basis_and_pen_inv = qr.solve(crossprod(B.new * Omega_vec, B.new) + PP, tol = 1e-18)
      } else{
        basis_and_pen_inv = qr.solve(cross.B.New + PP, tol = 1e-18)
      }
      
      # Inner loop convergence: Relative convergence of the linear predictor
      converge1 = sum((B.new %*% (par.old - par.new))^2) / sum((B.new %*% par.new)^2)
      
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
    converge2 <- if (denominator < 1e-8) Inf else numerator / denominator
    
    # Stop if convergence ratio is below threshold
    if (converge2 < tol2) break
    
    
    # Update the selection index
    old_sel <- sel
  }
  
  if (l == maxiter1 || k == maxiter2){
    # Iteration tolerance warning
    warning('Maximum number of iterations reached')
  }
  
  ## After-convergence knot selection ------------------------------------------
  
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
  
  ## Output --------------------------------------------------------------------
  list(K_sel = K_sel,                       # Selected knots
       New_Design_list = New_Design_list,   # Design matrices
       alpha.new = alpha.new,               # Coefficients
       lambda_sel = lambda.new              # Selected penalization
  )
}