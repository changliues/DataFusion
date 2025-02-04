library(pracma)


minimax_solve <- function(gram1, gram2, g1, g2, lambda1, lambda2){
  
  """
  This function returns the coefficients of the solution function h(r1),
  which, along with f(r2), solves the following minimax optimization problem:

  h = \argmin_{h} \max_f [
      \mathbb{E} \big( f(r2) \cdot (h(r1) \cdot g1 + g2) - f(r2)^2 \big)
      - \lambda_2 \|f\|_{RKHS}^2 + \lambda_1 \|h\|_{RKHS}^2
  ]

  Here, h belongs to the RKHS defined by the kernel function associated with 
  the Gram matrix `gram1`, while f belongs to the RKHS defined by the kernel 
  function associated with the Gram matrix `gram2`.

  Returns:
    list: A list `[alpha, beta]`, where:
        - `alpha` is the coefficient vector of h in its RKHS.
        - `beta` is the coefficient vector of f in its RKHS.
  """

  n <- nrow(gram1)

  gamma_matrix <- 1/4 * gram2 %*% solve(1/n * gram2 + lambda2 * diag(1, nrow = n) , tol = 1e-30)
  minimax_1 <- gram1 %*% diag(g1) %*% gamma_matrix %*% diag(g1) %*% gram1 + n^2*lambda1*gram1

  alpha <- as.vector(-pinv(minimax_1, tol = 1e-30) %*% gram1 %*% diag(g1) %*% gamma_matrix %*% g2)
  beta <- as.vector(1/2*solve(1/n * gram2 + lambda2 * diag(1, nrow = n), tol = 1e-30) %*% t(1/n*(alpha %*% gram1 * g1 + g2)))

  return(list(alpha=alpha, beta=beta))
}

score_nuisance_function <- function(h_values, gram2_score, g1, g2, lambda_score) {

  """
  Scores a fitted nuisance function h by solving the minimax problem with respect 
  to the coefficient vector `beta` of function f, while keeping h(r1) fixed.  

  The function f belongs to the RKHS associated with the Gram matrix `gram_score`, 
  and the regularization parameter `lambda_2` is specified as the function input `lambda_score`.  

  Returns:
    float: The computed score, where a higher value indicates a better fit.
  """

  n <- nrow(gram2_score)
  
  beta <- as.vector(1/2*solve(1/n * gram2_score + lambda_score * diag(1, nrow = n), tol = 1e-30) %*% t(1/n*(h_values * g1 + g2)))
  f_values <- beta %*% gram2_score
  
  metric <- mean((g1 * h_values + g2) * f_values - f_values^2)
  return(-metric)
}
