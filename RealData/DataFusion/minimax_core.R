library(pracma)


minimax_solve <- function(gram1, gram2, g1, g2, lambda1, lambda2){
  n <- nrow(gram1)

  gamma_matrix <- 1/4 * gram2 %*% solve(1/n * gram2 + lambda2 * diag(1, nrow = n) , tol = 1e-30)
  minimax_1 <- gram1 %*% diag(g1) %*% gamma_matrix %*% diag(g1) %*% gram1 + n^2*lambda1*gram1

  alpha <- as.vector(-pinv(minimax_1, tol = 1e-30) %*% gram1 %*% diag(g1) %*% gamma_matrix %*% g2)
  beta <- as.vector(1/2*solve(1/n * gram2 + lambda2 * diag(1, nrow = n), tol = 1e-30) %*% t(1/n*(alpha %*% gram1 * g1 + g2)))

  return(list(alpha=alpha, beta=beta))
}

score_nuisance_function <- function(h_values, gram2_score, g1, g2, lambda_score) {
  n <- nrow(gram2_score)
  
  beta <- as.vector(1/2*solve(1/n * gram2_score + lambda_score * diag(1, nrow = n), tol = 1e-30) %*% t(1/n*(h_values * g1 + g2)))
  f_values <- beta %*% gram2_score
  
  metric <- mean((g1 * h_values + g2) * f_values - f_values^2)
  return(-metric)
}