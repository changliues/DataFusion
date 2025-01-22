# Utility Functions
unvec <- function(a) {
  if (is.null(nrow(a))) {
    matrix(a, ncol = 1)
  } else {
    a
  }
}

join_data <- function(...) {
  args <- list(...)
  do.call(cbind, lapply(args, unvec))
}

# ProxCIData Class
ProxCIData <- setRefClass(
  "ProxCIData",
  fields = list(
    X = "matrix",
    Z = "numeric",
    M = "numeric",
    Y = "numeric",
    A = "numeric",
    G = "numeric",
    Xs = "matrix",
    n = "numeric",
    split_indices = "list"
  ),
  
  methods = list(
    initialize = function(..., X = matrix(), Z = numeric(), M = numeric(), Y = numeric(), A = numeric(), G = numeric(), Xs = matrix()) {
      .self$X <<- X
      .self$Z <<- Z
      .self$M <<- M
      .self$Y <<- Y
      .self$A <<- A
      .self$G <<- G
      .self$Xs <<- Xs
      .self$n <<- nrow(X)
      if(length(A)){
        stopifnot(length(A) == .self$n)
      }
      if(length(M)){
        stopifnot(length(M) == .self$n)
      }
      if(length(Y)){
        stopifnot(length(Y) == .self$n)
      }
      if(length(A)){
        stopifnot(length(A) == .self$n)
      }
      if(length(G)){
        stopifnot(length(G) == .self$n)
      }
    },
    
    r1 = function() {
      return(cbind(unvec(.self$M), unvec(.self$A), unvec(.self$X)))
    },
    
    r2 = function() {
      return(cbind(unvec(.self$Z), unvec(.self$A), unvec(.self$X)))
    },
    
    create_crossfit_split = function(n_splits, cache = TRUE) {
      if (n_splits <= 1) {
        idx <- seq_len(.self$n)
        split_indices <<- list(list(train = idx, eval = idx))
      } else {
        idx <- sample(seq_len(.self$n))
        split_size <- floor(.self$n / n_splits)
        split_indices <<- list()
        for (i in seq_len(n_splits)) {
          start <- (i - 1) * split_size + 1
          end <- i * split_size
          if (i == n_splits && end < .self$n) {
            end <- .self$n
          }
          eval_idx <- sort(idx[start:end])
          train_idx <- sort(idx[-(start:end)])
          stopifnot(length(eval_idx) + length(train_idx) == .self$n)
          split_indices[[i]] <<- list(train = train_idx, eval = eval_idx)
        }
      }
      
      if (cache) {
        .self$split_indices <- split_indices
      }
      
      return(split_indices)
    }
  )
)

