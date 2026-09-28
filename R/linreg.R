#' Linear regression using QR decomposition
#'
#' A reference class that fits a multiple linear regression model. The
#' coefficients and their variances are computed with the QR decomposition
#' of the design matrix.
#'
#' @field formula The model formula
#' @field data The data frame used to fit the model.
#' @field data_name The name of the data object, as written in the call
#' @field beta Named vector of regression coefficients
#' @field y_hat Fitted values.
#' @field e Residuals.
#' @field df Degrees of freedom.
#' @field sigma2 Residual variance.
#' @field var_beta Variance matrix of the coefficients
#' @field t_values The t-values of the coefficients.
#' @field p_values The p-values of the coefficients.
#'
#' @importFrom methods setRefClass
#' @importFrom stats model.matrix pt
#'
#' @examples
#' mod <- linreg$new(Petal.Length ~ Species, data = iris)
#' mod$beta
#'
#' @export linreg
linreg <- setRefClass("linreg",
                      fields = list(
                        formula = "formula",
                        data = "data.frame",
                        data_name = "character",
                        X = "matrix",
                        y = "numeric",
                        beta = "numeric",
                        y_hat = "numeric",
                        e = "numeric",
                        df = "numeric",
                        sigma2 = "numeric",
                        var_beta = "matrix",
                        t_values = "numeric",
                        p_values = "numeric"
                      ),
methods = list(
  initialize = function(formula, data, ...){
    callSuper(...)
    #empty calls, e.g $copy() should not count anything
    if(!missing(formula) && !missing(data)){
      X_mat           <- model.matrix(formula, data)
      y_vec           <- data[[all.vars(formula)[1]]]

      QR              <- qr(X_mat)
      Q               <- qr.Q(QR)
      R               <- qr.R(QR)

      b               <- as.vector(backsolve(R, t(Q) %*% y_vec))
      names(b)        <- colnames(X_mat)

      yh              <- as.vector(X_mat %*% b)
      res             <- y_vec - yh
      dof             <- nrow(X_mat) - ncol(X_mat)
      s2              <- sum(res^2) / dof

      R_inv           <- backsolve(R, diag(ncol(R)))
      vb              <- s2 * R_inv %*% t(R_inv)

      tv              <- b / sqrt(diag(vb))
      pv              <- 2 * pt(-abs(tv), dof)

      .self$formula   <- formula
      .self$data      <- data
      .self$data_name <- deparse(substitute(data))
      .self$X         <- X_mat
      .self$y         <- y_vec
      .self$beta      <- b
      .self$y_hat     <- yh
      .self$e         <- res
      .self$df        <- dof
      .self$sigma2    <- s2
      .self$var_beta  <- vb
      .self$t_values  <- tv
      .self$p_values  <- pv
    }
  }
  )
)
