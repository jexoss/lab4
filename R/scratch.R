formula <- Petal.Length ~ Species
data <- iris

X <- model.matrix(formula, data)
y <- data[[all.vars(formula)[1]]]


QR <- qr(X)
Q <- qr.Q(QR)
R <- qr.R(QR)

beta <- backsolve(R, t(Q) %*% y)

beta <- as.vector(beta)
names(beta) <- colnames(X)
beta

y_hat <- as.vector(X %*% beta)
e <- y - y_hat
df <- nrow(X) - ncol(X)
sigma2 <- sum(e^2)/df


R_inv <- backsolve(R, diag(ncol(R)))
var_beta <- sigma2 * R_inv %*% t(R_inv)

mod <- lm(Petal.Length ~ Species, iris)
all.equal(sigma2, summary(mod)$sigma^2)
all.equal(var_beta, vcov(mod), check.attributes = FALSE)

se <- sqrt(diag(var_beta))
t_values <- beta / se
p_values <- 2 * pt(-abs(t_values), df)

coefs <- summary(mod)$coefficients
all.equal(unname(t_values), unname(coefs[, "t value"]))
all.equal(unname(p_values), unname(coefs[, "Pr(>|t|)"]))
