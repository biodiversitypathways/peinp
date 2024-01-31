library(spOccupancy)

set.seed(500)
# Sites
J.x <- 10
J.y <- 10
J <- J.x * J.y
# Primary time periods
n.time <- sample(5:10, J, replace = TRUE)
n.time.max <- max(n.time)
# Replicates
n.rep <- matrix(NA, J, max(n.time))
for (j in 1:J) {
  n.rep[j, 1:n.time[j]] <- sample(1:4, n.time[j], replace = TRUE)
}
# Occurrence --------------------------
beta <- c(0.4, 0.5, -0.9)
trend <- TRUE
sp.only <- 0
psi.RE <- list()
# Detection ---------------------------
alpha <- c(-1, 0.7, -0.5)
p.RE <- list()
# Temporal parameters -----------------
rho <- 0.7
sigma.sq.t <- 0.6

# Get all the data
dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = n.time, n.rep = n.rep,
               beta = beta, alpha = alpha, sp.only = sp.only, trend = trend,
               psi.RE = psi.RE, p.RE = p.RE, sp = FALSE, ar1 = TRUE,
               sigma.sq.t = sigma.sq.t, rho = rho)

# Package all data into a list
# Occurrence
occ.covs <- list(int = dat$X[, , 1],
                 trend = dat$X[, , 2],
                 occ.cov.1 = dat$X[, , 3])
# Detection
det.covs <- list(det.cov.1 = dat$X.p[, , , 2],
                 det.cov.2 = dat$X.p[, , , 3])
# Data list bundle
data.list <- list(y = dat$y,
                  occ.covs = occ.covs,
                  det.covs = det.covs)
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   alpha.normal = list(mean = 0, var = 2.72),
                   rho.unif = c(-1, 1),
                   sigma.sq.t.ig = c(2, 0.5))

# Starting values
z.init <- apply(dat$y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(beta = 0, alpha = 0, z = z.init)

# Tuning
tuning.list <- list(rho = 0.5)

n.batch <- 20
batch.length <- 25
n.samples <- n.batch * batch.length
n.burn <- 100
n.thin <- 1

# Run the model
out <- tPGOcc(occ.formula = ~ trend + occ.cov.1,
              det.formula = ~ det.cov.1 + det.cov.2,
              data = data.list,
              inits = inits.list,
              priors = prior.list,
              tuning = tuning.list,
              n.batch = n.batch,
              batch.length = batch.length,
              verbose = TRUE,
              ar1 = TRUE,
              n.report = 25,
              n.burn = n.burn,
              n.thin = n.thin,
              n.chains = 1)

summary(out)
