
pei_dyn_occ2 <- pei_occu |>
  mutate(occupied = case_when(species_code == 'BTNW' ~ 1, TRUE ~ 0)) |>
  group_by(location, year) |>
  mutate(recs = n_distinct(recording_date_time)) |>
  ungroup() |>
  filter(recs == 4) |>
  group_by(location, year) |>
  crossing(sample = 1:4) |>
  ungroup()

pei_dyn_occ <- pei_occu |>
  inner_join(distance_to_coast) |>
  mutate(julian = yday(recording_date_time),
         presence = case_when(species_code == 'BTNW' ~ 1, TRUE ~ 0)) |>
  select(location, year, recording_date_time, julian, hour, species_code, habitat_nesting, coast_dist, presence) |>
  distinct() |>
  group_by(location) |>
  mutate(survey = dense_rank(recording_date_time)) |>
  ungroup() |>
  mutate(survey_name = paste0("det",survey),
         survey_date = str_replace(survey_name,"det","date")) |>
  filter(species_code == "BTNW") |>
  pivot_wider(names_from = survey_name, values_from = presence, values_fill = 0) |>
  pivot_wider(names_from = survey_date, values_from = julian, values_fill = 0)


pei_y_cross <- as.matrix(pei_dyn_occ[,9:28])
pei_DATE <- as.matrix(pei_dyn_occ[,29:48])
pei_y_cross[is.na(pei_DATE) != is.na(pei_y_cross)] <- NA

sd_date <- sd(c(pei_DATE), na.rm = T)
mean_date <- mean(pei_DATE, na.rm = T)
pei_DATE <- (pei_DATE - mean_date) / sd_date

pei_years <- pei_dyn_occ |>
  select(year) |>
  distinct() |>
  mutate(year = as.character(year)) |>
  pull()

covs <- pei_dyn_occ[,7] |> as.data.frame()

umf_pei <- unmarkedMultFrame(y = pei_y_cross,
                             siteCovs = covs,
                             yearlySiteCovs = list(year = pei_years),
                             obsCovs = list(date=pei_DATE),
                             numPrimary = 5)


umf <- unmarkedMultFrame(
  y = pei_dyn_occ$presence,
  siteCovs = pei_dyn_occ[, c("coast_dist")],
  obsCovs = pei_dyn_occ[, c("hour")],
  numPrimary = 5
)

# Constant parameters
fm0 <- colext(~1, ~1, ~1, ~1, umf_pei)
# Year-dependent detection
fm1 <- colext(~1, ~1, ~1, ~year, umf_pei)
# Like fm0, but with year-dependent colonization and extinction
fm2 <- colext(~1, ~year-1, ~year-1, ~1, umf_pei)
# A fully time-dependent model
fm3 <- colext(~1, ~year-1, ~year-1, ~year, umf_pei)
# Like fm3 with ocean-dependence of 1st-year occupancy
fm4 <- colext(~coast_dist, ~year-1, ~year-1, ~year, umf_pei)




# Like fm4 with date- and year-dependence of detection
fm5 <- colext(~coast, ~year-1, ~year-1, ~year + date + I(date^2), umf_pei, starts=c(coef(fm4), 0, 0))
# Same as fm5, but with detection in addition depending on forest cover
fm6 <- colext(~forest, ~year-1, ~year-1, ~year + date + I(date^2) + forest, umf)
























data(crossbill)
colnames(crossbill)

DATE <- as.matrix(crossbill[,32:58])
y.cross <- as.matrix(crossbill[,5:31])
y.cross[is.na(DATE) != is.na(y.cross)] <- NA


sd.DATE <- sd(c(DATE), na.rm=TRUE)
mean.DATE <- mean(DATE, na.rm=TRUE)
DATE <- (DATE - mean.DATE) / sd.DATE
years <- as.character(1999:2007)
years <- matrix(years, nrow(crossbill), 9, byrow=TRUE)
umf <- unmarkedMultFrame(y=y.cross,
                         siteCovs=crossbill[,2:3], yearlySiteCovs=list(year=years),
                         obsCovs=list(date=DATE),
                         numPrimary=9)




fm0 <- colext(~1, ~1, ~1, ~1, umf)

# Like fm0, but with year-dependent detection
fm1 <- colext(~1, ~1, ~1, ~year, umf)

# Like fm0, but with year-dependent colonization and extinction
fm2 <- colext(~1, ~year-1, ~year-1, ~1, umf)

# A fully time-dependent model
fm3 <- colext(~1, ~year-1, ~year-1, ~year, umf)

# Like fm3 with forest-dependence of 1st-year occupancy
fm4 <- colext(~forest, ~year-1, ~year-1, ~year, umf)

# Like fm4 with date- and year-dependence of detection
fm5 <- colext(~forest, ~year-1, ~year-1, ~year + date + I(date^2),
              umf, starts=c(coef(fm4), 0, 0))

# Same as fm5, but with detection in addition depending on forest cover
fm6 <- colext(~forest, ~year-1, ~year-1, ~year + date + I(date^2) +
                forest, umf)






years <- as.character(2019:2023)
years <- matrix(years, nrow(pei_occu), 5, byrow = TRUE)




dyn_DATE <- pei_dyn_occ |>
  mutate(survey_number = row_number()) |>
  pivot_wider(names_from = survey_number, values_from = presence, values_fill = 0) |>
  filter(species_code == "BTNW")








M <- 30                            # Number of sites
J <- 4                                  # num secondary sample periods
T <- 10
psi <- rep(NA, T)                       # Occupancy probability
muZ <- z <- array(dim = c(M, T))        # Expected and realized occurrence
y <- array(NA, dim = c(M, J, T))        # Detection histories

set.seed(13973)
psi[1] <- 0.4                           # Initial occupancy probability
p <- c(0.3,0.4,0.5,0.5,0.1,0.3,0.5,0.5,0.6,0.2)
phi <- runif(n=T-1, min=0.6, max=0.8)   # Survival probability (1-epsilon)
gamma <- runif(n=T-1, min=0.1, max=0.2) # Colonization probability

# Generate latent states of occurrence
# First year
z[,1] <- rbinom(M, 1, psi[1])           # Initial occupancy state
# Later years
for(i in 1:M){                          # Loop over sites
  for(k in 2:T){                        # Loop over years
    muZ[k] <- z[i, k-1]*phi[k-1] + (1-z[i, k-1])*gamma[k-1]
    z[i,k] <- rbinom(1, 1, muZ[k])
  }
}

# Generate detection/non-detection data
for(i in 1:M){
  for(k in 1:T){
    prob <- z[i,k] * p[k]
    for(j in 1:J){
      y[i,j,k] <- rbinom(1, 1, prob)
    }
  }
}

# Compute annual population occupancy
for (k in 2:T){
  psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1]
}

plot(1:T, colMeans(z), type = "b", xlab = "Year",
     ylab = "Proportion of sites occupied",
     col = "black", xlim=c(0.5, 10.5), xaxp=c(1,10,9),
     ylim = c(0,0.6), lwd = 2, lty = 1,
     frame.plot = FALSE, las = 1, pch=16)

psi.app <- colMeans(apply(y, c(1,3), max))
lines(1:T, psi.app, type = "b", col = "blue", lty=3, lwd = 2)
legend(1, 0.6, c("truth", "observed"),
       col=c("black", "blue"), lty=c(1,3), pch=c(16,1))











