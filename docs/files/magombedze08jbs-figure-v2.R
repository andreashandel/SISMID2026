# Load necessary libraries
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

# Define the system of differential equations from the paper
coinfection_model <- function(t, y, p) {
  # State variables
  Mr <- y[1]
  Mib <- y[2]
  Miv <- y[3]
  Tb <- y[4]
  T_cell <- y[5]
  T_star <- y[6]
  Cb <- y[7]
  Chv <- y[8]
  V <- y[9]
  
  # Parameters
  with(as.list(p), {
    dMr_dt <- sm + pm * (ro * V + Tb) * Mr - (kv * V * Mr) / (1 + ao * Chv) - ki * Tb * Mr - mu_m * Mr
    dMib_dt <- ki * Tb * Mr - kb * Mib - ka * (Mib / (Mib + A)) * T_cell - kl * Mib * Cb - mu_mi * Mib
    dMiv_dt <- (kv * V * Mr) / (1 + ao * Chv) - hm * Miv * Chv - mb * Miv
    dTb_dt <- kb * Nb * Mib - gamma1 * Tb * Mr - gamma2 * Tb * Cb + rm * Tb + kl * Nc * Mib * Cb
    dT_dt <- st + r1 * (T_cell * Tb / (Tb + B)) + r * (T_cell * V / (V + Aa)) - (beta * V * T_cell) / (1 + a1 * Chv) - r2 * (T_cell * V / (T_cell + R)) - mu_t * T_cell
    dT_star_dt <- (beta * V * T_cell) / (1 + a1 * Chv) - hv * T_star * Chv - alpha_t * T_star
    dCb_dt <- sb + pb * (Mib / (Mib + Rb)) * T_cell * Cb - mu_cb * Cb
    dChv_dt <- sv + pv * V * T_cell * Chv - mu_cv * Chv
    dV_dt <- (Nv * alpha_t * T_star) / (1 + a2 * Chv) + (Nm * mb * Miv) / (1 + a3 * Chv) - mu_v * V
    
    # Return the rates of change
    return(list(c(dMr_dt, dMib_dt, dMiv_dt, dTb_dt, dT_dt, dT_star_dt, dCb_dt, dChv_dt, dV_dt)))
  })
}

# Parameters from Table 2 of the paper
parameters <- c(
  sm = 5, pm = 0.000575, ro = 0.1, kv = 0.000025, mu_m = 0.011,
  kb = 0.00002, ka = 0.000125, A = 350, mu_mi = 0.011, hm = 0.000045,
  mb = 0.015, kl = 0.000125, Nc = 20, gamma1 = 0.000165, gamma2 = 0.00032,
  rm = 0.0025, st = 10, r1 = 0.00015, B = 500, r = 0.000125, Aa = 350,
  beta = 0.00025, a1 = 0.005, r2 = 0.000125, R = 500, mu_t = 0.01,
  hv = 0.0025, alpha_t = 0.027, sb = 0.05, pb = 0.00000125, Rb = 25,
  mu_cb = 0.01, sv = 0.05, pv = 0.00000125, mu_cv = 0.01, Nv = 500,
  a2 = 0.005, Nm = 50, a3 = 0.005, mu_v = 2, Nb = 20, ki = 0.00002475,
  ao = 0.05
)

# Initial conditions for Figure 4
# Coinf-1: Latent TB patient co-infected with HIV
y0_1 <- c(Mr = 282.43, Mib = 29.17, Miv = 0.0, Tb = 557.0, T_cell = 1056.0, T_star = 0.0, Cb = 5.23, Chv = 5.0, V = 2.0)
# Coinf-2: HIV patient co-infected with Mtb
y0_2 <- c(Mr = 480.0, Mib = 0.0, Miv = 0.8, Tb = 20.0, T_cell = 200.0, T_star = 5.0, Cb = 5.0, Chv = 10.0, V = 2500.0)

# Time points
t <- seq(from = 0, to = 2000, by = 1)

# Solve the ODEs for both scenarios
solution_1 <- ode(y = y0_1, times = t, func = coinfection_model, parms = parameters) %>% as.data.frame() %>% mutate(Scenario = "Coinf-1")
solution_2 <- ode(y = y0_2, times = t, func = coinfection_model, parms = parameters) %>% as.data.frame() %>% mutate(Scenario = "Coinf-2")

# Combine the results and prepare for plotting
combined_solution <- bind_rows(solution_1, solution_2)

# Define the labels for each plot facet
plot_labels <- c(
  `Mr` = "(a) Resting Macrophages",
  `Mib` = "(b) Mtb Infected Macrophages",
  `Miv` = "(c) HIV Infected Macrophages",
  `Tb` = "(d) Bacterial Load",
  `T_cell` = "(e) CD4+ T Count",
  `T_star` = "(f) HIV Infected CD4+ T Count",
  `Cb` = "(g) Mtb Specific CTLs",
  `Chv` = "(h) HIV Specific CTLs",
  `V` = "(i) Viral Load"
)

# Reshape the data into a long format for ggplot2
solution_long <- combined_solution %>%
  pivot_longer(
    cols = -c(time, Scenario),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  mutate(Variable = factor(Variable, levels = names(plot_labels))) # Ensure correct order

# Create the plot
p1 <- ggplot(solution_long, aes(x = time, y = Value, color = Scenario)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Variable, scales = "free_y", labeller = as_labeller(plot_labels), ncol = 3) +
  labs(
    title = "Figure 4 Recreation: Co-infection Dynamics of HIV and Mtb",
    x = "Time (Days)",
    y = "Population / Concentration",
    color = "Scenario"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

plot(p1)
