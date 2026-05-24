############################################################
##  Magombedze et al. (2008) – Equations 25 – 33
##  Parameter set: Table 2        ICs: Figure 4 caption
############################################################

# ---------- packages ----------
library(deSolve)
library(ggplot2)
library(patchwork)
library(tidyr)

# ---------- parameters (Table 2) ----------
pars <- c(
  # resting‑macrophage dynamics
  s_m   = 5,          # Mr source
  p_m   = 0.000575,   # Mr proliferation
  r_o   = 0.1,        # HIV‑induced prolif.
  k_v   = 0.000025,   # HIV→Mr infection rate
  mu_m  = 0.011,      # Mr death
  k_i   = 0.000025,   # Mtb→Mr infection rate  (Table has kv=0.000025; ki identical value)
  k_l   = 0.005,      # CTL lysis of Mib
  a_o   = 1,          # scaling (assumed 1; not in table but needed)
  
  # infected macrophages
  k_b   = 0.00002,    # Mib burst
  k_a   = 0.000125,   # apoptosis of Mib
  A     = 350,        # CD4 sat limit (apoptosis term)
  mu_mi = 0.011,      # infected‑mac death
  h_m   = 0.000045,   # lytic kill of Miv by CTL
  m_b   = 0.015,      # Miv burst (virus)
  
  # bacilli & burst sizes
  N_b   = 50,         # burst size (Mib → Tb)
  N_c   = 10,         # CTL‑lysis burst size **(not in Table 2; assumed equal to N_b)**  ## <<< adjust if known
  
  # extracellular Mtb
  gamma1 = 0.00065,   # Tb kill by Mr
  gamma2 = 0.02,      # Tb kill by Cb
  r_m    = 0.1,       # Tb multiplication
  
  # CD4 T‑cell kinetics
  s_t  = 10,
  r    = 0.1,   # HIV‑induced CD4 prolif. **(not listed; placeholder)**  ## <<< adjust if known
  r1   = 0.00001,   # Tb‑induced prolif.
  r2   = 0.00001,   # apoptosis term
  A_a  = 350,       # HIV sat limit
  B    = 500,       # Tb sat limit
  beta = 0.00025,   # HIV infection rate
  a1   = 0.075,     # infection inhibition factor
  R    = 350,       # apoptosis sat limit
  mu_t = 0.02,      # CD4 death
  
  # HIV‑infected CD4
  h_v     = 0.0025,
  alpha_t = 0.025,
  
  # HIV‑specific CTLs
  s_v  = 5,
  p_v  = 0.00001,
  mu_cv = 1.5,
  
  # Mtb‑specific CTLs
  s_b  = 5,
  p_b  = 0.000001,
  mu_cb = 0.95,
  R_b   = 500,
  
  # free virion production / clearance
  N_v  = 1000,
  N_m  = 800,
  a2   = 0.05,
  a3   = 0.85,
  mu_v = 1.5
)

# ---------- initial conditions (Figure 4 caption) ----------
state <- c(
  M_r    = 500.0,
  M_ib   = 0.0,
  M_iv   = 0.0,
  T_b    = 20.0,
  T      = 1000.0,
  T_star = 0.0,      # HIV‑infected CD4
  C_hv   = 5.0,
  C_b    = 5.0,
  V      = 10.0
)

# ---------- model ODEs ----------
coinf_ode <- function(t, y, p) {
  with(as.list(c(y, p)), {
    
    # resting macrophages
    dM_r <- s_m + p_m * (r_o*V + T_b) * M_r -
      k_v*V*M_r/(1 + a_o*C_hv) -
      k_i*T_b*M_r - mu_m*M_r
    
    # Mtb‑infected macrophages
    dM_ib <- k_i*T_b*M_r - k_b*M_ib -
      k_a*(M_ib/(M_ib + A))*T -
      k_l*M_ib*C_b - mu_mi*M_ib
    
    # HIV‑infected macrophages
    dM_iv <- k_v*V*M_r/(1 + a_o*C_hv) -
      h_m*M_iv*C_hv - m_b*M_iv
    
    # extracellular bacilli
    dT_b <- k_b*N_b*M_ib -
      gamma1*T_b*M_r -
      gamma2*T_b*C_b +
      r_m*T_b +
      k_l*N_c*M_ib*C_b
    
    # healthy CD4 T cells
    dT <- s_t +
      r  *T*V /(V + A_a) +
      r1 *T*T_b/(T_b + B) -
      beta*V*T/(1 + a1*C_hv) -
      r2 *T*V/(T + R) -
      mu_t*T
    
    # HIV‑infected CD4 T cells
    dT_star <- beta*V*T/(1 + a1*C_hv) -
      h_v*T_star*C_hv -
      alpha_t*T_star
    
    # HIV‑specific CTLs
    dC_hv <- s_v + p_v*V*T*C_hv - mu_cv*C_hv
    
    # Mtb‑specific CTLs
    dC_b <- s_b + p_b*M_ib/(M_ib + R_b)*T*C_b - mu_cb*C_b
    
    # free HIV virions
    dV <- N_v*alpha_t*T_star/(1 + a2*C_hv) +
      N_m*m_b*M_iv     /(1 + a3*C_hv) -
      mu_v*V
    
    list(c(dM_r, dM_ib, dM_iv, dT_b, dT, dT_star, dC_hv, dC_b, dV))
  })
}

# ---------- integration ----------
times <- seq(0, 600, by = 1)   # match ~ Fig 4 horizon
out   <- as.data.frame(ode(state, times, coinf_ode, pars))

# ---------- plotting helper ----------
labels <- c(
  M_r    = "(a) Resting macrophages",
  M_ib   = "(b) Mtb‑inf. macrophages",
  M_iv   = "(c) HIV‑inf. macrophages",
  T_b    = "(d) Extracellular Mtb",
  T      = "(e) CD4 T cells",
  T_star = "(f) HIV‑inf. CD4",
  C_hv   = "(g) HIV‑CTLs",
  C_b    = "(h) Mtb‑CTLs",
  V      = "(i) Free HIV virions"
)

plot_df <- out |>
  pivot_longer(-time, names_to = "var", values_to = "value") |>
  mutate(label = factor(labels[var], levels = labels))

# log‑scale for high‑range species (adjust as needed)
log_vars <- c("(d) Extracellular Mtb",
              "(g) HIV‑CTLs", "(h) Mtb‑CTLs",
              "(i) Free HIV virions")

plot_df <- plot_df |>
  mutate(value_plot = ifelse(label %in% log_vars,
                             log10(pmax(value, 1e-6)),
                             value))

panel <- function(df) {
  ggplot(df, aes(time, value_plot)) +
    geom_line(size = 0.35) +
    labs(x = "Time (days)", y = NULL) +
    ggtitle(unique(df$label)) +
    theme_minimal(base_size = 8) +
    theme(
      plot.title = element_text(size = 9, face = "bold"),
      axis.text  = element_text(size = 6),
      axis.title.x = element_text(size = 7),
      panel.grid.minor = element_blank()
    )
}

plots <- plot_df |>
  group_split(label) |>
  lapply(panel)

fig4 <- (plots[[1]] | plots[[2]] | plots[[3]]) /
  (plots[[4]] | plots[[5]] | plots[[6]]) /
  (plots[[7]] | plots[[8]] | plots[[9]])

# ---------- save ----------
ggsave("figure4_reproduction.png", fig4,
       width = 180, height = 180, units = "mm", dpi = 300)

# View in RStudio / R GUI
print(fig4)
