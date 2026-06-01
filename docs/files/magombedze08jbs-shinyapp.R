library(shiny)
library(deSolve)
library(ggplot2)
library(reshape2) # For data reshaping for ggplot2

# Define the ODE system for Mtb model (Equations 1-5)
# Source: "magombedze08jbs.pdf" [1]
mtb_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Equations from source [1]
    # Equation (1): Dynamics of resting macrophages (Mr)
    dMr <- sm + pm * Tb * Mr - ki * Tb * Mr - mu_m * Mr
    
    # Equation (2): Dynamics of Mtb infected macrophages (Mib)
    dMib <- ki * Tb * Mr - kb * Mib - ka * (Mib / (Mib + A)) * T - kl * Mib * Cb - mu_mi * Mib
    
    # Equation (3): Dynamics of free bacterial particles (Tb)
    # Nc parameter is used here as per source equation, but not listed in Table 2 [1, 4, 5].
    # Assuming Nc is a bacterial burst size factor for CTL lysis, similar to Nb.
    # Defaulting Nc to Nb for this implementation.
    dTb <- kb * Nb * Mib - gamma1 * Tb * Mr - gamma2 * Tb * Cb + rm * Tb + kl * Nc * Mib * Cb
    
    # Equation (4): Dynamics of health CD4+ T cells (T)
    dT <- st + r1 * (T * Tb / (Tb + B)) - mu_t * T
    
    # Equation (5): Dynamics of Mtb specific CTLs (Cb)
    dCb <- sb + pb * (Mib / (Mib + Rb)) * T * Cb - mu_cb * Cb
    
    list(c(dMr, dMib, dTb, dT, dCb))
  })
}

# Default parameter values based on Table 2 and simulation conditions from source [2-5]
# Note: 'ki' value is from simulation sections, not Table 2.
# 'Nc' is estimated as it's not in Table 2 but present in Eq. (3).
default_parameters <- c(
  sm = 5.0, # Mr source [4]
  pm = 0.000575, # Mr proliferation [4]
  ki = 0.00002, # Mr rate of infection by Mtb (from latency/disease simulations) [2, 3]
  mu_m = 0.011, # Mr death rate [4]
  kb = 0.00001, # Mib burst rate (defaulting to latency value) [2]
  ka = 0.000125, # Apoptosis [4]
  A = 350.0, # CD4+ T Sat Limit (for apoptosis) [4]
  mu_mi = 0.011, # Mib death rate [4]
  kl = 0.02, # Lysis of Mib by CTLs (defaulting to latency value) [2]
  Nb = 50.0, # Mib burst size [5]
  gamma1 = 0.00065, # Tb killing by Mr (defaulting to latency value) [2]
  gamma2 = 0.02, # Tb killing by Cb (defaulting to latency value) [2]
  rm = 0.0625, # Tb Multiplication (defaulting to latency value) [2]
  st = 10.0, # CD4+T supply rate [5]
  r1 = 0.00001, # CD4+ T Proliferation due to Tb [5]
  B = 500.0, # CD4+ Sat Limit due to Tb [5]
  mu_t = 0.02, # CD4+ T death rate [5]
  sb = 5.0, # Mtb CTL supply rate [5]
  pb = 0.000001, # Mtb CTL proliferation [5]
  Rb = 500.0, # CTL proliferation Sat Limit [5]
  mu_cb = 0.95, # Mtb CTL death rate [5]
  Nc = 50.0 # Estimated, not in Table 2. Assumed similar to Nb [1, 10, 17]
)

# Default initial conditions based on Fig. 1/2 in source [2, 3]
default_initial_conditions <- c(
  Mr = 500.0, # Resting macrophages
  Mib = 0.0, # Mtb infected macrophages
  Tb = 20.0, # Mtb pathogen
  T = 1000.0, # Health CD4+ T cells
  Cb = 5.0 # Mtb specific CTLs
)

# UI definition
ui <- fluidPage(
  titlePanel("Mtb Infection Dynamics Model (Equations 1-5)"),
  sidebarLayout(
    sidebarPanel(
      h3("Initial Conditions"),
      numericInput("Mr_init", "Resting Macrophages (Mr)", value = default_initial_conditions["Mr"], min = 0),
      numericInput("Mib_init", "Mtb Infected Macrophages (Mib)", value = default_initial_conditions["Mib"], min = 0),
      numericInput("Tb_init", "Mtb Pathogen (Tb)", value = default_initial_conditions["Tb"], min = 0),
      numericInput("T_init", "Health CD4+ T Cells (T)", value = default_initial_conditions["T"], min = 0),
      numericInput("Cb_init", "Mtb Specific CTLs (Cb)", value = default_initial_conditions["Cb"], min = 0),
      
      h3("Simulation Settings"),
      sliderInput("time_end", "Simulation Duration (Days)", min = 100, max = 5000, value = 2000, step = 100),
      actionButton("run_sim", "Run Simulation"),
      
      h3("Model Parameters (See Table 2 in source for details)"),
      fluidRow(
        column(6,
               numericInput("sm", "sm (Mr source)", value = default_parameters["sm"]),
               numericInput("pm", "pm (Mr proliferation)", value = default_parameters["pm"], step = 1e-6, format = "f"),
               numericInput("ki", "ki (Mr infection rate by Mtb)", value = default_parameters["ki"], step = 1e-6, format = "f"),
               numericInput("mu_m", "µm (Mr death rate)", value = default_parameters["mu_m"]),
               numericInput("kb", "kb (Mib burst rate)", value = default_parameters["kb"], step = 1e-6, format = "f"),
               numericInput("ka", "ka (Apoptosis)", value = default_parameters["ka"], step = 1e-6, format = "f"),
               numericInput("A", "A (CD4+ T Sat Limit)", value = default_parameters["A"]),
               numericInput("mu_mi", "µmi (Mib death rate)", value = default_parameters["mu_mi"]),
               numericInput("kl", "kl (Lysis of Mib by CTLs)", value = default_parameters["kl"]),
               numericInput("Nb", "Nb (Mib burst size)", value = default_parameters["Nb"])
        ),
        column(6,
               numericInput("gamma1", "γ1 (Tb killing by Mr)", value = default_parameters["gamma1"], step = 1e-6, format = "f"),
               numericInput("gamma2", "γ2 (Tb killing by Cb)", value = default_parameters["gamma2"]),
               numericInput("rm", "rm (Tb Multiplication)", value = default_parameters["rm"]),
               numericInput("st", "st (CD4+T supply rate)", value = default_parameters["st"]),
               numericInput("r1", "r1 (CD4+ T Proliferation by Tb)", value = default_parameters["r1"], step = 1e-6, format = "f"),
               numericInput("B", "B (CD4+ Sat Limit due to Tb)", value = default_parameters["B"]),
               numericInput("mu_t", "µt (CD4+ T death rate)", value = default_parameters["mu_t"]),
               numericInput("sb", "sb (Mtb CTL supply rate)", value = default_parameters["sb"]),
               numericInput("pb", "pb (Mtb CTL proliferation)", value = default_parameters["pb"], step = 1e-7, format = "f"),
               numericInput("Rb", "Rb (CTL proliferation Sat Limit)", value = default_parameters["Rb"]),
               numericInput("mu_cb", "µcb (Mtb CTL death rate)", value = default_parameters["mu_cb"]),
               # Nc is from equation (3) but not listed in Table 2, adding with assumed default
               numericInput("Nc", "Nc (Bacteria burst size from CTL lysis)", value = default_parameters["Nc"])
        )
      ),
      p("Note: 'Nc' parameter is present in Eq. (3) of the source [1] but not explicitly defined or given a value in Table 2 [4, 5]. It is estimated here to be similar to 'Nb' (Mib burst size).")
      
    ),
    mainPanel(
      h3("Simulated Population Dynamics Over Time"),
      plotOutput("plot_dynamics", height = "800px")
    )
  )
)

# Server logic
server <- function(input, output) {
  
  # Reactive expression to gather all parameters
  params <- reactive({
    c(
      sm = input$sm, pm = input$pm, ki = input$ki, mu_m = input$mu_m,
      kb = input$kb, ka = input$ka, A = input$A, mu_mi = input$mu_mi,
      kl = input$kl, Nb = input$Nb, gamma1 = input$gamma1, gamma2 = input$gamma2,
      rm = input$rm, st = input$st, r1 = input$r1, B = input$B, mu_t = input$mu_t,
      sb = input$sb, pb = input$pb, Rb = input$Rb, mu_cb = input$mu_cb,
      Nc = input$Nc
    )
  })
  
  # Reactive expression to gather initial conditions
  inits <- reactive({
    c(
      Mr = input$Mr_init, Mib = input$Mib_init, Tb = input$Tb_init,
      T = input$T_init, Cb = input$Cb_init
    )
  })
  
  # Run the simulation when the "Run Simulation" button is clicked
  sim_results <- eventReactive(input$run_sim, {
    times <- seq(0, input$time_end, by = 1)
    ode(y = inits(), times = times, func = mtb_model, parms = params())
  })
  
  # Render the plot
  output$plot_dynamics <- renderPlot({
    # Convert results to a data frame and melt for ggplot2
    df <- as.data.frame(sim_results())
    names(df)[18] <- "Time" # Rename first column to Time for clarity
    df_melted <- melt(df, id.vars = "Time", variable.name = "Population", value.name = "Count")
    
    # Create faceted plot
    ggplot(df_melted, aes(x = Time, y = Count, color = Population)) +
      geom_line(size = 1) +
      facet_wrap(~ Population, scales = "free_y", ncol = 2) + # Use free_y for different scales
      labs(
        title = "Mtb Infection Dynamics (Equations 1-5)",
        x = "Time (Days)",
        y = "Population Count"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
      ) +
      scale_color_brewer(palette = "Set1")
  })
}

# Run the application
shinyApp(ui = ui, server = server)
