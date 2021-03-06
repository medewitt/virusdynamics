#' # Chapter 4 - Anti-viral Drug Therapy
#' pg 35-37

library(deSolve)
library(tidyverse)

#' Set Up the ODE
#' Represents infected cells becoming:
#' * Latently infected: do not produce new virions, but
#'    contain replication competent virus that can
#'    be re-activated to become virus producing
#' * long-lived chronic producers: produce small amounts
#'    of virus over long periods
#' * cells which harbor defective provirus
base_ode <- function(time, state, parameters){
  with(as.list(c(state, parameters)),{
    
    dx <- lambda - d*x - beta * x * v
    dy1 <- q1*beta * x * v - a1 * y1 + alpha * y2
    dy2 <- q2*beta * x * v - a2 * y2 - alpha * y2
    dy3 <- q3*beta * x * v - a3 * y3
    dv <- k * y1 - u * v
    
    return(list(c(dx, dy1, dy2, dy3, dv)))
  })
  
}

t <- seq(0,30,.1)

params <- c(
  lambda =1e7,# Uninfected cell production rate
  d = .1, # Cell Death Rate
  a1 = .5, # Infected Cell Death Rate
  a2 = .01, # Infected Cell Death Rate
  a3 = .008, # Infected Cell Death Rate
  beta = 5e-10, # "Rate Constant"
  alpha = .3, # Virus Production rate of reactivated latent
  k = 500, # Virus productin from Infected cell
  u = 5, # Free Virus lifestapn
  q1 = .55, # P(Infected State | Infected)
  q2 = .05, # P(Latent State | Infected)
  q3 = .4 # P(Defective Provirus | Infected)
)

#' Guessing Initial values from a graph
x0 <- params["lambda"][1]/params["d"][1]
init <- c(x = unname(x0), 
          y1 = 1,y2=0,y3=0, v = 1)

out <- ode(init, t, base_ode, params)

out_df <- as_tibble(as.data.frame(out))

#' Funny side note is that the differences in the scales
#' are almost immediately reproduced as in the book

compartment_names <- c("Uninfected Cells", 
"Actively Infected Cells", 
"Latently Infected Cells",
"Cells Infected with Defective Virus",
"Free Virus")

p <- out_df %>% 
  setNames(c("time", compartment_names)) %>% 
  gather(compartment, value, -time) %>% 
  mutate(compartment = factor(compartment,
                              compartment_names)) %>% 
  ggplot(aes(time, y = value))+
  geom_line()+
  facet_wrap(~compartment, ncol = 1,scales = "free_y")+
  theme_classic()+
  labs(
    title = "Viral Dynamics of Primary Phase of HIV or SIV Infection",
    x = "Days",
    y = "Count"
  )+
  scale_y_continuous(labels = scales::comma_format(accuracy = 1000))+
  theme(strip.background = element_blank())

p	

cowplot::save_plot(filename = here::here("Ch4", "Figure4-5.png"),
                   plot = p)
