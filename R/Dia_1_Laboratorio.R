print("Hola mundo!")

# Día 1: “Costo-efectividad incremental en R”

# scripts/01_icer_ejemplo.R
# Ejemplo simple de costo-efectividad incremental (ICER) en R

# Cargamos tidyverse (útil para manipular data frames)
library(tidyverse)

# 1. Definir datos de ejemplo -------------------------------

# Supongamos dos estrategias:
# A = estándar
# B = nueva intervención

estrategias <- tribble(
  ~estrategia, ~costo,    ~qaly,
  "A",         1000000,   10.0,
  "B",         1200000,   10.3
)

estrategias

# 2. Ordenar por costo y calcular diferencias incrementales ----

resultado_icer <- estrategias %>%
  arrange(costo) %>%                  # por si acaso, ordenamos por costo
  mutate(
    costo_inc = costo - first(costo), # ΔC = C_k - C_estrategia_más_barata
    qaly_inc  = qaly  - first(qaly),  # ΔE = E_k - E_estrategia_más_barata
    icer      = costo_inc / qaly_inc  # ICER = ΔC / ΔE
  )

resultado_icer


# 3. Calcular NMB para distintos umbrales -------------------

# Definimos dos valores hipotéticos de λ (disposición a pagar por QALY)
lambdas <- c(500000, 1000000)

# Función pequeña para calcular NMB dado un λ
calcular_nmb <- function(df, lambda) {
  df %>%
    mutate(
      lambda = lambda,
      nmb    = lambda * qaly - costo
    )
}

# Aplicamos la función a cada λ y combinamos resultados
resultado_nmb <- map_dfr(lambdas, ~calcular_nmb(estrategias, .x))

resultado_nmb

# 4. Ver qué estrategia es óptima para cada λ ----------------

resultado_mejor_por_lambda <- resultado_nmb %>%
  group_by(lambda) %>%
  slice_max(nmb, n = 1, with_ties = FALSE) %>%
  ungroup()

resultado_mejor_por_lambda
