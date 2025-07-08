
# Librerias
library(dplyr)
library(car)
library(DHARMa)
library(ggplot2)

#--------------------------

#Abro el DF
datoscrudo <- read.csv("TablaCompleta.csv")
head(datoscrudo)

# Renombrar columnas
datoscrudo <- datoscrudo %>%
  rename(trat = ScTrat, test = ScTest, tiempo = Tiempo)
head (datoscrudo)
#Variables cuantitativas a cualitativas
datoscrudo$trat <- as.factor(datoscrudo$trat)
datoscrudo$test <- as.factor(datoscrudo$test)

#-------------------------
#MODELADO CON CRITERIOS DE ACEPTACIÓN
#-------------------------


#Columna de aceptación con criterio de 4 segundos
datosacepta4 <- datoscrudo
datosacepta4$acepta <- datosacepta4$tiempo > 4
head(datosacepta4)
acepta4 <- datosacepta4[datosacepta4$acepta == TRUE, ]
head(acepta4)

#Visualización
#gráfico de tiempo de ingesta medio con nubes de puntos
library(ggplot2)
library(dplyr)

# Paso 1: Calcular medias y errores estándar
resumen <- acepta3 %>%
  group_by(test, trat) %>%
  summarise(
    media = mean(tiempo, na.rm = TRUE),
    sd = sd(tiempo, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n)
  ) %>%
  ungroup()

# Paso 2: Gráfico
ggplot() +
  # Nubes de puntos
  geom_jitter(
    data = acepta3,
    aes(x = factor(test), y = tiempo, color = factor(trat)),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
    alpha = 0.5,
    size = 2
  ) +
  # Puntos de medias
  geom_point(
    data = resumen,
    aes(x = factor(test), y = media, group= factor(trat)),
    position = position_dodge(width = 0.6),
    size = 3,
    color = "black"
  ) +
  # Barras de error
  geom_errorbar(
    data = resumen,
    aes(
      x = factor(test),
      ymin = media - se,
      ymax = media + se,
      group = factor(trat)
    ),
    position = position_dodge(width = 0.6),
    width = 0.2,
    color = "black"
  ) +
    labs(
      x = "Solución de testeo",
      y = "Tiempo de ingesta (s)",
      color = "Tratamiento (%)",
    ) +
    theme_minimal(base_size = 14)


# Calculamos el porcentaje de aceptación por cada combinación de tratamiento y test
resultados <- datosacepta4 %>%
  group_by(trat, test) %>%
  summarise(porcentaje_aceptacion = mean(acepta == TRUE) * 100)

# Ver resultados
print(resultados)

#Gráfico
ggplot(resultados, aes(x = interaction(trat, test), y = porcentaje_aceptacion, fill = trat)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Combinación de Tratamiento y Test", y = "Porcentaje de aceptación", 
       title = "Porcentaje de aceptación por combinación de Tratamiento y Test") +
  scale_x_discrete(labels = c("trat1%-test1%", "trat1%-test5%", "trat1%-test10%", "trat1%-test20%", 
                              "trat20%-test1%", "trat20%-test5%", "trat20%-test10%", "trat20%-test20%")) +
  theme_minimal()


#Subset de los datos solo para los que aceptan
acepta4 <- datosacepta4[datosacepta4$acepta == TRUE, ]
head(acepta4)

#Gráfico de tiempo para todos los datos
datosacepta4_summary <- datosacepta4 %>%
  group_by(trat, test) %>%
  summarise(
    promedio = mean(tiempo, na.rm = TRUE),
    desv_est = sd(tiempo, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(datosacepta4_summary, aes(x = interaction(trat, test), y = promedio, fill = trat)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = promedio - desv_est, ymax = promedio + desv_est),
                position = position_dodge(width = 0.7), width = 0.25) +
  labs(x = "Combinación de Tratamiento y Test", y = "Promedio de Tiempo", 
       title = "Promedio de Tiempo por Tratamiento y Test") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = function(x) gsub("\\.", "-", x))  

#Gráfico de tiempo solo para los que aceptan
acepta4_summary <- acepta4 %>%
  group_by(trat, test) %>%
  summarise(
    promedio = mean(tiempo, na.rm = TRUE),
    desv_est = sd(tiempo, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(acepta4_summary, aes(x = interaction(trat, test), y = promedio, fill = trat)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = promedio - desv_est, ymax = promedio + desv_est),
                position = position_dodge(width = 0.7), width = 0.25) +
  labs(x = "Combinación de Tratamiento y Test", y = "Promedio de Tiempo", 
       title = "Promedio de Tiempo por Tratamiento y Test") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = function(x) gsub("\\.", "-", x))  

#Chisq 
chisq_test <- chisq.test(table(datosacepta4$trat, datosacepta4$test, datosacepta4$acepta))
print(chisq_test)

#Modelo GLM
#Distribuciones que no son: Norma

hist(acepta4$tiempo, breaks = 30, main = "Histograma de tiempos", xlab = "Tiempo")

#Distribución gamma
m_gamma <- glm(tiempo ~ trat + test, data = acepta4, family = Gamma(link = "log"))
#Sobre o subdispersión
GL.resid <- df.residual(m_gamma)
(dispersion<-sum(resid(m_gamma, type="pearson")^2/df.residual(m_gamma)))
#Chequeo de supuestos con DHARMa
simulationOutput <- simulateResiduals(fittedModel = m_gamma, refit=T, plot = T)
hist(simulationOutput)
# Graficar residuos vs. valores ajustados
plot(m_gamma$fitted.values, residuals(m_gamma))
abline(h = 0, col = "red")

#Gráfico del modelo
library(ggeffects)
library(ggplot2)

# Efectos marginales para trat y test
efectos <- ggpredict(m_gamma, terms = c("test", "trat"))

# Graficar
plot(efectos) +
  labs(
    x = "Solución (test)",
    y = "Tiempo predicho",
    color = "Tratamiento (%)",
  ) +
  theme_minimal()

#---------------------
#Columna de aceptación con criterio de 2 segundos
datosacepta2 <- datoscrudo
datosacepta2$acepta <- datosacepta2$tiempo > 2
head(datosacepta2)

#Subset de los datos solo para los que aceptan
acepta2 <- datosacepta2[datosacepta2$acepta == TRUE, ]
head(acepta2)

#Visualización

# Calculamos el porcentaje de aceptación por cada combinación de tratamiento y test
resultados <- datosacepta2 %>%
  group_by(trat, test) %>%
  summarise(porcentaje_aceptacion = mean(acepta == TRUE) * 100)

# Ver resultados
print(resultados)

#Gráfico
ggplot(resultados, aes(x = interaction(trat, test), y = porcentaje_aceptacion, fill = trat)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Combinación de Tratamiento y Test", y = "Porcentaje de aceptación", 
       title = "Porcentaje de aceptación por combinación de Tratamiento y Test") +
  scale_x_discrete(labels = c("trat1%-test1%", "trat1%-test5%", "trat1%-test10%", "trat1%-test20%", 
                              "trat20%-test1%", "trat20%-test5%", "trat20%-test10%", "trat20%-test20%")) +
  theme_minimal()

#Gráfico de tiempo para todos los datos
datosacepta2_summary <- datosacepta2 %>%
  group_by(trat, test) %>%
  summarise(
    promedio = mean(tiempo, na.rm = TRUE),
    desv_est = sd(tiempo, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(datosacepta2_summary, aes(x = interaction(trat, test), y = promedio, fill = trat)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = promedio - desv_est, ymax = promedio + desv_est),
                position = position_dodge(width = 0.7), width = 0.25) +
  labs(x = "Combinación de Tratamiento y Test", y = "Promedio de Tiempo", 
       title = "Promedio de Tiempo por Tratamiento y Test") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = function(x) gsub("\\.", "-", x))  

#Gráfico de tiempo solo para los que aceptan
acepta2_summary <- acepta2 %>%
  group_by(trat, test) %>%
  summarise(
    promedio = mean(tiempo, na.rm = TRUE),
    desv_est = sd(tiempo, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(acepta2_summary, aes(x = interaction(trat, test), y = promedio, fill = trat)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = promedio - desv_est, ymax = promedio + desv_est),
                position = position_dodge(width = 0.7), width = 0.25) +
  labs(x = "Combinación de Tratamiento y Test", y = "Promedio de Tiempo", 
       title = "Promedio de Tiempo por Tratamiento y Test") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = function(x) gsub("\\.", "-", x))

#Modelo GLM
#Distribuciones que no son: Normal

#Distribución gamma
m_gamma <- glm(tiempo ~ trat + test, data = acepta2, family = Gamma(link = "log"))
#Sobre o subdispersión
GL.resid <- df.residual(m_gamma)
(dispersion<-sum(resid(m_gamma, type="pearson")^2/df.residual(m_gamma)))
#Chequeo de supuestos con DHARMa
simulationOutput <- simulateResiduals(fittedModel = m_gamma, refit=T, plot = T)
hist(simulationOutput)



#---------------------
#Columna de aceptación con criterio de 3 segundos
datosacepta3 <- datoscrudo
datosacepta3$acepta <- datosacepta3$tiempo > 3
head(datosacepta3)

#Subset de los datos solo para los que aceptan
acepta3 <- datosacepta3[datosacepta3$acepta == TRUE, ]
head(acepta3)

#Modelo GLM
#Distribuciones que no son: Normal

#Distribución gamma
m_gamma <- glm(tiempo ~ trat + test, data = acepta3, family = Gamma(link = "log"))
#Sobre o subdispersión
GL.resid <- df.residual(m_gamma)
(dispersion<-sum(resid(m_gamma, type="pearson")^2/df.residual(m_gamma)))
#Chequeo de supuestos con DHARMa
simulationOutput <- simulateResiduals(fittedModel = m_gamma, refit=T, plot = T)
hist(simulationOutput)


#Visualización

# Calculamos el porcentaje de aceptación por cada combinación de tratamiento y test
resultados <- datosacepta3 %>%
  group_by(trat, test) %>%
  summarise(porcentaje_aceptacion = mean(acepta == TRUE) * 100)

# Ver resultados
print(resultados)

#Gráfico
ggplot(resultados, aes(x = interaction(trat, test), y = porcentaje_aceptacion, fill = trat)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Combinación de Tratamiento y Test", y = "Porcentaje de aceptación", 
       title = "Porcentaje de aceptación por combinación de Tratamiento y Test") +
  scale_x_discrete(labels = c("trat1%-test1%", "trat1%-test5%", "trat1%-test10%", "trat1%-test20%", 
                              "trat20%-test1%", "trat20%-test5%", "trat20%-test10%", "trat20%-test20%")) +
  theme_minimal()

#Gráfico de tiempo para todos los datos
datosacepta3_summary <- datosacepta3 %>%
  group_by(trat, test) %>%
  summarise(
    promedio = mean(tiempo, na.rm = TRUE),
    desv_est = sd(tiempo, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(datosacepta3_summary, aes(x = interaction(trat, test), y = promedio, fill = trat)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = promedio - desv_est, ymax = promedio + desv_est),
                position = position_dodge(width = 0.7), width = 0.25) +
  labs(x = "Combinación de Tratamiento y Test", y = "Promedio de Tiempo", 
       title = "Promedio de Tiempo por Tratamiento y Test") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = function(x) gsub("\\.", "-", x))  

#Gráfico de tiempo solo para los que aceptan
acepta3_summary <- acepta3 %>%
  group_by(trat, test) %>%
  summarise(
    promedio = mean(tiempo, na.rm = TRUE),
    desv_est = sd(tiempo, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(acepta3_summary, aes(x = interaction(trat, test), y = promedio, fill = trat)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = promedio - desv_est, ymax = promedio + desv_est),
                position = position_dodge(width = 0.7), width = 0.25) +
  labs(x = "Combinación de Tratamiento y Test", y = "Promedio de Tiempo", 
       title = "Promedio de Tiempo por Tratamiento y Test") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = function(x) gsub("\\.", "-", x))


#---------------------------
#Columna de aceptación con criterio de 5 segundos
datosacepta5 <- datoscrudo
datosacepta5$acepta <- datosacepta5$tiempo > 5
head(datosacepta5)

#Subset de los datos solo para los que aceptan
acepta5 <- datosacepta5[datosacepta5$acepta == TRUE, ]
head(acepta5)

#Visualización

# Calculamos el porcentaje de aceptación por cada combinación de tratamiento y test
resultados <- datosacepta5 %>%
  group_by(trat, test) %>%
  summarise(porcentaje_aceptacion = mean(acepta == TRUE) * 100)

# Ver resultados
print(resultados)

#Gráfico
ggplot(resultados, aes(x = interaction(trat, test), y = porcentaje_aceptacion, fill = trat)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Combinación de Tratamiento y Test", y = "Porcentaje de aceptación", 
       title = "Porcentaje de aceptación por combinación de Tratamiento y Test") +
  scale_x_discrete(labels = c("trat1%-test1%", "trat1%-test5%", "trat1%-test10%", "trat1%-test20%", 
                              "trat20%-test1%", "trat20%-test5%", "trat20%-test10%", "trat20%-test20%")) +
  theme_minimal()

#Gráfico de tiempo para todos los datos
datosacepta5_summary <- datosacepta5 %>%
  group_by(trat, test) %>%
  summarise(
    promedio = mean(tiempo, na.rm = TRUE),
    desv_est = sd(tiempo, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(datosacepta5_summary, aes(x = interaction(trat, test), y = promedio, fill = trat)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = promedio - desv_est, ymax = promedio + desv_est),
                position = position_dodge(width = 0.7), width = 0.25) +
  labs(x = "Combinación de Tratamiento y Test", y = "Promedio de Tiempo", 
       title = "Promedio de Tiempo por Tratamiento y Test") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = function(x) gsub("\\.", "-", x)) 

#Gráfico de tiempo para los que aceptan
acepta5_summary <- acepta5 %>%
  group_by(trat, test) %>%
  summarise(
    promedio = mean(tiempo, na.rm = TRUE),
    desv_est = sd(tiempo, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(acepta5_summary, aes(x = interaction(trat, test), y = promedio, fill = trat)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = promedio - desv_est, ymax = promedio + desv_est),
                position = position_dodge(width = 0.7), width = 0.25) +
  labs(x = "Combinación de Tratamiento y Test", y = "Promedio de Tiempo", 
       title = "Promedio de Tiempo por Tratamiento y Test") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = function(x) gsub("\\.", "-", x)) 

#Modelo GLM
#Distribuciones que no son: Normal

#Distribución gamma
m_gamma <- glm(tiempo ~ trat + test, data = acepta5, family = Gamma(link = "log"))
#Sobre o subdispersión
GL.resid <- df.residual(m_gamma)
(dispersion<-sum(resid(m_gamma, type="pearson")^2/df.residual(m_gamma)))
#Chequeo de supuestos con DHARMa
simulationOutput <- simulateResiduals(fittedModel = m_gamma, refit=T, plot = T)
hist(simulationOutput)

#---------
### Test no parametrico? ULTIMA OPCION

kruskal.test(tiempo ~ trat, data = acepta4)
kruskal.test(tiempo ~ test, data = acepta4)


##Comparaciones

pairwise.wilcox.test(datos$Tiempo, datosGamma$ScTrat,
                     p.adjust.method = "BH")

#--------
#Columna de aceptación con criterio de 7 segundos
datosacepta7 <- datoscrudo
datosacepta7$acepta <- datosacepta7$tiempo > 7
head(datosacepta7)

#Subset de los datos solo para los que aceptan
acepta7 <- datosacepta7[datosacepta7$acepta == TRUE, ]
head(acepta7)

#Visualización

# Calculamos el porcentaje de aceptación por cada combinación de tratamiento y test
resultados <- datosacepta7 %>%
  group_by(trat, test) %>%
  summarise(porcentaje_aceptacion = mean(acepta == TRUE) * 100)

# Ver resultados
print(resultados)

#Gráfico
ggplot(resultados, aes(x = interaction(trat, test), y = porcentaje_aceptacion, fill = trat)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Combinación de Tratamiento y Test", y = "Porcentaje de aceptación", 
       title = "Porcentaje de aceptación por combinación de Tratamiento y Test") +
  scale_x_discrete(labels = c("trat1%-test1%", "trat1%-test5%", "trat1%-test10%", "trat1%-test20%", 
                              "trat20%-test1%", "trat20%-test5%", "trat20%-test10%", "trat20%-test20%")) +
  theme_minimal()

#Gráfico de tiempo para todos los datos
datosacepta7_summary <- datosacepta7 %>%
  group_by(trat, test) %>%
  summarise(
    promedio = mean(tiempo, na.rm = TRUE),
    desv_est = sd(tiempo, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(datosacepta7_summary, aes(x = interaction(trat, test), y = promedio, fill = trat)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = promedio - desv_est, ymax = promedio + desv_est),
                position = position_dodge(width = 0.7), width = 0.25) +
  labs(x = "Combinación de Tratamiento y Test", y = "Promedio de Tiempo", 
       title = "Promedio de Tiempo por Tratamiento y Test") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = function(x) gsub("\\.", "-", x)) 

#Gráfico de tiempo para los que aceptan
acepta7_summary <- acepta7 %>%
  group_by(trat, test) %>%
  summarise(
    promedio = mean(tiempo, na.rm = TRUE),
    desv_est = sd(tiempo, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(acepta7_summary, aes(x = interaction(trat, test), y = promedio, fill = trat)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = promedio - desv_est, ymax = promedio + desv_est),
                position = position_dodge(width = 0.7), width = 0.25) +
  labs(x = "Combinación de Tratamiento y Test", y = "Promedio de Tiempo", 
       title = "Promedio de Tiempo por Tratamiento y Test") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = function(x) gsub("\\.", "-", x)) 



#Modelo GLM
#Distribuciones que no son: Normal

#Distribución gamma
m_gamma <- glm(tiempo ~ trat + test, data = acepta7, family = Gamma(link = "log"))
#Sobre o subdispersión
GL.resid <- df.residual(m_gamma)
(dispersion<-sum(resid(m_gamma, type="pearson")^2/df.residual(m_gamma)))
#Chequeo de supuestos con DHARMa
simulationOutput <- simulateResiduals(fittedModel = m_gamma, refit=T, plot = T)
hist(simulationOutput)


#-----------------------------
#PRUEBAS DE CHI CUADRADO
#-----------------------------

#Para criterio de aceptación 3

library(dplyr)
library(purrr)

# Calcular p-valores del chisq.test por cada concentración
pvalores_por_test3 <- datosacepta3 %>%
  group_by(test) %>%
  summarise(
    chisq = list(chisq.test(table(trat, acepta))),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = map_dbl(chisq, ~ .x$p.value),
    p_text = ifelse(p_value < 0.001, "p < 0.001", paste0("p = ", signif(p_value, 3)))
  ) %>%
  select(test, p_value, p_text)

#Para er las tablas de contingencia que usa
tablas_y_pvalores <- datosacepta3 %>%
  group_by(test) %>%
  summarise(
    tabla = list(table(trat, acepta)),  # Guarda la tabla
    chisq = list(chisq.test(tabla[[1]])),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = map_dbl(chisq, ~ .x$p.value),
    p_text = ifelse(p_value < 0.001, "p < 0.001", paste0("p = ", signif(p_value, 3)))
  )

#Por si acaso fisher
tabla_20 <- matrix(c(1, 30, 2, 29), nrow = 2, byrow = TRUE)
fisher.test(tabla_20)

# Ver resultados
print(pvalores_por_test3)

# Calcular chisq por cada tratamiento
pvalores_por_trat <- datosacepta3 %>%
  group_by(trat) %>%
  summarise(
    chisq = list(chisq.test(table(test, acepta))),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = map_dbl(chisq, ~ .x$p.value),
    p_text = ifelse(p_value < 0.001, "p < 0.001", paste0("p = ", signif(p_value, 3)))
  ) %>%
  select(trat, p_value, p_text)

#Haciendo las tablas de contingencia
# Crear las tablas y guardarlas junto al p-valor
tablas_por_trat3 <- datosacepta3 %>%
  group_by(trat) %>%
  summarise(
    tabla = list(table(test, acepta)),      # guarda la tabla
    chisq = list(chisq.test(tabla[[1]])),   # hace el test con esa tabla
    .groups = "drop"
  ) %>%
  mutate(
    p_value = map_dbl(chisq, ~ .x$p.value),
    p_text = ifelse(p_value < 0.001, "p < 0.001", paste0("p = ", signif(p_value, 3)))
  )
#Immprimir las tablas
for (i in seq_along(tablas_por_trat3$trat)) {
  cat("\n--- Tratamiento:", tablas_por_trat3$trat[i], "---\n")
  print(tablas_por_trat3$tabla[[i]])
  cat("p-value:", tablas_por_trat3$p_text[i], "\n")
}


# Ver resultado
print(pvalores_por_trat)

#Para 4 segundos 

# Calcular p-valores del chisq.test por cada concentración
pvalores_por_test4 <- datosacepta4 %>%
  group_by(test) %>%
  summarise(
    chisq = list(chisq.test(table(trat, acepta))),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = map_dbl(chisq, ~ .x$p.value),
    p_text = ifelse(p_value < 0.001, "p < 0.001", paste0("p = ", signif(p_value, 3)))
  ) %>%
  select(test, p_value, p_text)

# Ver resultados
print(pvalores_por_test4)

tablas_y_pvalores4 <- datosacepta4 %>%
  group_by(test) %>%
  summarise(
    tabla = list(table(trat, acepta)),  # Guarda la tabla
    chisq = list(chisq.test(tabla[[1]])),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = map_dbl(chisq, ~ .x$p.value),
    p_text = ifelse(p_value < 0.001, "p < 0.001", paste0("p = ", signif(p_value, 3)))
  )

tablas_y_pvalores4$tabla[[which(tablas_y_pvalores$test == 20)]]

#Por si acaso fisher
tabla_20 <- matrix(c(1, 30, 5, 26), nrow = 2, byrow = TRUE)
fisher.test(tabla_20)

# Calcular chisq por cada tratamiento
pvalores_por_trat4 <- datosacepta4 %>%
  group_by(trat) %>%
  summarise(
    chisq = list(chisq.test(table(test, acepta))),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = map_dbl(chisq, ~ .x$p.value),
    p_text = ifelse(p_value < 0.001, "p < 0.001", paste0("p = ", signif(p_value, 3)))
  ) %>%
  select(trat, p_value, p_text)

# Ver resultado
print(pvalores_por_trat4)

#--------------------------
#Gráficos de aceptación con chisq

library(dplyr)

resumen_grafico <- datosacepta4 %>%
  group_by(test, trat) %>%
  summarise(
    n_total = n(),
    n_acepta = sum(acepta),
    porcentaje_acepta = 100 * n_acepta / n_total,
    .groups = "drop"
  )

resumen_grafico <- resumen_grafico %>%
  left_join(pvalores_por_test4, by = "test")
resumen_grafico <- resumen_grafico %>%
  mutate(p_asteriscos = case_when(
    p_value <= 0.001 ~ "***",
    p_value <= 0.01 ~ "**",
    p_value <= 0.05 ~ "*",
    TRUE ~ ""
  ))

library(ggplot2)

ggplot(resumen_grafico, aes(x = factor(test), y = porcentaje_acepta, fill = factor(trat))) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    aes(label = p_asteriscos, y = porcentaje_acepta + 3),  # +3 para que esté justo arriba de la barra
    position = position_dodge(width = 0.8),
    size = 6,
    color = "black"
  ) +
  scale_fill_brewer(palette = "Set2", name = "Tratamiento (%)") +
  labs(x = "solución de testeo", y = "% de aceptación") +
  theme_minimal(base_size = 14)

#--------------------------
#Tablade aceptación por tratamientos 
tabla_resumen <- datosacepta3 %>%
  group_by(trat, test) %>%
  summarise(
    total = n(),
    aceptan = sum(acepta == TRUE),
    porcentaje = round(100 * aceptan / total, 1),
    .groups = "drop"
  )

tabla_resumen
