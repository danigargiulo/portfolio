# Efectos del estradiol prepuberal en comportamiento obsesivo-compulsivo en ratones

Este proyecto fue realizado en el 2023 como colaboración con un proyecto de investigación en neurociencias.
Analiza el efecto del tratamiento hormonal con estradiol durante la pubertad en ratones *Mus musculus*, evaluando su comportamiento obsesivo-compulsivo a través del **Marble Burying Test**. Se aplicaron modelos estadísticos mixtos para evaluar el impacto del tratamiento, el sexo y la camada sobre la dinámica del entierro de canicas.

---

## Herramientas estadísticas utilizadas

- **Lenguaje:** R
- **Librerías principales:**
  - lme4, glmmTMB – modelos lineales mixtos
  - car, nlme, performance, DHARMa – evaluación de supuestos y efectos
  - ggplot2, corrplot, ggcorrplot – visualización
  - pastecs, reshape2 – manipulación y resumen de datos

---

## Diseño experimental - Dataset

- **Variable respuesta:** Tiempo (min) en que cada ratón entierra el 35% o 50% de las canicas
- **Variables explicativas:**
  - **Tratamiento:** Control vs. Estradiol
  - **Sexo:** Macho vs. Hembra
  - **Camada:** Factor aleatorio (6 niveles)
- **Diseño:** Bloques al azar, con **Camada como bloque aleatorio**
- **Unidad experimental:** cada ratón
- **Test conductual utilizado:** Marble Burying Test

---

## Modelos aplicados

- **GLMM** con distribución Gamma y enlace logarítmico
- Comparaciones de medias y análisis de interacciones
- Evaluación de supuestos con DHARMa y Shapiro-Wilk
- Análisis del efecto fijo de tratamiento y sexo, y componente aleatorio de camada
