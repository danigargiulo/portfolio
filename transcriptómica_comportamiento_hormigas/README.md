
# Proyectos de análisis transcriptómico y de comportamiento 

Este repositorio contiene dos análisis independientes realizados en el marco de proyectos de investigación en biología experimental y molecular, analizando perfil transcriptómico y comportamental de *Linepithema humile*.  Ambos estudios están en desarrollo y vinculados a trabajos en proceso de publicación científica.

---
## Parte 1: Análisis transcriptómico en *Linepithema humile*

**Resumen:**  
Análisis de expresión génica diferencial entre condiciones experimentales mediante RNA-seq. Incluye filtrado de genes, transformación rlog, análisis con DESeq2 y exploración de enriquecimiento funcional con GO terms (usando ermineR).

**Herramientas utilizadas:**  
- R (DESeq2, HTSFilter, pheatmap, ggpubr, ermineR, rJava)
- Transformación rlog, PCA, MA-plot, análisis GO

**Datos:**  
Matriz de conteos crudos (Raw_counts.txt) y diseño experimental (experimental_design.txt).

---
## Parte 2: Análisis de comportamiento de ingesta

**Resumen:**  
Se analizó el comportamiento de aceptación de soluciones azucaradas en función de distintos tratamientos y concentraciones de testeo. El análisis incluyó visualizaciones, modelado estadístico (GLM con distribución Gamma), validación de supuestos y análisis de proporciones.

**Herramientas utilizadas:**  
- R (dplyr, ggplot2, car, DHARMa, ggeffects, purrr)
- Modelos GLM, test de Chi-cuadrado, tests no paramétricos

**Datos:**  
Tabla con registros individuales de tiempo de ingesta bajo distintas combinaciones de tratamiento y concentración de testeo.

---

## Estado del proyecto

Ambos análisis están en curso y forman parte de colaboraciones en investigación científica. Algunos archivos pueden no estar completamente depurados ya que los estudios están siendo preparados para publicación.
