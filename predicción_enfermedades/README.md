# Predicción de Enfermedades Cardiovasculares

Este proyecto fue realizado en el marco del curso de Data Science de la UTN durante el año 2024.  
Tiene como objetivo construir un modelo predictivo a través de machine learning capaz de estimar la presencia de enfermedades cardiovasculares a partir de variables clínicas fácilmente accesibles, como presión arterial, colesterol, edad, y resultados de ECG.

Las enfermedades cardiovasculares son una de las principales causas de muerte en el mundo, representando aproximadamente un 32% de las muertes a nivel global. Modelos predictivos como el desarrollado aquí pueden ayudar a identificar individuos en riesgo y permitir intervenciones preventivas más eficaces.

---

## Tecnologías y herramientas utilizadas

- **Lenguaje:** Python
- **Librerías:**  
  - Manipulación de datos: pandas, numpy  
  - Visualización: matplotlib, seaborn  
  - Modelado: scikit-learn (logistic regression, decision tree, random forest)
  - Evaluación: confusion_matrix, classification_report, roc_auc_score
- **Entorno:** Jupyter Notebook

---

## Dataset utilizado

El dataset contiene registros clínicos de pacientes, con variables como:
- Edad, sexo
- Tipo de dolor de pecho
- Presión sanguínea en reposo
- Colesterol sérico
- Frecuencia cardíaca máxima
- Resultados de electrocardiogramas
- Ejercicio inducido y síntomas asociados
- Talasemia
- Diagnóstico final (presencia o ausencia de enfermedad)

---

## Objetivos del proyecto

- Preprocesamiento y análisis exploratorio de datos clínicos
- Análisis de correlaciones y visualización de variables más relevantes
- Implementación de modelos de clasificación
- Evaluación comparativa de desempeño
- Interpretación de resultados en contexto médico-preventivo
