# Proyecto: Análisis de Similitud entre gen o CDS y Genoma

Este repositorio contiene un código desarrollado en Python para analizar la similitud entre un gen o CDS (secuencia codificante de ADN) y un genoma completo. El código permite realizar comparaciones entre las secuencias en diferentes orientaciones y generar varias visualizaciones para interpretar los resultados.

---

## Características del Proyecto

### 1. **Funciones principales**

- **Cálculo de similitudes:**
  - La función `similarities` permite calcular la similitud entre una ventana deslizante del genoma y el CDS.
- **Transformaciones de ADN:**
  - Complemento (`complement`).
  - Reverso (`reverse`).
  - Reverso complementario (`reverse_complement`).

### 2. **Visualizaciones generadas**

El código genera las siguientes gráficas:

1. **Similitud Forward:** Comparación de la similitud entre el CDS y el genoma en su orientación directa.
2. **Similitud Reverse:** Comparación de la similitud entre el CDS y el genoma en su orientación inversa.
3. **Comparación de Similitudes Forward y Reverse:** Gráfica que superpone ambas orientaciones para su análisis conjunto.
4. **Destacando máximos:** Visualización de las posiciones con similitud máxima en cada orientación.
5. **Resaltando valores sobre un umbral:** Muestra las posiciones donde la similitud supera un umbral definido.
6. **Distribución de Similitudes:** Histogramas que comparan las distribuciones forward y reverse.
7. **Mapa de calor:** Representa las similitudes como un mapa de calor para facilitar su interpretación visual.
8. **Similitud acumulativa:** Análisis de distribuciones acumulativas para forward y reverse.

---

## Requisitos

- **Python 3.7 o superior**
- Bibliotecas necesarias:
  - `pandas`
  - `seaborn`
  - `matplotlib`


Puedes instalarlas con:

```bash
pip install -r requirements.txt
```

---

## Archivos

1. **`genome_sequence_matching.py`** : Archivo principal con el código para cargar los datos, calcular similitudes y generar las visualizaciones.
2. **`solanum_phureja_chromosome1.fasta`**: Archivo FASTA con el genoma de referencia.
3. **`CDS_ribosomal_protein_L2_cromosome1_S_phureja.fasta`**: Archivo FASTA con la secuencia del CDS a comparar.
4. **Resultados:**
   - Las gráficas se guardan en formato PDF en el directorio principal con nombres descriptivos.

---

## Uso

1. Coloca los archivos FASTA del genoma y el CDS en el directorio principal.
2. Ejecuta el script principal:
   ```bash
   python genome_sequence_matching.py
   ```
3. Las gráficas se generarán y guardarán en el directorio principal.
4. Consulta la consola para ver los valores  y otros cálculos relevantes.

---


## Licencia

Este proyecto está licenciado bajo la [MIT License](LICENSE).



