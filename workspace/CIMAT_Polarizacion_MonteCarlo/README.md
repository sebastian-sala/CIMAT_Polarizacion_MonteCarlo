# Modelado Monte Carlo de la polarización del cielo diurno

Repositorio de apoyo para el póster presentado en la Escuela de Verano CIMAT
Mérida 2026.

## Resumen

Se simuló mediante el método de Monte Carlo la propagación de fotones en una
atmósfera plano-paralela con dispersión Rayleigh-aerosoles, utilizando el
formalismo de Stokes-Mueller. El objetivo fue estudiar la aparición de una
componente circular débil, medida mediante la elipticidad neta \(V/I\), como
efecto estadístico asociado a dispersión múltiple.

## Contenido

```text
CIMAT_Polarizacion_MonteCarlo/
├── README.md
├── codigo/
│   └── simulacion_polarizacion_montecarlo.m
├── resultados/
│   ├── corridas_100.csv
│   └── resumen_estadistico.txt
├── figuras/
│   ├── histograma_colisiones.png
│   ├── esfera_poincare.png
│   └── elipse_polarizacion.png
├── referencias/
│   └── referencias.bib
├── poster/
│   └── poster_cimat_montecarlo.pdf
└── reporte/
    └── reporte_extendido.pdf
```

## Resultados principales

Para 100 corridas independientes de \(10\,000\) fotones:

- Promedio firmado: \(\overline{V/I}=2.20\times10^{-3}\).
- Desviación estándar entre corridas: \(\sigma_{V/I}=1.60\times10^{-2}\).
- Magnitud media: \(\overline{|V/I|}=1.01\times10^{-2}\).
- Mediana de la magnitud: \(\mathrm{mediana}(|V/I|)=6.68\times10^{-3}\).

El resultado robusto no es el signo de una corrida individual, sino la escala
de la componente circular: una elipticidad débil del orden de \(10^{-3}\) a
\(10^{-2}\), frente a una polarización lineal dominante.

## Cómo correr el código

1. Abrir `codigo/simulacion_polarizacion_montecarlo.m` en MATLAB.
2. Ejecutar el script.
3. El programa simula \(10\,000\) fotones, calcula el vector de Stokes promedio
   de los fotones transmitidos con dispersión múltiple y genera tres figuras:
   histograma de colisiones, esfera de Poincaré y elipse de polarización.

## Autores

- Alonso Mendoza Hernández
- Jafet Alani Rivera Torres
- Sebastián Sala Baltazar

Facultad de Ciencias, Universidad Nacional Autónoma de México.
