# Simulación ML-DSA (FIPS 204) en Python puro

Material ejecutable del informe «Post-Quantum Cryptography: ML-DSA».
Implementación educativa, sin dependencias externas (solo la biblioteca
estándar de Python 3.9+; `hashlib` con SHAKE-128/256), de los tres conjuntos
de parámetros de FIPS 204: ML-DSA-44, ML-DSA-65 y ML-DSA-87.

## Contenido

| Fichero | Descripción |
|---|---|
| `mldsa_puro.py` | Implementación completa comentada en español: NTT, muestreo por rechazo, `ExpandA`, `SampleInBall`, `Power2Round`, `Decompose`, `MakeHint`/`UseHint`, Algoritmos 6-8 (KeyGen/Sign/Verify internos), API externa con `ctx`, modos **hedged** y **determinista**, **HashML-DSA** (pre-hash con OID DER de FIPS 204 §5.4) y **ExternalMu** (μ externo, Algoritmo 7 línea 6 / Algoritmo 8 línea 7). Incluye una demo en `__main__`. |
| `validar_acvp.py` | Descarga (si faltan) los 6 JSON oficiales del repositorio ACVP-Server del NIST y valida **bit a bit los 615 casos**: keyGen 75, sigGen 360 (24 grupos × 15), sigVer 180 (12 grupos × 15). Cubre interfaces externa pura, externa pre-hash e interna, ExternalMu, y modos determinista/hedged. |
| `experimento_rechazo.py` | Experimento 1: número de iteraciones del bucle de firma frente a la media teórica (4,25 / 5,1 / 3,85) y a la distribución geométrica. Semilla fija: resultados reproducibles bit a bit. |
| `experimento_seguridad.py` | Experimento 2: 33 pruebas negativas y de consistencia (bit alterado en c̃/z/h, mensaje o clave pública distintos, `ctx` incorrecto, determinista reproducible, hedged variable, tamaños exactos). |
| `experimento_tiempos.py` | Experimento 3: tiempos de keygen/firma/verificación (medias de 25 repeticiones) y coeficiente de variación. |
| `fig1_tamanos.png` … `fig5_flujo.png` | Figuras del informe (tamaños, rechazo, benchmarks, cronología, flujo de algoritmos). |

## Uso rápido

```bash
python3 mldsa_puro.py              # demo: firma pura, HashML-DSA y ExternalMu
python3 validar_acvp.py            # valida los 615 casos ACVP (descarga los JSON)
python3 experimento_rechazo.py     # experimento 1 (≈1 min)
python3 experimento_seguridad.py   # experimento 2 (≈30 s)
python3 experimento_tiempos.py     # experimento 3 (≈2 min)
```

## Resultados obtenidos (18 de octubre de 2026)

- **Validación ACVP** (vectores del repositorio `usnistgov/ACVP-Server`,
  `gen-val/json-files/ML-DSA-*-FIPS204`, descargados el 18-10-2026):
  **615/615 casos** — keyGen 75/75, sigGen 360/360, sigVer 180/180 —
  en ≈26 s en un núcleo de contenedor.
- **Rechazo**: medias empíricas 4,50 / 5,34 / 3,44 iteraciones
  (n = 200/400/120) frente a las teóricas 4,25 / 5,1 / 3,85; máximos 27/27/14.
- **Seguridad**: 33/33 pruebas superadas.
- **Tiempos** (Python puro, 1 núcleo): ML-DSA-44 → 7,6 / 31,4 / 8,5 ms;
  ML-DSA-65 → 12,4 / 68,4 / 12,2 ms; ML-DSA-87 → 17,4 / 73,3 / 18,3 ms
  (keygen / firma / verificación; CV de la firma 53-60 %).

> Nota: esta es una implementación **didáctica y de validación**, no de
> producción: no garantiza tiempo constante ni borrado seguro de secretos.
> Para uso real, utilice una biblioteca certificada (p. ej., OpenSSL 3.5+,
> liboqs, Bouncy Castle).
