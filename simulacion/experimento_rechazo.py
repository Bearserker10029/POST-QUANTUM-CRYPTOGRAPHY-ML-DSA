# -*- coding: utf-8 -*-
"""
experimento_rechazo.py — Experimento 1: muestreo por rechazo en la firma ML-DSA
================================================================================
Mide empíricamente el número de iteraciones del bucle de firma de FIPS 204
(Algoritmo 7) y lo compara con las repeticiones esperadas del estándar
(Tabla 1 de FIPS 204: 4,25 / 5,1 / 3,85) y con la distribución geométrica
teórica de parámetro p = 1/media.

Uso:  python3 experimento_rechazo.py [n44 n65 n87]
Por defecto: 200 firmas para ML-DSA-44, 400 para ML-DSA-65 y 120 para ML-DSA-87
(en Python puro cada firma tarda decenas de ms; ajuste a su máquina).
La aleatoriedad (mensajes y rnd) se deriva de una semilla fija, por lo que
los resultados son reproducibles bit a bit en cualquier máquina.
"""

import random
import sys
from collections import Counter

import mldsa_puro as mld

ESPERADAS = {"ML-DSA-44": 4.25, "ML-DSA-65": 5.1, "ML-DSA-87": 3.85}
SEMILLA = 20261018          # fija ⇒ resultados reproducibles


def correr(nombre, n, azar):
    pk, sk = mld.keygen_internal(azar.randbytes(32), mld.PARAMS[nombre])
    stats = []
    for i in range(n):
        msg = azar.randbytes(32) + i.to_bytes(4, "big")
        mld.sign(sk, msg, b"EXPERIMENTO", nombre, rnd=azar.randbytes(32),
                 stats=stats)
    media = sum(stats) / len(stats)
    maximo = max(stats)
    frec = Counter(stats)
    print(f"{nombre}: n={len(stats)}  media={media:.2f} "
          f"(teórica {ESPERADAS[nombre]})  máx={maximo}")
    # Comparación con la geométrica teórica
    p = 1.0 / ESPERADAS[nombre]
    print("  iter | empírica | geométrica teórica")
    for k in range(1, 13):
        emp = frec.get(k, 0) / len(stats)
        geo = (1 - p) ** (k - 1) * p
        print(f"  {k:4d} |  {emp:7.3f} |  {geo:7.3f}")
    return stats


if __name__ == "__main__":
    ns = [int(x) for x in sys.argv[1:4]] or [200, 400, 120]
    azar = random.Random(SEMILLA)
    for nombre, n in zip(mld.PARAMS, ns):
        correr(nombre, n, azar)
