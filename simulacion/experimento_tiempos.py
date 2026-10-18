# -*- coding: utf-8 -*-
"""
experimento_tiempos.py — Experimento 3: rendimiento de la implementación educativa
====================================================================================
Mide tiempos de generación de claves, firma y verificación de mldsa_puro.py
(Python puro, sin optimización vectorial). Estas cifras NO son representativas
de implementaciones de producción (C/Rust/AVX2), que son entre 100 y 500
veces más rápidas; sirven para ilustrar el coste relativo de cada operación
y la variabilidad de la firma por el bucle de rechazo.

Uso:  python3 experimento_tiempos.py [repeticiones]
"""

import statistics
import sys
import time

import mldsa_puro as mld


def medir(func, n):
    ts = []
    for _ in range(n):
        t0 = time.perf_counter()
        func()
        ts.append((time.perf_counter() - t0) * 1000)   # ms
    return ts


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    msg = b"Benchmark ML-DSA en Python puro"
    print(f"repeticiones: {n} por operación y nivel\n")
    for params in mld.PARAMS:
        pk, sk = mld.keygen(params)
        t_kg = medir(lambda: mld.keygen(params), n)
        t_sg = medir(lambda: mld.sign(sk, msg, b"BENCH", params), n)
        sig = mld.sign(sk, msg, b"BENCH", params)
        t_vf = medir(lambda: mld.verify(pk, msg, sig, b"BENCH", params), n)
        cv = statistics.stdev(t_sg) / statistics.mean(t_sg) * 100
        print(f"{params}: keygen={statistics.mean(t_kg):.1f} ms  "
              f"firma={statistics.mean(t_sg):.1f} ms (CV {cv:.0f}%)  "
              f"verif={statistics.mean(t_vf):.1f} ms")
