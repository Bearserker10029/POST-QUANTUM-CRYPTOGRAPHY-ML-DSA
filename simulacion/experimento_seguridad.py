# -*- coding: utf-8 -*-
"""
experimento_seguridad.py — Experimento 2: propiedades de seguridad de ML-DSA
==============================================================================
Batería de pruebas negativas y de consistencia sobre mldsa_puro.py:

  a) alteración de un solo bit de la firma (en c̃, en z y en h) -> rechazo
  b) mismo mensaje con firma de otro mensaje -> rechazo
  c) verificación con una clave pública ajena -> rechazo
  d) cadena de contexto distinta -> rechazo; contexto correcto -> aceptación
  e) modo determinista: reproducible bit a bit; modo hedged: firma distinta
     cada vez; ambas verifican
  f) tamaños exactos de claves y firma (FIPS 204, Tabla 2)

Uso:  python3 experimento_seguridad.py
"""

import mldsa_puro as mld

FALLOS = []


def comprobar(nombre, condicion):
    print(f"  [{'OK' if condicion else 'FALLO'}] {nombre}")
    if not condicion:
        FALLOS.append(nombre)


def correr(params):
    print(f"\n== {params} ==")
    pk, sk = mld.keygen(params)
    pk2, sk2 = mld.keygen(params)
    msg = b"Transferir 1000 USD a la cuenta 123-456"
    ctx = b"BANCA-2026"
    sig = mld.sign(sk, msg, ctx, params)

    # a) un solo bit alterado en cada componente de la firma
    lam4 = mld.PARAMS[params][7] // 4          # bytes de c̃
    pos_z = lam4 + 10                           # dentro de z
    pos_h = len(sig) - 5                        # dentro de h
    for nombre, pos in (("c̃", 3), ("z", pos_z), ("h", pos_h)):
        mala = bytearray(sig)
        mala[pos] ^= 0x01
        comprobar(f"bit alterado en {nombre} -> rechazada",
                  not mld.verify(pk, msg, bytes(mala), ctx, params))

    # b) firma de otro mensaje
    sig_otro = mld.sign(sk, b"Transferir 1 USD a la cuenta 123-456", ctx, params)
    comprobar("firma de otro mensaje -> rechazada",
              not mld.verify(pk, msg, sig_otro, ctx, params))

    # c) clave pública ajena
    comprobar("clave pública ajena -> rechazada",
              not mld.verify(pk2, msg, sig, ctx, params))

    # d) contexto
    comprobar("contexto distinto -> rechazada",
              not mld.verify(pk, msg, sig, b"BANCA-2025", params))
    comprobar("contexto correcto -> aceptada",
              mld.verify(pk, msg, sig, ctx, params))

    # e) determinista vs hedged
    s1 = mld.sign(sk, msg, ctx, params, rnd=b"\x00" * 32)
    s2 = mld.sign(sk, msg, ctx, params, rnd=b"\x00" * 32)
    comprobar("determinista reproducible bit a bit", s1 == s2)
    h1 = mld.sign(sk, msg, ctx, params)          # hedged (rnd aleatorio)
    h2 = mld.sign(sk, msg, ctx, params)
    comprobar("hedged produce firmas distintas", h1 != h2)
    comprobar("ambas verifican",
              mld.verify(pk, msg, s1, ctx, params)
              and mld.verify(pk, msg, h1, ctx, params))

    # f) tamaños exactos
    s = mld.SIZES[params]
    comprobar(f"tamaños exactos pk={s['pk']} sk={s['sk']} sig={s['sig']}",
              len(pk) == s["pk"] and len(sk) == s["sk"] and len(sig) == s["sig"])


if __name__ == "__main__":
    for p in mld.PARAMS:
        correr(p)
    print("\nRESULTADO:", "todo correcto" if not FALLOS else f"{len(FALLOS)} fallos")
    raise SystemExit(1 if FALLOS else 0)
