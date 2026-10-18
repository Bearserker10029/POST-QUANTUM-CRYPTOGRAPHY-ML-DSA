# -*- coding: utf-8 -*-
"""
validar_acvp.py — Validación de mldsa_puro.py contra los vectores ACVP del NIST
=================================================================================
Descarga (si no existen) los ficheros oficiales de generación y validación del
repositorio ACVP-Server del NIST y ejecuta todos los casos de prueba de
ML-DSA keyGen, sigGen y sigVer (FIPS 204), comparando bit a bit.

Uso:
    python3 validar_acvp.py [directorio_con_los_json]

Si no se indica directorio, usa ./acvp_vectores (los descarga si faltan).
Se registran la fecha de descarga y el total de casos ejecutados.

Referencia: https://github.com/usnistgov/ACVP-Server
            gen-val/json-files/ML-DSA-*-FIPS204/
"""

import json
import os
import sys
import time
import urllib.request

import mldsa_puro as mld

BASE = ("https://raw.githubusercontent.com/usnistgov/ACVP-Server/master/"
        "gen-val/json-files")
ALGS = ["ML-DSA-keyGen-FIPS204", "ML-DSA-sigGen-FIPS204", "ML-DSA-sigVer-FIPS204"]


def descargar(directorio):
    os.makedirs(directorio, exist_ok=True)
    for alg in ALGS:
        for nombre, remoto in (("prompt", "internalProjection.json"),
                               ("expected", "expectedResults.json")):
            ruta = os.path.join(directorio, f"{alg}_{nombre}.json")
            if os.path.exists(ruta):
                continue
            print(f"descargando {alg}/{remoto} ...", flush=True)
            urllib.request.urlretrieve(f"{BASE}/{alg}/{remoto}", ruta)


def cargar(directorio, alg, nombre):
    with open(os.path.join(directorio, f"{alg}_{nombre}.json")) as f:
        return json.load(f)


def _m_prime_externa(tc, ctx_hex):
    """M′ = 0x00 ‖ len(ctx) ‖ ctx ‖ M (interfaz externa pura)."""
    ctx = bytes.fromhex(ctx_hex) if ctx_hex else b""
    return b"\x00" + bytes([len(ctx)]) + ctx + bytes.fromhex(tc["message"])


def _m_prime_prehash(tc, ctx_hex, hash_alg):
    """M′ = 0x01 ‖ len(ctx) ‖ ctx ‖ OID ‖ PH(M) (HashML-DSA)."""
    ctx = bytes.fromhex(ctx_hex) if ctx_hex else b""
    return (b"\x01" + bytes([len(ctx)]) + ctx + mld.PH_OID[hash_alg]
            + mld._HASHLIB[hash_alg](bytes.fromhex(tc["message"])))


def validar_keygen(directorio):
    prompt = cargar(directorio, ALGS[0], "prompt")
    esperado = cargar(directorio, ALGS[0], "expected")
    total = ok = 0
    for gp, ge in zip(prompt["testGroups"], esperado["testGroups"]):
        params = mld.PARAMS[gp["parameterSet"]]
        for tp, te in zip(gp["tests"], ge["tests"]):
            total += 1
            pk, sk = mld.keygen_internal(bytes.fromhex(tp["seed"]), params)
            if pk.hex().upper() == te["pk"].upper() and \
                    sk.hex().upper() == te["sk"].upper():
                ok += 1
            else:
                print(f"  FALLO keyGen tgId={gp['tgId']} tcId={tp['tcId']}")
    return total, ok


def _es_prehash(gp):
    """preHash viene como cadena 'pure'/'preHash' (o booleano en borradores)."""
    v = gp.get("preHash", False)
    return v is True or v == "preHash"


def _hash_alg(gp, tp):
    return tp.get("hashAlg") or gp.get("hashAlg")


def _contexto(gp, tp):
    return tp.get("context") or gp.get("context") or ""


def _firma_caso(gp, tp):
    """Ejecuta un caso sigGen y devuelve la firma generada."""
    params = mld.PARAMS[gp["parameterSet"]]
    sk = bytes.fromhex(tp["sk"])
    determinista = gp.get("deterministic", False)
    rnd = b"\x00" * 32 if determinista else bytes.fromhex(tp["rnd"])
    interfaz = gp.get("signatureInterface", "external")
    prehash = _es_prehash(gp)
    external_mu = gp.get("externalMu", False) or "mu" in tp

    if external_mu:                                # interfaz ExternalMu
        return mld.sign_internal(sk, b"", rnd, params,
                                 mu=bytes.fromhex(tp["mu"]))
    if interfaz == "internal":                     # M′ entra tal cual
        return mld.sign_internal(sk, bytes.fromhex(tp["message"]), rnd, params)
    if prehash:                                    # HashML-DSA (externa)
        m_prime = _m_prime_prehash(tp, _contexto(gp, tp), _hash_alg(gp, tp))
    else:                                          # ML-DSA pura (externa)
        m_prime = _m_prime_externa(tp, _contexto(gp, tp))
    return mld.sign_internal(sk, m_prime, rnd, params)


def validar_siggen(directorio):
    prompt = cargar(directorio, ALGS[1], "prompt")
    esperado = cargar(directorio, ALGS[1], "expected")
    total = ok = 0
    desglose = {}
    for gp, ge in zip(prompt["testGroups"], esperado["testGroups"]):
        clave = (gp["parameterSet"], gp.get("signatureInterface", "external"),
                 gp.get("preHash", "none"),
                 gp.get("externalMu", False) or "mu" in gp["tests"][0],
                 gp.get("deterministic", False))
        for tp, te in zip(gp["tests"], ge["tests"]):
            total += 1
            sig = _firma_caso(gp, tp)
            if sig.hex().upper() == te["signature"].upper():
                ok += 1
                desglose[clave] = desglose.get(clave, 0) + 1
            else:
                print(f"  FALLO sigGen tgId={gp['tgId']} tcId={tp['tcId']}")
    return total, ok, desglose


def validar_sigver(directorio):
    prompt = cargar(directorio, ALGS[2], "prompt")
    esperado = cargar(directorio, ALGS[2], "expected")
    total = ok = 0
    for gp, ge in zip(prompt["testGroups"], esperado["testGroups"]):
        params = mld.PARAMS[gp["parameterSet"]]
        interfaz = gp.get("signatureInterface", "external")
        prehash = _es_prehash(gp)
        external_mu = gp.get("externalMu", False)
        for tp, te in zip(gp["tests"], ge["tests"]):
            total += 1
            pk = bytes.fromhex(tp["pk"])
            sig = bytes.fromhex(tp["signature"])
            if external_mu:
                res = mld.verify_internal(pk, b"", sig, params,
                                          mu=bytes.fromhex(tp["mu"]))
            elif interfaz == "internal":
                res = mld.verify_internal(pk, bytes.fromhex(tp["message"]),
                                          sig, params)
            elif prehash:
                res = mld.verify_internal(
                    pk, _m_prime_prehash(tp, _contexto(gp, tp), _hash_alg(gp, tp)),
                    sig, params)
            else:
                res = mld.verify_internal(
                    pk, _m_prime_externa(tp, _contexto(gp, tp)), sig, params)
            if res == te["testPassed"]:
                ok += 1
            else:
                print(f"  FALLO sigVer tgId={gp['tgId']} tcId={tp['tcId']}")
    return total, ok


def main():
    directorio = sys.argv[1] if len(sys.argv) > 1 else "acvp_vectores"
    descargar(directorio)
    t0 = time.time()

    t_kg, ok_kg = validar_keygen(directorio)
    print(f"keyGen : {ok_kg}/{t_kg}")
    t_sg, ok_sg, desglose = validar_siggen(directorio)
    print(f"sigGen : {ok_sg}/{t_sg}")
    for clave, n in sorted(desglose.items()):
        print(f"         {clave}: {n}")
    t_sv, ok_sv = validar_sigver(directorio)
    print(f"sigVer : {ok_sv}/{t_sv}")

    total, oks = t_kg + t_sg + t_sv, ok_kg + ok_sg + ok_sv
    print(f"TOTAL  : {oks}/{total} casos en {time.time() - t0:.1f} s")
    sys.exit(0 if oks == total else 1)


if __name__ == "__main__":
    main()
