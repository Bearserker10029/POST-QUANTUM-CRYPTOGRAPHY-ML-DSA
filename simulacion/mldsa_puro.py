# -*- coding: utf-8 -*-
"""
ML-DSA (FIPS 204) — Implementación educativa en Python puro
============================================================
Implementación de referencia didáctica del estándar NIST FIPS 204
(Module-Lattice-Based Digital Signature Standard, agosto 2024),
derivado de CRYSTALS-Dilithium.

Incluye los tres conjuntos de parámetros (ML-DSA-44, ML-DSA-65 y ML-DSA-87),
la API externa con cadena de contexto (Algoritmos 1-3), la variante pre-hash
HashML-DSA (Algoritmos 4-5, con la tabla de OIDs de la Sección 5.4), los
modos 'hedged' y determinista, y la interfaz ExternalMu (representante del
mensaje μ calculado fuera del módulo firmante).

Validada bit a bit con los vectores oficiales ACVP del NIST; véase
validar_acvp.py para reproducir la comprobación.

ADVERTENCIA: código con fines académicos/demostrativos. No es de tiempo
constante ni resistente a canales laterales. Para producción use una
biblioteca validada (p. ej., OpenSSL 3.5+, liboqs, mldsa-native).
"""

import hashlib
import os

# ----------------------------------------------------------------------------
# Parámetros globales (FIPS 204, Sección 2-4)
# ----------------------------------------------------------------------------
Q = 8380417          # primo q = 2^23 - 2^13 + 1
N = 256              # grado del anillo R_q = Z_q[X]/(X^256 + 1)
ZETA = 1753          # raíz 512-ésima primitiva de la unidad mod q
D = 13               # bits descartados de t (Power2Round)

PARAMS = {
    # nombre:      (k, l, eta, tau, gamma1,      gamma2,        omega, lam)
    "ML-DSA-44": (4, 4, 2, 39, 2**17, (Q - 1) // 88, 80, 128),
    "ML-DSA-65": (6, 5, 4, 49, 2**19, (Q - 1) // 32, 55, 192),
    "ML-DSA-87": (8, 7, 2, 60, 2**19, (Q - 1) // 32, 75, 256),
}

# ----------------------------------------------------------------------------
# Utilidades de bits/bytes (FIPS 204, Sección 7.1)
# ----------------------------------------------------------------------------

def bits_to_bytes(b):
    """Algorithm 12: BitsToBytes."""
    assert len(b) % 8 == 0
    out = bytearray(len(b) // 8)
    for i, bit in enumerate(b):
        out[i // 8] |= (bit & 1) << (i % 8)
    return bytes(out)


def bytes_to_bits(z):
    """Algorithm 13: BytesToBits."""
    b = []
    for byte in z:
        for i in range(8):
            b.append((byte >> i) & 1)
    return b


def integer_to_bits(x, a):
    """Algorithm 9: IntegerToBits (little-endian, a bits)."""
    return [(x >> i) & 1 for i in range(a)]


def integer_to_bytes(x, a):
    """Algorithm 11: IntegerToBytes (little-endian, a bytes)."""
    return bytes((x >> (8 * i)) & 0xFF for i in range(a))


def bitlen(a):
    return a.bit_length()


def shake128(data, outlen):
    return hashlib.shake_128(data).digest(outlen)


def shake256(data, outlen):
    return hashlib.shake_256(data).digest(outlen)


# ----------------------------------------------------------------------------
# Empaquetado de coeficientes (FIPS 204, Sección 7.2)
# ----------------------------------------------------------------------------

def simple_bit_pack(w, b):
    """Algorithm 16: empaqueta coeficientes de [0, b] con bitlen(b) bits."""
    bl = bitlen(b)
    bits = []
    for c in w:
        bits.extend(integer_to_bits(c, bl))
    return bits_to_bytes(bits)


def bit_pack(w, a, b):
    """Algorithm 17: empaqueta w en [-a, b] como b - w
    (FIPS 204: 'BitPack works by merely subtracting w from b')."""
    bl = bitlen(a + b)
    bits = []
    for c in w:
        bits.extend(integer_to_bits(b - c, bl))
    return bits_to_bytes(bits)


def simple_bit_unpack(v, b):
    """Algorithm 18."""
    bl = bitlen(b)
    bits = bytes_to_bits(v)
    w = []
    for i in range(N):
        c = 0
        for j in range(bl):
            c |= bits[i * bl + j] << j
        w.append(c)
    return w


def bit_unpack(v, a, b):
    """Algorithm 19."""
    bl = bitlen(a + b)
    bits = bytes_to_bits(v)
    w = []
    for i in range(N):
        c = 0
        for j in range(bl):
            c |= bits[i * bl + j] << j
        w.append(b - c)   # FIPS 204 Alg. 19: w_i <- b - BitsToInteger(...)
    return w


def hint_bit_pack(h, omega, k):
    """Algorithm 20: empaqueta el vector de pistas h en omega+k bytes."""
    y = bytearray(omega + k)
    idx = 0
    for i in range(k):
        for j in range(N):
            if h[i][j] != 0:
                y[idx] = j
                idx += 1
        y[omega + i] = idx
    return bytes(y)


def hint_bit_unpack(y, omega, k):
    """Algorithm 21: devuelve h o None si está malformado."""
    if len(y) != omega + k:
        return None
    h = [[0] * N for _ in range(k)]
    idx = 0
    for i in range(k):
        if y[omega + i] < idx or y[omega + i] > omega:
            return None
        first = idx
        while idx < y[omega + i]:
            j = y[idx]
            if idx > first and y[idx - 1] >= j:
                return None          # los índices deben ser crecientes
            h[i][j] = 1
            idx += 1
    for i in range(idx, omega):
        if y[i] != 0:
            return None              # bytes sobrantes deben ser cero
    return h


# ----------------------------------------------------------------------------
# Aritmética modular y NTT (FIPS 204, Secciones 7.5-7.6 y Apéndice B)
# ----------------------------------------------------------------------------

def mod_pm(r, m):
    """Reducción centrada r mod^± m, resultado en (-m/2, m/2]."""
    s = r % m
    if s > m // 2:
        s -= m
    return s


def _bitrev8(x):
    """Algorithm 44: BitRev8."""
    r = 0
    for i in range(8):
        r |= ((x >> i) & 1) << (7 - i)
    return r


ZETAS = [pow(ZETA, _bitrev8(k), Q) for k in range(256)]


def ntt(w):
    """Algorithm 41: transformada NTT, w en R_q -> T_q."""
    w = [x % Q for x in w]
    k = 1
    m = 128
    while m >= 1:
        for start in range(0, 256, 2 * m):
            z = ZETAS[k]
            k += 1
            for j in range(start, start + m):
                t = (z * w[j + m]) % Q
                w[j + m] = (w[j] - t) % Q
                w[j] = (w[j] + t) % Q
        m //= 2
    return w


def ntt_inv(w_hat):
    """Algorithm 42: NTT inversa, T_q -> R_q."""
    w = [x % Q for x in w_hat]
    k = 256
    m = 1
    while m < 256:
        for start in range(0, 256, 2 * m):
            k -= 1
            z = ZETAS[k]
            for j in range(start, start + m):
                t = w[j]
                w[j] = (t + w[j + m]) % Q
                w[j + m] = (z * (w[j + m] - t)) % Q
        m *= 2
    inv256 = pow(256, -1, Q)
    return [(x * inv256) % Q for x in w]


def poly_add(a, b):
    return [(x + y) % Q for x, y in zip(a, b)]


def poly_sub(a, b):
    return [(x - y) % Q for x, y in zip(a, b)]


def poly_mul_ntt(a_hat, b_hat):
    """Algorithm 45: MultiplyNTT (multiplicación punto a punto en T_q)."""
    return [(x * y) % Q for x, y in zip(a_hat, b_hat)]


# ----------------------------------------------------------------------------
# Muestreo pseudoaleatorio (FIPS 204, Sección 7.3)
# ----------------------------------------------------------------------------

def coeff_from_three_bytes(b0, b1, b2):
    """Algorithm 14."""
    z = ((b2 & 0x7F) << 16) | (b1 << 8) | b0
    return z if z < Q else None


def coeff_from_half_byte(b, eta):
    """Algorithm 15."""
    if eta == 2:
        if b < 15:
            return 2 - (b % 5)
        return None
    if b < 9:
        return 4 - b
    return None


def rej_ntt_poly(seed34):
    """Algorithm 30: RejNTTPoly — polinomio uniforme en T_q desde SHAKE128."""
    coeffs = []
    stream = shake128(seed34, 840)          # 840 bytes >= necesidad esperada
    pos = 0
    while len(coeffs) < N:
        if pos + 3 > len(stream):
            stream += shake128(seed34, len(stream) + 840)[len(stream):]
        z = coeff_from_three_bytes(stream[pos], stream[pos + 1], stream[pos + 2])
        pos += 3
        if z is not None:
            coeffs.append(z)
    return coeffs


def rej_bounded_poly(seed66, eta):
    """Algorithm 31: RejBoundedPoly — coeficientes pequeños en [-eta, eta]."""
    coeffs = []
    stream = shake256(seed66, 520)
    pos = 0
    while len(coeffs) < N:
        if pos >= len(stream):
            stream += shake256(seed66, len(stream) + 520)[len(stream):]
        byte = stream[pos]
        pos += 1
        for b in (byte & 0x0F, byte >> 4):
            z = coeff_from_half_byte(b, eta)
            if z is not None and len(coeffs) < N:
                coeffs.append(z)
    return coeffs


def expand_a(rho, k, l):
    """Algorithm 32: ExpandA — matriz Â k×l en dominio NTT."""
    A = []
    for r in range(k):
        row = []
        for s in range(l):
            row.append(rej_ntt_poly(rho + integer_to_bytes(s, 1) + integer_to_bytes(r, 1)))
        A.append(row)
    return A


def expand_s(rho2, eta, k, l):
    """Algorithm 33: ExpandS — vectores secretos s1 (l) y s2 (k)."""
    s1 = [rej_bounded_poly(rho2 + integer_to_bytes(r, 2), eta) for r in range(l)]
    s2 = [rej_bounded_poly(rho2 + integer_to_bytes(l + r, 2), eta) for r in range(k)]
    return s1, s2


def expand_mask(rho2, kappa, gamma1, l):
    """Algorithm 34: ExpandMask — vector de máscaras y."""
    n_poly_bytes = 32 * bitlen(2 * gamma1 - 1)   # 576 B (2^17) o 640 B (2^19)
    y = []
    for r in range(l):
        v = shake256(rho2 + integer_to_bytes(kappa + r, 2), n_poly_bytes)
        y.append(bit_unpack(v, gamma1 - 1, gamma1))
    return y


def sample_in_ball(seed32, tau):
    """Algorithm 29: SampleInBall — polinomio c con tau coeficientes ±1."""
    c = [0] * N
    stream = shake256(seed32, 8 + 136 + 64)   # 8 bytes de signos + margen
    signs = stream[:8]
    pos = 8
    for i in range(N - tau, N):
        while True:
            if pos >= len(stream):
                stream += shake256(seed32, len(stream) + 136)[len(stream):]
            j = stream[pos]
            pos += 1
            if j <= i:
                break
        c[i] = c[j]
        sign_bit = (signs[(i - (N - tau)) // 8] >> ((i - (N - tau)) % 8)) & 1
        c[j] = -1 if sign_bit else 1
    return c


# ----------------------------------------------------------------------------
# Redondeo y pistas (FIPS 204, Sección 7.4)
# ----------------------------------------------------------------------------

def power2round(r):
    """Algorithm 35: r = r1·2^d + r0 con r0 en (-2^(d-1), 2^(d-1)]."""
    r = r % Q
    r0 = mod_pm(r, 2**D)
    return (r - r0) >> D, r0


def decompose(r, gamma2):
    """Algorithm 36."""
    r = r % Q
    r0 = mod_pm(r, 2 * gamma2)
    if r - r0 == Q - 1:
        return 0, r0 - 1
    return (r - r0) // (2 * gamma2), r0


def high_bits(r, gamma2):
    return decompose(r, gamma2)[0]


def low_bits(r, gamma2):
    return decompose(r, gamma2)[1]


def make_hint(z, r, gamma2):
    """Algorithm 39: 1 si HighBits(r) != HighBits(r+z)."""
    return 1 if high_bits(r, gamma2) != high_bits(r + z, gamma2) else 0


def use_hint(h, r, gamma2):
    """Algorithm 41 (MakeHint/UseHint): ajusta los bits altos con la pista."""
    m = (Q - 1) // (2 * gamma2)
    r1, r0 = decompose(r, gamma2)
    if h == 1:
        if r0 > 0:
            return (r1 + 1) % m
        return (r1 - 1) % m
    return r1


# ----------------------------------------------------------------------------
# Codificación de claves y firmas (FIPS 204, Sección 7.2)
# ----------------------------------------------------------------------------

def pk_encode(rho, t1):
    """Algorithm 22."""
    out = bytearray(rho)
    for poly in t1:
        out += simple_bit_pack([c % Q for c in poly], 2**(bitlen(Q - 1) - D) - 1)
    return bytes(out)


def pk_decode(pk, k):
    """Algorithm 23."""
    bl = bitlen(Q - 1) - D
    nbytes = N * bl // 8   # 256 coeficientes x 10 bits = 320 bytes
    rho = pk[:32]
    t1 = [simple_bit_unpack(pk[32 + i * nbytes: 32 + (i + 1) * nbytes],
                            2**bl - 1) for i in range(k)]
    return rho, t1


def sk_encode(rho, K, tr, s1, s2, t0, eta):
    """Algorithm 24."""
    out = bytearray(rho + K + tr)
    for poly in s1:
        out += bit_pack(poly, eta, eta)
    for poly in s2:
        out += bit_pack(poly, eta, eta)
    for poly in t0:
        out += bit_pack(poly, 2**(D - 1) - 1, 2**(D - 1))
    return bytes(out)


def sk_decode(sk, k, l, eta):
    """Algorithm 25."""
    rho = sk[0:32]
    K = sk[32:64]
    tr = sk[64:128]
    pos = 128
    nb_eta = 32 * bitlen(2 * eta)
    s1 = [bit_unpack(sk[pos + i * nb_eta: pos + (i + 1) * nb_eta], eta, eta)
          for i in range(l)]
    pos += l * nb_eta
    s2 = [bit_unpack(sk[pos + i * nb_eta: pos + (i + 1) * nb_eta], eta, eta)
          for i in range(k)]
    pos += k * nb_eta
    nb_t0 = 32 * D
    t0 = [bit_unpack(sk[pos + i * nb_t0: pos + (i + 1) * nb_t0],
                     2**(D - 1) - 1, 2**(D - 1)) for i in range(k)]
    return rho, K, tr, s1, s2, t0


def sig_encode(c_tilde, z, h, gamma1, omega, k):
    """Algorithm 26."""
    out = bytearray(c_tilde)
    for poly in z:
        out += bit_pack([mod_pm(c, Q) for c in poly], gamma1 - 1, gamma1)
    out += hint_bit_pack(h, omega, k)
    return bytes(out)


def sig_decode(sig, gamma1, omega, k, l, lam):
    """Algorithm 27."""
    cbytes = lam // 4
    zbytes = 32 * bitlen(2 * gamma1 - 1)
    expected = cbytes + l * zbytes + omega + k
    if len(sig) != expected:
        return None
    c_tilde = sig[:cbytes]
    z = [bit_unpack(sig[cbytes + i * zbytes: cbytes + (i + 1) * zbytes],
                    gamma1 - 1, gamma1) for i in range(l)]
    h = hint_bit_unpack(sig[cbytes + l * zbytes:], omega, k)
    if h is None:
        return None
    return c_tilde, z, h


def w1_encode(w1, gamma2, k):
    """Algorithm 28."""
    max_w1 = (Q - 1) // (2 * gamma2) - 1
    out = bytearray()
    for poly in w1:
        out += simple_bit_pack(poly, max_w1)
    return bytes(out)


# ----------------------------------------------------------------------------
# Núcleo matemático: generación, firma y verificación internas
# ----------------------------------------------------------------------------

def _mat_vec_ntt(A, v_hat):
    """Multiplicación matriz-vector en dominio NTT."""
    k = len(A)
    out = []
    for r in range(k):
        acc = [0] * N
        for s in range(len(v_hat)):
            prod = poly_mul_ntt(A[r][s], v_hat[s])
            acc = poly_add(acc, prod)
        out.append(acc)
    return out


def keygen_internal(xi, params):
    """Algorithm 6: ML-DSA.KeyGen_internal(ξ)."""
    k, l, eta, tau, gamma1, gamma2, omega, lam = params
    seed = shake256(xi + integer_to_bytes(k, 1) + integer_to_bytes(l, 1), 128)
    rho, rho2, K = seed[:32], seed[32:96], seed[96:128]

    A = expand_a(rho, k, l)
    s1, s2 = expand_s(rho2, eta, k, l)

    s1_hat = [ntt(p) for p in s1]
    t = [poly_add(ntt_inv(mv), e) for mv, e in zip(_mat_vec_ntt(A, s1_hat), s2)]

    t1, t0 = [], []
    for poly in t:
        p1, p0 = [], []
        for c in poly:
            r1, r0 = power2round(c)
            p1.append(r1)
            p0.append(r0)
        t1.append(p1)
        t0.append(p0)

    pk = pk_encode(rho, t1)
    tr = shake256(pk, 64)
    sk = sk_encode(rho, K, tr, s1, s2, t0, eta)
    return pk, sk


def sign_internal(sk, M_prime, rnd, params, stats=None, mu=None):
    """Algorithm 7: ML-DSA.Sign_internal(sk, M′, rnd).

    `stats` (opcional): lista donde se anota el número de iteraciones del
    bucle de muestreo por rechazo (para fines demostrativos).
    `mu` (opcional): representante del mensaje ya calculado (64 B); si se
    proporciona, se omite su cómputo — es la interfaz "ExternalMu", en la
    que μ puede calcularse en un módulo criptográfico distinto (FIPS 204,
    nota de la línea 6 del Algoritmo 7).
    """
    k, l, eta, tau, gamma1, gamma2, omega, lam = params
    beta = tau * eta
    rho, K, tr, s1, s2, t0 = sk_decode(sk, k, l, eta)

    s1_hat = [ntt([c % Q for c in p]) for p in s1]
    s2_hat = [ntt([c % Q for c in p]) for p in s2]
    t0_hat = [ntt([c % Q for c in p]) for p in t0]
    A = expand_a(rho, k, l)

    if mu is None:
        mu = shake256(tr + M_prime, 64)
    rho2 = shake256(K + rnd + mu, 64)

    kappa = 0
    while True:
        y = expand_mask(rho2, kappa, gamma1, l)
        kappa += l

        y_hat = [ntt([c % Q for c in p]) for p in y]
        w = [ntt_inv(mv) for mv in _mat_vec_ntt(A, y_hat)]
        w1 = [[high_bits(c, gamma2) for c in p] for p in w]

        c_tilde = shake256(mu + w1_encode(w1, gamma2, k), lam // 4)
        c = sample_in_ball(c_tilde, tau)
        c_hat = ntt([x % Q for x in c])

        cs1 = [ntt_inv(poly_mul_ntt(c_hat, s)) for s in s1_hat]
        z = [poly_add(p, cp) for p, cp in zip(y, cs1)]
        if max(max(abs(mod_pm(x, Q)) for x in p) for p in z) >= gamma1 - beta:
            continue                        # rechazo: z filtraría s1

        cs2 = [ntt_inv(poly_mul_ntt(c_hat, s)) for s in s2_hat]
        w_minus_cs2 = [poly_sub(p, cp) for p, cp in zip(w, cs2)]
        r0 = [[low_bits(x, gamma2) for x in p] for p in w_minus_cs2]
        if max(max(abs(x) for x in p) for p in r0) >= gamma2 - beta:
            continue                        # rechazo: pista inconsistente

        ct0 = [ntt_inv(poly_mul_ntt(c_hat, t)) for t in t0_hat]
        if max(max(abs(mod_pm(x, Q)) for x in p) for p in ct0) >= gamma2:
            continue                        # rechazo: ct0 demasiado grande

        h = []
        n_ones = 0
        for i in range(k):
            hp = []
            for j in range(N):
                hint = make_hint(-ct0[i][j] % Q, (w_minus_cs2[i][j] + ct0[i][j]) % Q,
                                 gamma2)
                hp.append(hint)
                n_ones += hint
            h.append(hp)
        if n_ones > omega:
            continue                        # rechazo: demasiadas pistas

        if stats is not None:
            stats.append(kappa // l)
        return sig_encode(c_tilde, z, h, gamma1, omega, k)


def verify_internal(pk, M_prime, sig, params, mu=None):
    """Algorithm 8: ML-DSA.Verify_internal(pk, M′, σ).

    `mu` (opcional): representante del mensaje ya calculado (64 B); si se
    proporciona, se omite su cómputo (interfaz "ExternalMu", nota de la
    línea 7 del Algoritmo 8).
    """
    k, l, eta, tau, gamma1, gamma2, omega, lam = params
    beta = tau * eta

    decoded = sig_decode(sig, gamma1, omega, k, l, lam)
    if decoded is None:
        return False
    c_tilde, z, h = decoded
    if max(max(abs(mod_pm(x, Q)) for x in p) for p in z) >= gamma1 - beta:
        return False

    rho, t1 = pk_decode(pk, k)
    A = expand_a(rho, k, l)
    if mu is None:
        tr = shake256(pk, 64)
        mu = shake256(tr + M_prime, 64)

    c = sample_in_ball(c_tilde, tau)
    c_hat = ntt([x % Q for x in c])

    z_hat = [ntt([mod_pm(x, Q) % Q for x in p]) for p in z]
    t1_hat = [ntt([(x << D) % Q for x in p]) for p in t1]

    Az = _mat_vec_ntt(A, z_hat)
    ct1 = [poly_mul_ntt(c_hat, t) for t in t1_hat]
    w_approx = [ntt_inv(poly_sub(a, b)) for a, b in zip(Az, ct1)]

    w1_prime = [[use_hint(h[i][j], w_approx[i][j], gamma2)
                 for j in range(N)] for i in range(k)]

    c_tilde2 = shake256(mu + w1_encode(w1_prime, gamma2, k), lam // 4)
    return c_tilde == c_tilde2


# ----------------------------------------------------------------------------
# API externa FIPS 204 (Algoritmos 1-3)
# ----------------------------------------------------------------------------

def keygen(params_name="ML-DSA-44"):
    """Algorithm 1: genera un par de claves con semilla aleatoria de 32 B."""
    xi = os.urandom(32)
    return keygen_internal(xi, PARAMS[params_name])


def _format_message(M, ctx):
    if len(ctx) > 255:
        raise ValueError("ctx demasiado largo (máx. 255 bytes)")
    return b"\x00" + bytes([len(ctx)]) + ctx + M


def sign(sk, M, ctx=b"", params_name="ML-DSA-44", rnd=None, stats=None):
    """Algorithm 2: firma 'hedged' (aleatorizada, por defecto) o determinista."""
    if rnd is None:
        rnd = os.urandom(32)              # variante hedged (recomendada)
    M_prime = _format_message(M, ctx)
    return sign_internal(sk, M_prime, rnd, PARAMS[params_name], stats)


def verify(pk, M, sig, ctx=b"", params_name="ML-DSA-44"):
    """Algorithm 3: verificación."""
    M_prime = _format_message(M, ctx)
    return verify_internal(pk, M_prime, sig, PARAMS[params_name])


# ----------------------------------------------------------------------------
# HashML-DSA: variante pre-hash (FIPS 204, Sección 5.4; Algoritmos 4 y 5)
# ----------------------------------------------------------------------------

# Funciones de pre-hash aprobadas: nombre -> (función, DER encoding del OID)
_HASHLIB = {
    "SHA2-224": lambda m: hashlib.sha224(m).digest(),
    "SHA2-256": lambda m: hashlib.sha256(m).digest(),
    "SHA2-384": lambda m: hashlib.sha384(m).digest(),
    "SHA2-512": lambda m: hashlib.sha512(m).digest(),
    "SHA2-512/224": lambda m: hashlib.new("sha512_224", m).digest(),
    "SHA2-512/256": lambda m: hashlib.new("sha512_256", m).digest(),
    "SHA3-224": lambda m: hashlib.sha3_224(m).digest(),
    "SHA3-256": lambda m: hashlib.sha3_256(m).digest(),
    "SHA3-384": lambda m: hashlib.sha3_384(m).digest(),
    "SHA3-512": lambda m: hashlib.sha3_512(m).digest(),
    "SHAKE-128": lambda m: hashlib.shake_128(m).digest(32),
    "SHAKE-256": lambda m: hashlib.shake_256(m).digest(64),
}

# DER encoding (con tag y longitud) de los OID de 2.16.840.1.101.3.4.2.x
def _der_oid(last_byte):
    return bytes([0x06, 0x09, 0x60, 0x86, 0x48, 0x01, 0x65, 0x03, 0x04, 0x02,
                  last_byte])

PH_OID = {
    "SHA2-224": _der_oid(0x04),
    "SHA2-256": _der_oid(0x01),   # 2.16.840.1.101.3.4.2.1 (ejemplo del FIPS)
    "SHA2-384": _der_oid(0x02),
    "SHA2-512": _der_oid(0x03),
    "SHA2-512/224": _der_oid(0x05),
    "SHA2-512/256": _der_oid(0x06),
    "SHA3-224": _der_oid(0x07),
    "SHA3-256": _der_oid(0x08),
    "SHA3-384": _der_oid(0x09),
    "SHA3-512": _der_oid(0x0A),
    "SHAKE-128": _der_oid(0x0B),  # ejemplo del FIPS 204 (salida de 256 bits)
    "SHAKE-256": _der_oid(0x0C),
}


def _format_prehash(M, ctx, ph):
    """Construye M′ = 0x01 ‖ len(ctx) ‖ ctx ‖ OID_PH ‖ PH(M) (Sección 5.4.1)."""
    if ph not in _HASHLIB:
        raise ValueError(f"función de pre-hash no soportada: {ph}")
    if len(ctx) > 255:
        raise ValueError("ctx demasiado largo (máx. 255 bytes)")
    return (b"\x01" + bytes([len(ctx)]) + ctx + PH_OID[ph] + _HASHLIB[ph](M))


def hash_sign(sk, M, ctx=b"", ph="SHAKE-256", params_name="ML-DSA-44",
              rnd=None, stats=None):
    """Algorithm 4: HashML-DSA.Sign — firma con pre-hash (hedged o determinista)."""
    if rnd is None:
        rnd = os.urandom(32)
    M_prime = _format_prehash(M, ctx, ph)
    return sign_internal(sk, M_prime, rnd, PARAMS[params_name], stats)


def hash_verify(pk, M, sig, ctx=b"", ph="SHAKE-256", params_name="ML-DSA-44"):
    """Algorithm 5: HashML-DSA.Verify — verificación con pre-hash."""
    M_prime = _format_prehash(M, ctx, ph)
    return verify_internal(pk, M_prime, sig, PARAMS[params_name])


# ----------------------------------------------------------------------------
# Interfaz ExternalMu: μ calculado fuera del módulo firmante
# (FIPS 204, notas del Algoritmo 7 línea 6 y Algoritmo 8 línea 7: el
# representante del mensaje μ puede computarse en un módulo distinto)
# ----------------------------------------------------------------------------

def compute_mu(pk, M, ctx=b"", params_name="ML-DSA-44", ph=None):
    """μ := H(H(pk) ‖ M′), con M′ en formato puro (0x00) o pre-hash (0x01)."""
    tr = shake256(pk, 64)
    M_prime = _format_prehash(M, ctx, ph) if ph else _format_message(M, ctx)
    return shake256(tr + M_prime, 64)


def sign_mu(sk, mu, rnd=None, params_name="ML-DSA-44", stats=None):
    """Firma a partir de μ precalculado (la clave privada no ve el mensaje)."""
    if rnd is None:
        rnd = os.urandom(32)
    return sign_internal(sk, b"", rnd, PARAMS[params_name], stats, mu=mu)


def verify_mu(pk, mu, sig, params_name="ML-DSA-44"):
    """Verifica una firma a partir de μ precalculado."""
    return verify_internal(pk, b"", sig, PARAMS[params_name], mu=mu)


# Tamaños oficiales (FIPS 204, Tabla 2) para comprobación
SIZES = {
    "ML-DSA-44": {"pk": 1312, "sk": 2560, "sig": 2420},
    "ML-DSA-65": {"pk": 1952, "sk": 4032, "sig": 3309},
    "ML-DSA-87": {"pk": 2592, "sk": 4896, "sig": 4627},
}


if __name__ == "__main__":
    for name in PARAMS:
        pk, sk = keygen(name)
        msg = b"Demostracion ML-DSA - FIPS 204"
        sig = sign(sk, msg, params_name=name)
        ok = verify(pk, msg, sig, params_name=name)
        s = SIZES[name]
        assert len(pk) == s["pk"] and len(sk) == s["sk"] and len(sig) == s["sig"]

        # Variante pre-hash (HashML-DSA)
        sig_ph = hash_sign(sk, msg, ph="SHA2-512", params_name=name)
        ok_ph = hash_verify(pk, msg, sig_ph, ph="SHA2-512", params_name=name)

        # Interfaz ExternalMu
        mu = compute_mu(pk, msg, params_name=name)
        sig_mu = sign_mu(sk, mu, params_name=name)
        ok_mu = verify_mu(pk, mu, sig_mu, params_name=name)

        print(f"{name}: pk={len(pk)}B sk={len(sk)}B sig={len(sig)}B "
              f"pura={'OK' if ok else 'FALLO'} "
              f"prehash={'OK' if ok_ph else 'FALLO'} "
              f"externalMu={'OK' if ok_mu else 'FALLO'}")
