"""
Elliptic Curve Cryptography (ECC)

Module to generate Public Key from elliptic curve secp256k1.

(Curve used in Bitcoin cryptocurrency)

- Source:

            https://github.com/gtmadureira/cryptography/ecc/secp256k1/blob/main/secp256k1.py

- Author:

            • Gustavo Madureira (gtmadureira@gmail.com)
            • https://gtmadureira.github.io/
"""


# Mathematical domain parameters of the elliptic curve 'secp256k1'.
# Source: https://www.secg.org/sec2-v2.pdf
_FP_CURVE = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
_A_CURVE = 0x0000000000000000000000000000000000000000000000000000000000000000
_B_CURVE = 0x0000000000000000000000000000000000000000000000000000000000000007
_GX_CURVE = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
_GY_CURVE = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
_N_CURVE = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
_H_CURVE = 0x0000000000000000000000000000000000000000000000000000000000000001


def modular_inverse(k: int, p: int) -> int:
    """
    Extended Euclidean algorithm/'division' in elliptic curve.

    Returns the multiplicative inverse of k modulo p. Where the only
    integer x is defined such that (x * k) % p == 1.

    k must be non-zero and p must be a prime.
    """
    if k == 0:
        raise ZeroDivisionError("Division by zero!")
    if k < 0:
        # k ** -1 = p - (-k) ** -1  (mod p)
        return p - modular_inverse(-k, p)
    # Extended Euclidean algorithm.
    s, old_s = 0, 1
    t, old_t = 1, 0
    r, old_r = p, k
    while r != 0:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t
    gcd, x, y = old_r, old_s, old_t
    assert gcd == 1
    assert (k * x) % p == 1
    assert (p * y) % k == 1
    return x % p


def ec_point_doubling(px: int, py: int) -> tuple:
    """
    Point doubling in elliptic curve.

    It doubles Point-P.
    """
    mn = 3 * px * px + _A_CURVE
    md = 2 * py
    mr = (mn * modular_inverse(md, _FP_CURVE)) % _FP_CURVE
    xr = (mr * mr - 2 * px) % _FP_CURVE
    yr = (mr * (px - xr) - py) % _FP_CURVE
    return (xr, yr)


def ec_point_addition(px: int, py: int, qx: int, qy: int) -> tuple:
    """
    Point addition in elliptic curve.

    It adds Point-P with Point-Q.
    """
    if px == qx:
        mr = ((3 * (px * px) + _A_CURVE) *
              modular_inverse(2 * py, _FP_CURVE)) % _FP_CURVE
    else:
        mr = ((qy - py) * modular_inverse(qx - px, _FP_CURVE)) % _FP_CURVE
    xr = (mr * mr - px - qx) % _FP_CURVE
    yr = (mr * (px - xr) - py) % _FP_CURVE
    return (xr, yr)


def ec_point_multiplication(scalar: int) -> tuple:
    """
    Point multiplication in elliptic curve.

    It doubles Point-P and adds Point-P with Point-Q.
    """
    if not 0 < scalar < _N_CURVE:
        raise Exception("Invalid Scalar/Private Key")
    scalarbin = bin(scalar)[2:]
    px, py = _GX_CURVE, _GY_CURVE
    qx, qy = px, py
    for i in range(1, len(scalarbin)):
        qx, qy = ec_point_doubling(qx, qy)
        if scalarbin[i] == "1":
            qx, qy = ec_point_addition(px, py, qx, qy)
    return (qx, qy)


if __name__ == "__main__":
    private_key = 1
    while True:
        print()
        print("             Point Number: " + str(private_key))
        print("              Private Key: " +
              hex(private_key)[2:].zfill(64).upper())
        public_key = ec_point_multiplication(private_key)
        if public_key[1] % 2 == 1:
            i = "03"
        else:
            i = "02"
        print("  Uncompressed Public Key: " + "04" +
              hex(public_key[0])[2:].zfill(64).upper() +
              hex(public_key[1])[2:].zfill(64).upper())
        print("    Compressed Public Key: " + i +
              hex(public_key[0])[2:].zfill(64).upper())
        print()
        private_key += 1
