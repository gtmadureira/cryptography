"""
Elliptic Curve Cryptography (ECC)

Module to generate Public Key from elliptic curve secp256k1.

(Curve used in Bitcoin cryptocurrency)

- Source:

            https://github.com/gtmadureira/cryptography/blob/main/ecc/secp256k1/secp256k1.py

- Author:

            • Gustavo Madureira (gtmadureira@gmail.com)
            • https://gtmadureira.github.io/
"""


from typing import Tuple


# Type Hints.
Point = Tuple[int, int]


"""
    Mathematical domain parameters of the elliptic curve secp256k1.
    Source: https://www.secg.org/sec2-v2.pdf
"""

# Finite field (Fp):
_FP_CURVE_ = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F

# The elliptic curve: (y^2 = x^3 + ax + b) over Fp is defined by:
_A_CURVE_ = 0x0000000000000000000000000000000000000000000000000000000000000000
_B_CURVE = 0x0000000000000000000000000000000000000000000000000000000000000007

# The generator point in uncompressed form is:
_GX_CURVE_ = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
_GY_CURVE_ = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
_GENERATOR_POINT_CURVE: Point = (_GX_CURVE_, _GY_CURVE_)

# Finally the order of generator point and the cofactor are:
_N_CURVE_ = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
_H_CURVE_ = 0x0000000000000000000000000000000000000000000000000000000000000001


def modular_inverse(k: int, p: int) -> int:
    """
    Extended Euclidean algorithm/division in elliptic curve.

    Returns the multiplicative inverse of k modulo p. Where the only
    integer x is defined such that (x * k) % p == 1.

    k must be non-zero and p must be a prime.
    """
    if k == 0:
        raise ZeroDivisionError("Division by zero!")
    if k < 0:
        return p - modular_inverse(-k, p)
    s, old_s = (0, 1)
    t, old_t = (1, 0)
    r, old_r = (p, k)
    while r != 0:
        quotient = old_r // r
        old_r, r = (r, old_r - quotient * r)
        old_s, s = (s, old_s - quotient * s)
        old_t, t = (t, old_t - quotient * t)
    gcd, x, y = (old_r, old_s, old_t)
    assert gcd == 1
    assert (k * x) % p == 1
    assert (p * y) % k == 1
    return x % p


def is_on_curve(point: Point) -> bool:
    """Returns True if the given point lies on the elliptic curve."""
    if point is None:
        return True
    px, py = point
    result = (py * py - px * px * px - _A_CURVE_ *
              px - _B_CURVE) % _FP_CURVE_ == 0
    return result


def ec_point_doubling(point_P: Point) -> Point:
    """
    Point doubling in elliptic curve.

    It doubles Point-P.
    """
    assert is_on_curve(point_P)
    if point_P is None:
        return None  # type: ignore
    px, py = point_P
    slope = ((3 * px * px + _A_CURVE_) *
             modular_inverse(2 * py, _FP_CURVE_)) % _FP_CURVE_
    rx = (slope * slope - 2 * px) % _FP_CURVE_
    ry = (slope * (px - rx) - py) % _FP_CURVE_
    new_point = (rx, ry)
    assert is_on_curve(new_point)
    return new_point


def ec_point_addition(point_P: Point, point_Q: Point) -> Point:
    """
    Point addition in elliptic curve.

    It adds Point-P with Point-Q.
    """
    assert is_on_curve(point_P)
    assert is_on_curve(point_Q)
    if point_P is None:
        return point_Q
    if point_Q is None:
        return point_P
    px, py = point_P
    qx, qy = point_Q
    if px == qx and py != qy:
        return None  # type: ignore
    if px == qx and py == qy:
        return ec_point_doubling(point_P)
    slope = ((py - qy) * modular_inverse(px - qx, _FP_CURVE_)) % _FP_CURVE_
    rx = (slope * slope - px - qx) % _FP_CURVE_
    ry = (slope * (px - rx) - py) % _FP_CURVE_
    new_point = (rx, ry)
    assert is_on_curve(new_point)
    return new_point


def ec_point_multiplication(scalar: int,
                            point: Point = _GENERATOR_POINT_CURVE) -> Point:
    """
    Point multiplication in elliptic curve.

    It doubles Point-P and adds Point-P with Point-Q.
    """
    assert is_on_curve(point)
    if not 0 < scalar < _N_CURVE_:
        return None  # type: ignore
    if point is None:
        return None  # type: ignore
    scalarbin = bin(scalar)[2:]
    result = None
    current = point
    for i in range(1, len(scalarbin)):
        current = ec_point_doubling(current)
        if scalarbin[i] == "1":
            current = ec_point_addition(point, current)
    result = current
    assert is_on_curve(result)
    return result


if __name__ == "__main__":
    private_key = 1
    while True:
        public_key = ec_point_multiplication(private_key)
        if public_key[1] % 2 == 1:
            i = "03"
        else:
            i = "02"
        print("\n\033[92m" +
              "           Point Number: " + str(private_key) + "\n"
              "            Private Key: " +
              hex(private_key)[2:].zfill(64).upper() + "\n"
              "Uncompressed Public Key: " + "04" +
              hex(public_key[0])[2:].zfill(64).upper() +
              hex(public_key[1])[2:].zfill(64).upper() + "\n"
              "  Compressed Public Key: " + i +
              hex(public_key[0])[2:].zfill(64).upper() + "\n" +
              "\033[0m")
        private_key += 1
