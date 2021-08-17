"""
[ LICENSE ] ( GNU GPLv3 ) :

SECP256K1 at Weierstrass Form - Elliptic Curve Cryptography (ECC).
Copyright (C) 2021  Gustavo Madureira

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see https://www.gnu.org/licenses/ .


[ ABOUT ] :

Elliptic Curve Cryptography (ECC)

Module to generate Public Key from elliptic curve secp256k1,
using Weierstrass Form.

(Curve used in Bitcoin cryptocurrency)

    Source:

            https://github.com/gtmadureira/cryptography/blob/main/ecc/secp256k1/secp256k1_weierstrass.py

    Author:

            • Gustavo Madureira (gtmadureira@gmail.com)
            • https://gtmadureira.github.io/
"""


from typing import Optional, Tuple


# Type Hints.
Point = Tuple[int, int]


def int_from_hex(hexadecimal_string: str) -> int:
    """Converts the hexadecimal string to integer."""
    result = int(
        "0x" + "".join(hexadecimal_string.replace(":", "").split()), 16)
    return result


#       Mathematical domain parameters of the elliptic curve secp256k1.
#       Source: https://www.secg.org/sec2-v2.pdf


# The finite field (Fp) is defined by:
_FP_CURVE_ = int_from_hex(
    "FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFE FFFFFC2F")


# The elliptic curve (y^2 = x^3 + ax + b) over Fp is defined by:
_A_CURVE_ = int_from_hex(
    "00000000 00000000 00000000 00000000 00000000 00000000 00000000 00000000")

_B_CURVE_ = int_from_hex(
    "00000000 00000000 00000000 00000000 00000000 00000000 00000000 00000007")


# The generator point is defined by:
_GX_CURVE_ = int_from_hex(
    "79BE667E F9DCBBAC 55A06295 CE870B07 029BFCDB 2DCE28D9 59F2815B 16F81798")

_GY_CURVE_ = int_from_hex(
    "483ADA77 26A3C465 5DA4FBFC 0E1108A8 FD17B448 A6855419 9C47D08F FB10D4B8")

_GENERATOR_POINT_CURVE_: Point = (_GX_CURVE_, _GY_CURVE_)


# The order of generator point and the cofactor are defined by:
_N_CURVE_ = int_from_hex(
    "FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFE BAAEDCE6 AF48A03B BFD25E8C D0364141")

_H_CURVE_ = int_from_hex(
    "00000000 00000000 00000000 00000000 00000000 00000000 00000000 00000001")


# The point that points to infinity in elliptic curve is defined by:
_POINT_INFINITY_CURVE_ = None


def modular_inverse(k: int, p: int) -> int:
    """
    Extended Euclidean algorithm/division in elliptic curve.

    Returns the multiplicative inverse of k modulo p. Where the only
    integer x is defined such that (k * x) % p == 1.

    k must be non-zero and p must be a prime.
    """
    if k == 0:
        raise ZeroDivisionError("Division by zero!")
    if k < 0:
        result = p - modular_inverse(-k, p)
        return result
    old_r, r = (k, p)
    old_s, s = (1, 0)
    old_t, t = (0, 1)
    while old_r != 0:
        quotient = r // old_r
        r, old_r = (old_r, r - quotient * old_r)
        s, old_s = (old_s, s - quotient * old_s)
        t, old_t = (old_t, t - quotient * old_t)
    gcd, x, y = (r, s, t)
    if k == 1:
        assert gcd == 1
        assert (k * x) % p == 1
        assert (p * y) % k == 0
    if k == p:
        assert gcd == p
        assert (k * x) % p == 0
        assert (p * y) % k == 0
    if k != 1 and k != p:
        assert gcd == 1
        assert (k * x) % p == 1
        assert (p * y) % k == 1
    result = x % p
    return result


def is_infinite(point: Optional[Point]) -> bool:
    """
    Returns whether or not it is the point at infinity in elliptic
    curve.
    """
    result = point is _POINT_INFINITY_CURVE_
    return result


def is_on_curve(point: Optional[Point]) -> bool:
    """Returns True if the given point lies on the elliptic curve."""
    if is_infinite(point):
        result = True
        return result
    xp, yp = point  # type: ignore
    result = (yp ** 2 - xp ** 3 - _A_CURVE_ * xp - _B_CURVE_) % _FP_CURVE_ == 0
    return result


def x(point: Point) -> int:
    """
    Refer to the x coordinate of a point (assuming it is not infinity),
    then returns this value.
    """
    assert not is_infinite(point)
    assert is_on_curve(point)
    result = point[0]
    return result


def y(point: Point) -> int:
    """
    Refer to the y coordinate of a point (assuming it is not infinity),
    then returns this value.
    """
    assert not is_infinite(point)
    assert is_on_curve(point)
    result = point[1]
    return result


def has_even_y(point: Point) -> bool:
    """
    Where the point is not is_infinite(point),
    it returns y(point) mod 2 = 0.
    """
    assert not is_infinite(point)
    assert is_on_curve(point)
    result = y(point) % 2 == 0
    return result


def ec_point_doubling(point_p: Optional[Point]) -> Optional[Point]:
    """
    Point doubling in elliptic curve.

    It doubles Point-P.
    """
    assert is_on_curve(point_p)
    if is_infinite(point_p):
        result = _POINT_INFINITY_CURVE_
        return result
    xp, yp = point_p  # type: ignore
    if 0 in {xp, yp}:
        result = _POINT_INFINITY_CURVE_
        return result
    slope = ((3 * xp ** 2 + _A_CURVE_) *
             modular_inverse(2 * yp, _FP_CURVE_)) % _FP_CURVE_
    xr = (slope ** 2 - 2 * xp) % _FP_CURVE_
    yr = (slope * (xp - xr) - yp) % _FP_CURVE_
    result = (xr, yr)
    assert is_on_curve(result)
    return result


def ec_point_addition(point_p: Optional[Point],
                      point_q: Optional[Point]) -> Optional[Point]:
    """
    Point addition in elliptic curve.

    It adds Point-P with Point-Q.
    """
    assert is_on_curve(point_p)
    assert is_on_curve(point_q)
    if is_infinite(point_p):
        result = point_q
        return result
    if is_infinite(point_q):
        result = point_p
        return result
    xp, yp = point_p  # type: ignore
    xq, yq = point_q  # type: ignore
    if xp == xq and yp != yq:
        result = _POINT_INFINITY_CURVE_
        return result
    if xp == xq and yp == yq:
        if 0 in {xp, yp}:
            result = _POINT_INFINITY_CURVE_
            return result
        else:
            result = ec_point_doubling(point_p)
            return result
    slope = ((yq - yp) * modular_inverse(xq - xp, _FP_CURVE_)) % _FP_CURVE_
    xr = (slope ** 2 - xp - xq) % _FP_CURVE_
    yr = (slope * (xp - xr) - yp) % _FP_CURVE_
    result = (xr, yr)
    assert is_on_curve(result)
    return result


def ec_point_multiplication(scalar: int, point: Optional[Point]) -> Point:
    """
    Scalar point multiplication in elliptic curve.

    It doubles Point-P and adds Point-P with Point-Q.
    """
    assert is_on_curve(point)
    if not 0 < scalar < _N_CURVE_:
        raise Exception("Invalid Scalar/Private Key!")
    if point is None or 0 in {point[0], point[1]}:
        raise Exception("""
            None (Generator/Base Point) has been provided
            or points to infinity on the elliptic curve!
            """)
    result = None
    scalarbin = bin(scalar)[2:]
    current = point
    for i in range(1, len(scalarbin)):
        current = ec_point_doubling(current)  # type: ignore
        if scalarbin[i] == "1":
            current = ec_point_addition(point, current)  # type: ignore
    result = current
    assert is_on_curve(result)
    return result


if __name__ == "__main__":

    # Elliptic curve point multiplication test.
    from platform import system
    from time import sleep, perf_counter
    from subprocess import check_call as run_command

    start = perf_counter()

    # Tests the operating system type and sets the screen clear command.
    if system() == "Windows":

        def clear() -> None:
            """Screen clear command for Windows operating system."""
            run_command("cls")

    elif system() == "Darwin" or system() == "Linux":

        def clear() -> None:
            """Screen clear command for macOS/Linux operating system."""
            run_command("clear")

    private_key = 1
    while private_key < 10001:
        public_key = ec_point_multiplication(
            private_key, _GENERATOR_POINT_CURVE_)
        if has_even_y(public_key):
            prefix = "02"
        else:
            prefix = "03"
        data = (str(private_key),
                hex(private_key)[2:].zfill(64).upper(),
                hex(x(public_key))[2:].zfill(64).upper(),
                hex(y(public_key))[2:].zfill(64).upper(),
                prefix)
        private_key += 1
        sleep(0.0)
        clear()
        print(f"""\033[92m
        SECP256K1 at Weierstrass Form  Copyright (C) 2021  Gustavo Madureira
        This program comes with ABSOLUTELY NO WARRANTY.
        This is free software, and you are welcome to redistribute it
        under certain conditions.


            Point Number: {data[0]}
             Private Key: {data[1]}
 Uncompressed Public Key: 04{data[2]}{data[3]}
   Compressed Public Key: {data[4]}{data[2]}
\033[0m""")
    elapsed = perf_counter() - start
    print(f"""\033[92m
        Finished in {elapsed:.02f} seconds.
\033[0m""")
