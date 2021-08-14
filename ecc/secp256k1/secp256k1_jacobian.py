"""
[ LICENSE ] ( GNU GPLv3 ) :

SECP256K1 at Jacobian Form - Elliptic Curve Cryptography (ECC).
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
using Jacobian Form.

(Curve used in Bitcoin cryptocurrency)

    Source:

            https://github.com/gtmadureira/cryptography/blob/main/ecc/secp256k1/secp256k1_jacobian.py

    Author:

            • Gustavo Madureira (gtmadureira@gmail.com)
            • https://gtmadureira.github.io/
"""


from typing import Optional, Tuple


# Type Hints.
Point = Tuple[int, int]
Jacobian_Coordinate = Tuple[int, int, int]


"""
        Mathematical domain parameters of the elliptic curve secp256k1.
        Source: https://www.secg.org/sec2-v2.pdf
"""

# Finite field (Fp):
_FP_CURVE_ = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F

# The elliptic curve (y^2 = x^3 + ax + b) over Fp is defined by:
_A_CURVE_ = 0x0000000000000000000000000000000000000000000000000000000000000000
_B_CURVE_ = 0x0000000000000000000000000000000000000000000000000000000000000007

# The generator point in uncompressed form is:
_GX_CURVE_ = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
_GY_CURVE_ = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
_GENERATOR_POINT_CURVE_: Point = (_GX_CURVE_, _GY_CURVE_)

# Order of generator point and the cofactor are:
_N_CURVE_ = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
_H_CURVE_ = 0x0000000000000000000000000000000000000000000000000000000000000001

# Definition for a point that points to infinity in elliptic curve:
_POINT_INFINITY_CURVE_ = None


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


def to_jacobian(point: Point) -> Jacobian_Coordinate:
    """
    Convert an affine point to Jacobian coordinate.

    A Jacobian coordinate is represented as (x, y, z).
    """
    xp, yp = point
    xp, yp, zp = (xp, yp, 1)
    result = (xp, yp, zp)
    return result


def from_jacobian(point: Jacobian_Coordinate) -> Optional[Point]:
    """
    Convert a Jacobian coordinate to affine point, or None if at
    infinity.

    An affine point is represented as (x, y).
    """
    xp, yp, zp = point
    if zp == 0:
        result = None
        return result
    if zp == 1:
        result = (xp, yp)
        return result
    inv_1 = modular_inverse(zp, _FP_CURVE_)
    inv_2 = (inv_1 ** 2) % _FP_CURVE_
    inv_3 = (inv_2 * inv_1) % _FP_CURVE_
    result = ((inv_2 * xp) % _FP_CURVE_, (inv_3 * yp) % _FP_CURVE_)
    return result


def jacobian_point_doubling(
        point_p: Jacobian_Coordinate) -> Jacobian_Coordinate:
    """
    Point doubling in elliptic curve, using Jacobian coordinate
    (x, y, z).

    It doubles Point-P.
    """
    xp, yp, zp = point_p
    if zp == 0:
        result = (0, 1, 0)
        return result
    xp_2 = (xp ** 2) % _FP_CURVE_
    yp_2 = (yp ** 2) % _FP_CURVE_
    yp_4 = (yp_2 ** 2) % _FP_CURVE_
    s = (4 * xp * yp_2) % _FP_CURVE_
    m = 3 * xp_2
    if _A_CURVE_:
        m += _A_CURVE_ * pow(zp, 4, _FP_CURVE_)
    m = m % _FP_CURVE_
    xr = (m ** 2 - 2 * s) % _FP_CURVE_
    yr = (m * (s - xr) - 8 * yp_4) % _FP_CURVE_
    zr = (2 * yp * zp) % _FP_CURVE_
    result = (xr, yr, zr)
    return result


def jacobian_point_addition_mixed(
        point_p: Jacobian_Coordinate,
        point_q: Jacobian_Coordinate) -> Jacobian_Coordinate:
    """
    Point addition (mixed) in elliptic curve, using Jacobian coordinate
    (x, y, z) and affine coordinate in Jacobian form (x, y, 1).

    It adds Point-P with Point-Q.
    """
    xp, yp, zp = point_p
    xq, yq, zq = point_q
    assert zq == 1
    if zp == 0:
        result = point_q
        return result
    zp_2 = (zp ** 2) % _FP_CURVE_
    zp_3 = (zp_2 * zp) % _FP_CURVE_
    uq = (xq * zp_2) % _FP_CURVE_
    sq = (yq * zp_3) % _FP_CURVE_
    if xp == uq:
        if yp != sq:
            result = (0, 1, 0)
            return result
        result = jacobian_point_doubling(point_p)
        return result
    h = uq - xp
    r = sq - yp
    h_2 = (h ** 2) % _FP_CURVE_
    h_3 = (h_2 * h) % _FP_CURVE_
    xp_h_2 = (xp * h_2) % _FP_CURVE_
    xr = (r ** 2 - h_3 - 2 * xp_h_2) % _FP_CURVE_
    yr = (r * (xp_h_2 - xr) - yp * h_3) % _FP_CURVE_
    zr = (h * zp) % _FP_CURVE_
    result = (xr, yr, zr)
    return result


def jacobian_point_addition(
        point_p: Jacobian_Coordinate,
        point_q: Jacobian_Coordinate) -> Jacobian_Coordinate:
    """
    Point addition in elliptic curve, using Jacobian coordinate
    (x, y, z).

    It adds Point-P with Point-Q.
    """
    xp, yp, zp = point_p
    xq, yq, zq = point_q
    if zp == 0:
        result = point_q
        return result
    if zq == 0:
        result = point_p
        return result
    if zp == 1:
        result = jacobian_point_addition_mixed(point_q, point_p)
        return result
    if zq == 1:
        result = jacobian_point_addition_mixed(point_p, point_q)
        return result
    zp_2 = (zp ** 2) % _FP_CURVE_
    zp_3 = (zp_2 * zp) % _FP_CURVE_
    zq_2 = (zq ** 2) % _FP_CURVE_
    zq_3 = (zq_2 * zq) % _FP_CURVE_
    up = (xp * zq_2) % _FP_CURVE_
    uq = (xq * zp_2) % _FP_CURVE_
    sp = (yp * zq_3) % _FP_CURVE_
    sq = (yq * zp_3) % _FP_CURVE_
    if up == uq:
        if sp != sq:
            result = (0, 1, 0)
            return result
        result = jacobian_point_doubling(point_p)
        return result
    h = uq - up
    r = sq - sp
    h_2 = (h ** 2) % _FP_CURVE_
    h_3 = (h_2 * h) % _FP_CURVE_
    up_h_2 = (up * h_2) % _FP_CURVE_
    xr = (r ** 2 - h_3 - 2 * up_h_2) % _FP_CURVE_
    yr = (r * (up_h_2 - xr) - sp * h_3) % _FP_CURVE_
    zr = (h * zp * zq) % _FP_CURVE_
    result = (xr, yr, zr)
    return result


def jacobian_point_multiplication(scalar: int,
                                  point: Optional[Point]) -> Point:
    """
    Point multiplication in elliptic curve, using Jacobian coordinate
    (x, y, z).

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
    jacobian = to_jacobian(point)
    current = jacobian
    for i in range(1, len(scalarbin)):
        current = jacobian_point_doubling(current)
        if scalarbin[i] == "1":
            current = jacobian_point_addition(jacobian, current)
    affine = from_jacobian(current)
    result = affine
    assert is_on_curve(result)
    return result  # type: ignore


if __name__ == "__main__":

    from time import sleep
    from platform import system
    from subprocess import check_call as run_command

    # Tests the operating system type and sets the screen clear command.
    if system() == "Windows":

        def clear() -> None:
            """Screen clear command for Windows operating system."""
            run_command("cls")

    elif system() == "Darwin" or system() == "Linux":

        def clear() -> None:
            """Screen clear command for macOS/Linux operating system."""
            run_command("clear")

    # Elliptic curve point multiplication test.
    private_key = 1
    while True:
        public_key = jacobian_point_multiplication(
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
        sleep(0.65)
        clear()
        print(f"""\033[92m
        SECP256K1 at Jacobian Form  Copyright (C) 2021  Gustavo Madureira
        This program comes with ABSOLUTELY NO WARRANTY.
        This is free software, and you are welcome to redistribute it
        under certain conditions.


            Point Number: {data[0]}
             Private Key: {data[1]}
 Uncompressed Public Key: 04{data[2]}{data[3]}
   Compressed Public Key: {data[4]}{data[2]}
\033[0m""")
