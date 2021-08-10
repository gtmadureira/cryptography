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


from time import sleep
from platform import system
from typing import Optional, Tuple
from subprocess import check_call as run_command


# Type Hints.
Point = Tuple[int, int]
Jacobian_Coordinate = Tuple[int, int, int]


# Tests the operating system type and sets the screen clear command.
if system() == "Windows":

    def clear() -> None:
        """Screen clear command for Windows operating system."""
        run_command("cls")

elif system() == "Darwin" or system() == "Linux":

    def clear() -> None:
        """Screen clear command for macOS/Linux operating system."""
        run_command("clear")


"""
        Mathematical domain parameters of the elliptic curve secp256k1.
        Source: https://www.secg.org/sec2-v2.pdf
"""

# Finite field (Fp):
_FP_CURVE_ = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F

# The elliptic curve: (y^2 = x^3 + ax + b) over Fp is defined by:
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
_INFINITE_POINT_CURVE_ = None


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
    r, old_r = (p, k)
    s, old_s = (0, 1)
    t, old_t = (1, 0)
    while r != 0:
        quotient = old_r // r
        old_r, r = (r, old_r - quotient * r)
        old_s, s = (s, old_s - quotient * s)
        old_t, t = (t, old_t - quotient * t)
    gcd, x, y = (old_r, old_s, old_t)
    assert gcd == 1
    assert (k * x) % p == 1
    if k == 1:
        assert (p * y) % k == 0
    if k != 1:
        assert (p * y) % k == 1
    result = x % p
    return result


def is_infinite(point: Optional[Point]) -> bool:
    """
    Returns whether or not it is the point at infinity in elliptic
    curve.
    """
    result = point is _INFINITE_POINT_CURVE_
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
    """Convert from Weierstrass form to Jacobian form."""
    o = (point[0], point[1], 1)
    result = o
    return result


def jacobian_doubling(point_p: Jacobian_Coordinate) -> Jacobian_Coordinate:
    """
    Point doubling in elliptic curve.

    It doubles Point-P.
    """
    if not point_p[1]:
        result = (0, 0, 0)
        return result
    ysq = (point_p[1] ** 2) % _FP_CURVE_
    s = (4 * point_p[0] * ysq) % _FP_CURVE_
    m = (3 * point_p[0] ** 2 + _A_CURVE_ * point_p[2] ** 4) % _FP_CURVE_
    nx = (m ** 2 - 2 * s) % _FP_CURVE_
    ny = (m * (s - nx) - 8 * ysq ** 2) % _FP_CURVE_
    nz = (2 * point_p[1] * point_p[2]) % _FP_CURVE_
    result = (nx, ny, nz)
    return result


def jacobian_addition(point_p: Jacobian_Coordinate,
                      point_q: Jacobian_Coordinate) -> Jacobian_Coordinate:
    """
    Point addition in elliptic curve.

    It adds Point-P with Point-Q.
    """
    if not point_p[1]:
        result = point_q
        return result
    if not point_q[1]:
        result = point_p
        return result
    u1 = (point_p[0] * point_q[2] ** 2) % _FP_CURVE_
    u2 = (point_q[0] * point_p[2] ** 2) % _FP_CURVE_
    s1 = (point_p[1] * point_q[2] ** 3) % _FP_CURVE_
    s2 = (point_q[1] * point_p[2] ** 3) % _FP_CURVE_
    if u1 == u2:
        if s1 != s2:
            result = (0, 0, 1)
            return result
        result = jacobian_doubling(point_p)
        return result
    h = u2 - u1
    r = s2 - s1
    h2 = (h * h) % _FP_CURVE_
    h3 = (h * h2) % _FP_CURVE_
    u1h2 = (u1 * h2) % _FP_CURVE_
    nx = (r ** 2 - h3 - 2 * u1h2) % _FP_CURVE_
    ny = (r * (u1h2 - nx) - s1 * h3) % _FP_CURVE_
    nz = (h * point_p[2] * point_q[2]) % _FP_CURVE_
    result = (nx, ny, nz)
    return result


def from_jacobian(point: Jacobian_Coordinate) -> Point:
    """Convert from Jacobian form to Weierstrass form."""
    z = modular_inverse(point[2], _FP_CURVE_)
    result = ((point[0] * z ** 2) % _FP_CURVE_,
              (point[1] * z ** 3) % _FP_CURVE_)
    return result


def jacobian_multiplication(point: Jacobian_Coordinate,
                            scalar: int) -> Jacobian_Coordinate:
    """
    Point multiplication in elliptic curve.

    It doubles Point-P and adds Point-P with Point-Q.
    """
    result = None
    if point[1] == 0 or scalar == 0:
        result = (0, 0, 1)
        return result
    if scalar == 1:
        result = point
        return result
    if (scalar % 2) == 0:
        result = jacobian_doubling(jacobian_multiplication(point, scalar // 2))
        return result
    if (scalar % 2) == 1:
        result = jacobian_addition(jacobian_doubling(
            jacobian_multiplication(point, scalar // 2)), point)
        return result
    return result  # type: ignore


def ec_point_multiplication(
        scalar: int,
        point: Optional[Point] = _GENERATOR_POINT_CURVE_) -> Point:
    """
    It prepares the generator point, from Weierstrass format converting
    to Jacobian coordinate before starting the multiplication, and at
    the end converts it again but from Jacobian coordinate to
    Weierstrass format.
    """
    assert is_on_curve(point)
    if not 0 < scalar < _N_CURVE_:
        raise Exception("Invalid Scalar/Private Key!")
    if point is None or point[0] == 0 or point[1] == 0:
        raise Exception("None (Generator/Base Point) has been provided " +
                        "or points to infinity on the elliptic curve!")
    result = None
    result = from_jacobian(jacobian_multiplication(to_jacobian(point), scalar))
    assert is_on_curve(result)
    return result


if __name__ == "__main__":
    # Elliptic curve point multiplication test.
    private_key = 1
    while True:
        public_key = ec_point_multiplication(private_key)
        if has_even_y(public_key):
            i = "02"
        else:
            i = "03"
        sleep(0.65)
        clear()
        print("\n\033[92m" +
              "\t\tSECP256K1 at Jacobian Form  Copyright (C) 2021  "
              "Gustavo Madureira" + "\n"
              "\t\tThis program comes with ABSOLUTELY NO WARRANTY." + "\n"
              "\t\tThis is free software, and you are welcome to "
              "redistribute it" + "\n"
              "\t\tunder certain conditions." +
              "\033[0m")
        print("\n\033[92m" +
              "           Point Number: " + str(private_key) + "\n"
              "            Private Key: " +
              hex(private_key)[2:].zfill(64).upper() + "\n"
              "Uncompressed Public Key: " + "04" +
              hex(x(public_key))[2:].zfill(64).upper() +
              hex(y(public_key))[2:].zfill(64).upper() + "\n"
              "  Compressed Public Key: " + i +
              hex(x(public_key))[2:].zfill(64).upper() + "\n" +
              "\033[0m")
        private_key += 1
