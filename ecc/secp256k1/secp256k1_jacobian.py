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
along with this program.  If not, see <https://www.gnu.org/licenses/>.


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
from typing import Tuple
from platform import system
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


def to_jacobian(point: Point) -> Jacobian_Coordinate:
    o = (point[0], point[1], 1)
    result = o
    return result


def jacobian_double(point_P: Jacobian_Coordinate) -> Jacobian_Coordinate:
    if not point_P[1]:
        result = (0, 0, 0)
        return result
    ysq = (point_P[1] ** 2) % _FP_CURVE_
    S = (4 * point_P[0] * ysq) % _FP_CURVE_
    M = (3 * point_P[0] ** 2 + _A_CURVE_ * point_P[2] ** 4) % _FP_CURVE_
    nx = (M ** 2 - 2 * S) % _FP_CURVE_
    ny = (M * (S - nx) - 8 * ysq ** 2) % _FP_CURVE_
    nz = (2 * point_P[1] * point_P[2]) % _FP_CURVE_
    result = (nx, ny, nz)
    return result


def jacobian_add(point_P: Jacobian_Coordinate,
                 point_Q: Jacobian_Coordinate) -> Jacobian_Coordinate:
    if not point_P[1]:
        result = point_Q
        return result
    if not point_Q[1]:
        result = point_P
        return result
    U1 = (point_P[0] * point_Q[2] ** 2) % _FP_CURVE_
    U2 = (point_Q[0] * point_P[2] ** 2) % _FP_CURVE_
    S1 = (point_P[1] * point_Q[2] ** 3) % _FP_CURVE_
    S2 = (point_Q[1] * point_P[2] ** 3) % _FP_CURVE_
    if U1 == U2:
        if S1 != S2:
            result = (0, 0, 1)
            return result
        result = jacobian_double(point_P)
        return result
    H = U2 - U1
    R = S2 - S1
    H2 = (H * H) % _FP_CURVE_
    H3 = (H * H2) % _FP_CURVE_
    U1H2 = (U1 * H2) % _FP_CURVE_
    nx = (R ** 2 - H3 - 2 * U1H2) % _FP_CURVE_
    ny = (R * (U1H2 - nx) - S1 * H3) % _FP_CURVE_
    nz = (H * point_P[2] * point_Q[2]) % _FP_CURVE_
    result = (nx, ny, nz)
    return result


def from_jacobian(point: Jacobian_Coordinate) -> Point:
    z = modular_inverse(point[2], _FP_CURVE_)
    result = ((point[0] * z ** 2) % _FP_CURVE_,
              (point[1] * z ** 3) % _FP_CURVE_)
    return result


def jacobian_multiply(point: Jacobian_Coordinate,
                      scalar: int) -> Jacobian_Coordinate:
    result = None
    if point[1] == 0 or scalar == 0:
        result = (0, 0, 1)
        return result
    if scalar == 1:
        result = point
        return result
    if scalar < 0 or scalar >= _N_CURVE_:
        result = jacobian_multiply(point, scalar % _N_CURVE_)
        return result
    if (scalar % 2) == 0:
        result = jacobian_double(jacobian_multiply(point, scalar // 2))
        return result
    if (scalar % 2) == 1:
        result = jacobian_add(jacobian_double(
            jacobian_multiply(point, scalar // 2)), point)
        return result
    return result  # type: ignore


def ec_point_multiplication(scalar: int,
                            point: Point = _GENERATOR_POINT_CURVE) -> Point:
    result = from_jacobian(jacobian_multiply(to_jacobian(point), scalar))
    return result


if __name__ == "__main__":
    private_key = 1
    while True:
        public_key = ec_point_multiplication(private_key)
        if public_key[1] % 2 == 1:
            i = "03"
        else:
            i = "02"
        sleep(0.0)
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
              hex(public_key[0])[2:].zfill(64).upper() +
              hex(public_key[1])[2:].zfill(64).upper() + "\n"
              "  Compressed Public Key: " + i +
              hex(public_key[0])[2:].zfill(64).upper() + "\n" +
              "\033[0m")
        private_key += 1
