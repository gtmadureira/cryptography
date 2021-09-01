"""
[ LICENSE ]  ( GNU GPLv3 )  :

SECP256K1 at Weierstrass Form  -  Elliptic Curve Cryptography (ECC).
Copyright (C) 2021  Gustavo Madureira

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see https://www.gnu.org/licenses/.


[ ABOUT ]  :

Elliptic Curve Cryptography (ECC).
Module for asymmetric cryptography with elliptic curve SECP256K1, using
the Weierstrass form.

- Used to perform public key generation;
- Used to perform digital signature generation and verification with
  ECDSA and Schnorr;
- Used to perform data encryption and decryption.

(Curve used in cryptocurrencies such as Bitcoin, Ethereum, etc...)

Works on Python 3.8 or higher.

    Source:

            https://github.com/gtmadureira/cryptography/blob/main/ecc/secp256k1/secp256k1_weierstrass.py

    Author:

            • Gustavo Madureira (gtmadureira@gmail.com)
            • https://gtmadureira.github.io/
"""


from typing import Final, Tuple


# Type Hints.
Point = Tuple[int, int]


#       Mathematical domain parameters of the elliptic curve SECP256K1.
#       Source: https://www.secg.org/sec2-v2.pdf


# The finite field (Fp) is defined by:
FP_CURVE: Final[int] = \
    0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F


# The elliptic curve (y^2 = x^3 + ax + b) over Fp is defined by:
A_CURVE: Final[int] = \
    0x0000000000000000000000000000000000000000000000000000000000000000

B_CURVE: Final[int] = \
    0x0000000000000000000000000000000000000000000000000000000000000007


# The generator point is defined by:
GENERATOR_POINT_CURVE: Final[Point] = \
    (0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
     0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8)


# The order of generator point and the cofactor are defined by:
N_CURVE: Final[int] = \
    0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

H_CURVE: Final[int] = \
    0x0000000000000000000000000000000000000000000000000000000000000001


# The point that points to infinity on the elliptic curve is defined by:
POINT_INFINITY_CURVE: Final[Point] = (0, 0)


def modular_inverse(k: int, p: int) -> int:
    """
    Extended Euclidean algorithm/division on the elliptic curve.

    Returns the multiplicative inverse of {k % p}. Where the only
    integer {x} is defined such that {(k * x) % p == 1}.

    {k} must be non-zero and {p} must be a prime.
    """
    if k == 0:
        raise ZeroDivisionError("Division by zero!")
    if k < 0:
        result = p - modular_inverse(- k, p)
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


def is_infinite(point: Point) -> bool:
    """
    Returns True if the point at infinity on the elliptic curve,
    otherwise it returns False.
    """
    result = point == POINT_INFINITY_CURVE or 0 in point
    return result


def is_on_curve(point: Point) -> bool:
    """
    Returns True if the point lies on the elliptic curve, otherwise it
    returns False.
    """
    if is_infinite(point):
        result = True
        return result
    xp, yp = point
    result = (pow(yp, 2, FP_CURVE) - pow(xp, 3, FP_CURVE) -
              A_CURVE * xp - B_CURVE) % FP_CURVE == 0
    return result


def x(point: Point) -> int:
    """
    Refers to {x} coordinate of point, assuming it is not at infinity,
    then returns {x}.
    """
    assert not is_infinite(point)
    assert is_on_curve(point)
    xp, _ = point
    result = xp
    return result


def y(point: Point) -> int:
    """
    Refers to {y} coordinate of point, assuming it is not at infinity,
    then returns {y}.
    """
    assert not is_infinite(point)
    assert is_on_curve(point)
    _, yp = point
    result = yp
    return result


def has_even_y(point: Point) -> bool:
    """
    Where point is not at infinity, it returns True if {yp mod 2 = 0},
    otherwise it returns False.
    """
    assert not is_infinite(point)
    assert is_on_curve(point)
    _, yp = point
    result = yp % 2 == 0
    return result


def ec_point_doubling(point_p: Point) -> Point:
    """
    Point doubling on the elliptic curve.

    It doubles Point-P.
    """
    assert is_on_curve(point_p)
    if is_infinite(point_p):
        result = POINT_INFINITY_CURVE
        return result
    xp, yp = point_p
    slope = ((3 * pow(xp, 2, FP_CURVE) + A_CURVE) *
             modular_inverse(2 * yp, FP_CURVE)) % FP_CURVE
    xr = (pow(slope, 2, FP_CURVE) - 2 * xp) % FP_CURVE
    yr = (slope * (xp - xr) - yp) % FP_CURVE
    result = (xr, yr)
    assert is_on_curve(result)
    return result


def ec_point_addition(point_p: Point, point_q: Point) -> Point:
    """
    Point addition on the elliptic curve.

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
    xp, yp = point_p
    xq, yq = point_q
    if xp == xq and yp != yq:
        result = POINT_INFINITY_CURVE
        return result
    if xp == xq and yp == yq:
        result = ec_point_doubling(point_p)
        return result
    slope = ((yq - yp) * modular_inverse(xq - xp, FP_CURVE)) % FP_CURVE
    xr = (pow(slope, 2, FP_CURVE) - xp - xq) % FP_CURVE
    yr = (slope * (xp - xr) - yp) % FP_CURVE
    result = (xr, yr)
    assert is_on_curve(result)
    return result


def ec_point_multiplication(scalar: int, point: Point) -> Point:
    """
    Scalar multiplication of point on the elliptic curve.

    It doubles Point-P and adds Point-P with Point-Q.
    """
    assert is_on_curve(point)
    if scalar == 0 or is_infinite(point):
        result = POINT_INFINITY_CURVE
        return result
    if scalar < 0 or scalar >= N_CURVE:
        result = ec_point_multiplication(scalar % N_CURVE, point)
        return result
    scalar_binary = bin(scalar)[2:]
    current = point
    for i in range(1, len(scalar_binary)):
        current = ec_point_doubling(current)
        if scalar_binary[i] == "1":
            current = ec_point_addition(point, current)
    result = current
    assert is_on_curve(result)
    return result


if __name__ == "__main__":

    # Elliptic curve scalar multiplication test.
    from platform import system
    from random import randrange
    from time import perf_counter
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

    start = perf_counter()
    elapsed = 0.0
    while elapsed < 900.0:
        private_key = randrange(1, N_CURVE)
        public_key = ec_point_multiplication(
            private_key, GENERATOR_POINT_CURVE)
        if has_even_y(public_key):
            prefix = "02"
        else:
            prefix = "03"
        data = (str(private_key),
                hex(private_key)[2:].zfill(64).upper(),
                hex(x(public_key))[2:].zfill(64).upper(),
                hex(y(public_key))[2:].zfill(64).upper(),
                prefix)
        clear()
        elapsed = perf_counter() - start
        print(f"""\033[92m
        SECP256K1 at Weierstrass Form  Copyright (C) 2021  Gustavo Madureira
        This program comes with ABSOLUTELY NO WARRANTY.
        This is free software, and you are welcome to redistribute it
        under certain conditions.


           Point Number: {data[0]}
            Private Key: {data[1]}
Uncompressed Public Key: 04{data[2]}{data[3]}
  Compressed Public Key: {data[4]}{data[2]}

           Elapsed time: {elapsed:.02f} seconds.
\033[0m""")
