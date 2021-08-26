"""
[ LICENSE ]  ( GNU GPLv3 )  :

SECP256K1 at Jacobian Form  -  Elliptic Curve Cryptography (ECC).
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
the Jacobian form.

- Used to perform public key generation;
- Used to perform digital signature generation and verification with
  ECDSA;
- Used to perform data encryption and decryption.

(Curve used in cryptocurrencies such as Bitcoin, Ethereum, etc...)

    Source:

            https://github.com/gtmadureira/cryptography/blob/main/ecc/secp256k1/secp256k1_jacobian.py

    Author:

            • Gustavo Madureira (gtmadureira@gmail.com)
            • https://gtmadureira.github.io/
"""


from typing import Final, Tuple


# Type Hints.
Point = Tuple[int, int]
Jacobian_Coordinate = Tuple[int, int, int]


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


# The point that points to infinity on the elliptic curve will be
# defined by:
POINT_INFINITY_CURVE: Final[Point] = (0, 0)


# The point that points to infinity on the elliptic curve over Jacobian
# coordinate will be defined by:
POINT_INFINITY_JACOBIAN: Final[Jacobian_Coordinate] = (1, 1, 0)


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
    result = point is POINT_INFINITY_CURVE or 0 in point
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


def is_infinite_jacobian(point: Jacobian_Coordinate) -> bool:
    """
    Returns True if the point at infinity on the elliptic curve over
    Jacobian coordinate, otherwise it returns False.
    """
    _, _, zp = point
    result = point is POINT_INFINITY_JACOBIAN or zp == 0
    return result


def is_on_curve_jacobian(point: Jacobian_Coordinate) -> bool:
    """
    Returns True if the point lies on the elliptic curve over Jacobian
    coordinate, otherwise it returns False.
    """
    if is_infinite_jacobian(point):
        result = True
        return result
    xp, yp, zp = point
    zp_2 = pow(zp, 2, FP_CURVE)
    zp_4 = pow(zp, 4, FP_CURVE)
    result = (pow(yp, 2, FP_CURVE) - pow(xp, 3, FP_CURVE) -
              A_CURVE * xp * zp_4 - B_CURVE * zp_2 * zp_4) % FP_CURVE == 0
    return result


def is_affine_jacobian(point: Jacobian_Coordinate) -> bool:
    """
    Returns True if the point is affine form in Jacobian coordinate
    (x, y, 1).
    """
    assert not is_infinite_jacobian(point)
    _, _, zp = point
    result = zp == 1
    return result


def to_jacobian(point: Point) -> Jacobian_Coordinate:
    """
    Convert an affine point to Jacobian coordinate, or returns point at
    infinity.

    A Jacobian coordinate is represented as (x, y, z).
    """
    assert is_on_curve(point)
    if is_infinite(point):
        result = POINT_INFINITY_JACOBIAN
        return result
    xp, yp = point
    xp, yp, zp = (xp, yp, 1)
    result = (xp, yp, zp)
    return result


def from_jacobian(point: Jacobian_Coordinate) -> Point:
    """
    Convert a Jacobian coordinate to affine point, or returns point at
    infinity.

    An affine point is represented as (x, y).
    """
    assert is_on_curve_jacobian(point)
    if is_infinite_jacobian(point):
        result = POINT_INFINITY_CURVE
        return result
    xp, yp, zp = point
    if is_affine_jacobian(point):
        result = (xp, yp)
        return result
    zp_inv = modular_inverse(zp, FP_CURVE)
    zp_inv_2 = pow(zp_inv, 2, FP_CURVE)
    zp_inv_3 = pow(zp_inv, 3, FP_CURVE)
    result = ((zp_inv_2 * xp) % FP_CURVE, (zp_inv_3 * yp) % FP_CURVE)
    return result


def jacobian_point_doubling(
        point_p: Jacobian_Coordinate) -> Jacobian_Coordinate:
    """
    Point doubling on the elliptic curve over Jacobian coordinate
    (x, y, z).

    It doubles Point-P.
    """
    assert is_on_curve_jacobian(point_p)
    if is_infinite_jacobian(point_p):
        result = POINT_INFINITY_JACOBIAN
        return result
    xp, yp, zp = point_p
    xp_2 = pow(xp, 2, FP_CURVE)
    yp_2 = pow(yp, 2, FP_CURVE)
    yp_4 = pow(yp, 4, FP_CURVE)
    s = (4 * xp * yp_2) % FP_CURVE
    m = (3 * xp_2) % FP_CURVE
    if A_CURVE:
        zp_4 = pow(zp, 4, FP_CURVE)
        m += (A_CURVE * zp_4) % FP_CURVE
    m = m % FP_CURVE
    xr = (pow(m, 2, FP_CURVE) - 2 * s) % FP_CURVE
    yr = (m * (s - xr) - 8 * yp_4) % FP_CURVE
    zr = (2 * yp * zp) % FP_CURVE
    result = (xr, yr, zr)
    assert is_on_curve_jacobian(result)
    return result


def jacobian_point_addition_mixed(
        point_p: Jacobian_Coordinate,
        point_q: Jacobian_Coordinate) -> Jacobian_Coordinate:
    """
    Point addition (mixed) on the elliptic curve over Jacobian
    coordinate (x, y, z) and affine form in Jacobian coordinate
    (x, y, 1).

    It adds Point-P with Point-Q.
    """
    assert is_on_curve_jacobian(point_p)
    assert is_on_curve_jacobian(point_q)
    assert is_affine_jacobian(point_q)
    if is_infinite_jacobian(point_p):
        result = point_q
        return result
    xp, yp, zp = point_p
    xq, yq, _ = point_q
    zp_2 = pow(zp, 2, FP_CURVE)
    zp_3 = pow(zp, 3, FP_CURVE)
    uq = (xq * zp_2) % FP_CURVE
    sq = (yq * zp_3) % FP_CURVE
    if xp == uq and yp != sq:
        result = POINT_INFINITY_JACOBIAN
        return result
    if xp == uq and yp == sq:
        result = jacobian_point_doubling(point_p)
        return result
    h = (uq - xp) % FP_CURVE
    r = (sq - yp) % FP_CURVE
    h_2 = pow(h, 2, FP_CURVE)
    h_3 = pow(h, 3, FP_CURVE)
    xp_h_2 = (xp * h_2) % FP_CURVE
    xr = (pow(r, 2, FP_CURVE) - h_3 - 2 * xp_h_2) % FP_CURVE
    yr = (r * (xp_h_2 - xr) - yp * h_3) % FP_CURVE
    zr = (h * zp) % FP_CURVE
    result = (xr, yr, zr)
    assert is_on_curve_jacobian(result)
    return result


def jacobian_point_addition(
        point_p: Jacobian_Coordinate,
        point_q: Jacobian_Coordinate) -> Jacobian_Coordinate:
    """
    Point addition on the elliptic curve over Jacobian coordinate
    (x, y, z).

    It adds Point-P with Point-Q.
    """
    assert is_on_curve_jacobian(point_p)
    assert is_on_curve_jacobian(point_q)
    if is_infinite_jacobian(point_p):
        result = point_q
        return result
    if is_infinite_jacobian(point_q):
        result = point_p
        return result
    if is_affine_jacobian(point_p):
        result = jacobian_point_addition_mixed(point_q, point_p)
        return result
    if is_affine_jacobian(point_q):
        result = jacobian_point_addition_mixed(point_p, point_q)
        return result
    xp, yp, zp = point_p
    xq, yq, zq = point_q
    zp_2 = pow(zp, 2, FP_CURVE)
    zp_3 = pow(zp, 3, FP_CURVE)
    zq_2 = pow(zq, 2, FP_CURVE)
    zq_3 = pow(zq, 3, FP_CURVE)
    up = (xp * zq_2) % FP_CURVE
    uq = (xq * zp_2) % FP_CURVE
    sp = (yp * zq_3) % FP_CURVE
    sq = (yq * zp_3) % FP_CURVE
    if up == uq and sp != sq:
        result = POINT_INFINITY_JACOBIAN
        return result
    if up == uq and sp == sq:
        result = jacobian_point_doubling(point_p)
        return result
    h = (uq - up) % FP_CURVE
    r = (sq - sp) % FP_CURVE
    h_2 = pow(h, 2, FP_CURVE)
    h_3 = pow(h, 3, FP_CURVE)
    up_h_2 = (up * h_2) % FP_CURVE
    xr = (pow(r, 2, FP_CURVE) - h_3 - 2 * up_h_2) % FP_CURVE
    yr = (r * (up_h_2 - xr) - sp * h_3) % FP_CURVE
    zr = (h * zp * zq) % FP_CURVE
    result = (xr, yr, zr)
    assert is_on_curve_jacobian(result)
    return result


def fast_jacobian_point_addition(point_p: Point, point_q: Point) -> Point:
    """
    Fast point addition on the elliptic curve over affine (x, y) to
    Jacobian coordinate (x, y, z).

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
    jacobian_p = to_jacobian(point_p)
    jacobian_q = to_jacobian(point_q)
    jacobian_r = jacobian_point_addition_mixed(jacobian_p, jacobian_q)
    point_r = from_jacobian(jacobian_r)
    result = point_r
    assert is_on_curve(result)
    return result


def fast_jacobian_point_multiplication(scalar: int, point: Point) -> Point:
    """
    Fast scalar multiplication of point on the elliptic curve over
    affine (x, y) to Jacobian coordinate (x, y, z).

    It doubles Point-P and adds Point-P with Point-Q.
    """
    assert is_on_curve(point)
    if scalar == 0 or is_infinite(point):
        result = POINT_INFINITY_CURVE
        return result
    if scalar < 0 or scalar >= N_CURVE:
        result = fast_jacobian_point_multiplication(scalar % N_CURVE, point)
        return result
    scalar_binary = bin(scalar)[2:]
    jacobian = to_jacobian(point)
    current = jacobian
    for i in range(1, len(scalar_binary)):
        current = jacobian_point_doubling(current)
        if scalar_binary[i] == "1":
            current = jacobian_point_addition(jacobian, current)
    affine = from_jacobian(current)
    result = affine
    assert is_on_curve(result)
    return result


if __name__ == "__main__":

    # Elliptic curve scalar multiplication test.
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

    private_key = \
        0xE05AF5BC208C749190567B921A0C28FE112CD8B54E9FF82F77FA58998B694D4C
    limit = private_key + 10001
    while private_key < limit:
        public_key = fast_jacobian_point_multiplication(
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
        private_key += 1
        sleep(0.0)
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
    elapsed = perf_counter() - start
    print(f"""\033[92m
        Finished in {elapsed:.02f} seconds.
\033[0m""")
