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


    Author Nickname:

            • Hattoshi Hanzōmoto (hanzomoto_hattoshi@protonmail.com)
"""


from typing import Final, Tuple


# Type Hints.
Point = Tuple[int, int]


#       Mathematical domain parameters of the elliptic curve SECP256K1.
#       Source: https://www.secg.org/sec2-v2.pdf


# The finite field (Fp) is defined by:
FP_CURVE: Final[int] = \
    0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F


# The elliptic curve (y^2 = x^3 + ax + b) over (Fp) is defined by the
# coefficients:
A_CURVE: Final[int] = \
    0x0000000000000000000000000000000000000000000000000000000000000000

B_CURVE: Final[int] = \
    0x0000000000000000000000000000000000000000000000000000000000000007


# The generator point (G) is defined by:
GENERATOR_POINT_CURVE: Final[Point] = \
    (0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
     0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8)


# The order (n) of generator point and the cofactor (h) are defined by:
N_CURVE: Final[int] = \
    0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

H_CURVE: Final[int] = \
    0x0000000000000000000000000000000000000000000000000000000000000001


# The point that points to infinity on the elliptic curve is defined by:
POINT_INFINITY_CURVE: Final[Point] = (0, 0)


def modular_inverse(_k: int, _p: int) -> int:
    """
    Using Fermat's little theorem on the elliptic curve to find the
    modular multiplicative inverse.

    Returns the modular multiplicative inverse of {_k % _p}. Where the
    only integer {_x} is defined such that {(_k * _x) % _p == 1}.

    {_p} must be a prime number.
    """
    if _k % _p == 0:
        result = None
        return result  # type: ignore
    _x = pow(_k, _p - 2, _p)
    assert (_k * _x) % _p == 1
    result = _x
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
    _xp, _yp = point
    result = (pow(_yp, 2, FP_CURVE) - pow(_xp, 3, FP_CURVE) -
              A_CURVE * _xp - B_CURVE) % FP_CURVE == 0
    return result


def x_coordinate(point: Point) -> int:
    """
    Refers to {x} coordinate of point, assuming it is not at infinity,
    then returns {x}.
    """
    assert not is_infinite(point)
    assert is_on_curve(point)
    _xp, _ = point
    result = _xp
    return result


def y_coordinate(point: Point) -> int:
    """
    Refers to {y} coordinate of point, assuming it is not at infinity,
    then returns {y}.
    """
    assert not is_infinite(point)
    assert is_on_curve(point)
    _, _yp = point
    result = _yp
    return result


def has_even_y(point: Point) -> bool:
    """
    Where point is not at infinity, it returns True if {_yp % 2 == 0},
    otherwise it returns False.
    """
    assert not is_infinite(point)
    assert is_on_curve(point)
    _yp = y_coordinate(point)
    result = _yp % 2 == 0
    return result


def ec_point_doubling(point_p: Point) -> Point:
    """
    Point doubling on the elliptic curve.

    It doubles Point-P.
    """
    assert is_on_curve(point_p)
    if is_infinite(point_p):
        point_r = POINT_INFINITY_CURVE
        result = point_r
        assert is_on_curve(result)
        return result
    _xp, _yp = point_p
    slope = ((3 * pow(_xp, 2, FP_CURVE) + A_CURVE) *
             modular_inverse(2 * _yp, FP_CURVE)) % FP_CURVE
    _xr = (pow(slope, 2, FP_CURVE) - 2 * _xp) % FP_CURVE
    _yr = (slope * (_xp - _xr) - _yp) % FP_CURVE
    point_r = (_xr, _yr)
    result = point_r
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
        point_r = point_q
        result = point_r
        return result
    if is_infinite(point_q):
        point_r = point_p
        result = point_r
        return result
    _xp, _yp = point_p
    _xq, _yq = point_q
    if _xp == _xq and _yp != _yq:
        point_r = POINT_INFINITY_CURVE
        result = point_r
        assert is_on_curve(result)
        return result
    if _xp == _xq and _yp == _yq:
        point_r = ec_point_doubling(point_p)
        result = point_r
        return result
    slope = ((_yq - _yp) * modular_inverse(_xq - _xp, FP_CURVE)) % FP_CURVE
    _xr = (pow(slope, 2, FP_CURVE) - _xp - _xq) % FP_CURVE
    _yr = (slope * (_xp - _xr) - _yp) % FP_CURVE
    point_r = (_xr, _yr)
    result = point_r
    assert is_on_curve(result)
    return result


def ec_point_multiplication(scalar: int, point: Point) -> Point:
    """
    Scalar multiplication of point on the elliptic curve.

    It doubles Point-P and adds Point-P with Point-Q.

    - Point-P is defined as Point and Point-Q is defined as
      Current.
    """
    assert is_on_curve(point)
    if scalar == 0 or scalar == N_CURVE or is_infinite(point):
        new_point = POINT_INFINITY_CURVE
        result = new_point
        assert is_on_curve(result)
        return result
    if scalar == 1:
        result = point
        return result
    if scalar < 0 or scalar > N_CURVE:
        new_point = ec_point_multiplication(scalar % N_CURVE, point)
        result = new_point
        return result
    scalar_binary = bin(scalar)[2:]
    current = point
    for i in range(1, len(scalar_binary)):
        current = ec_point_doubling(current)
        if scalar_binary[i] == "1":
            current = ec_point_addition(point, current)
    new_point = current
    result = new_point
    assert is_on_curve(result)
    return result


if __name__ == "__main__":

    # Elliptic curve scalar multiplication test.
    from random import randrange
    from time import perf_counter
    from colorama import just_fix_windows_console  # type: ignore
    from os import name as system_type, system as run_command

    # Define clear_screen function.
    def clear_screen() -> None:
        """
        Tests the operating system type and sets the screen clear command.
        """
        # Screen clear command for Windows operating system.
        if system_type == "nt":
            _ = run_command("cls")

        # Screen clear command for macOS/Linux operating system.
        elif system_type == "posix":
            _ = run_command("clear")

    # Get ANSI escapes from color scheme to work on Windows.
    just_fix_windows_console()

    start = perf_counter()
    ELAPSED = 0.0
    while ELAPSED < 900.0:
        private_key = randrange(1, N_CURVE)
        public_key = ec_point_multiplication(
            private_key, GENERATOR_POINT_CURVE)
        if has_even_y(public_key):
            PREFIX = "02"
        else:
            PREFIX = "03"
        data = (str(private_key),
                hex(private_key)[2:].zfill(64).upper(),
                hex(x_coordinate(public_key))[2:].zfill(64).upper(),
                hex(y_coordinate(public_key))[2:].zfill(64).upper(),
                PREFIX)
        clear_screen()
        ELAPSED = perf_counter() - start
        print(f"""\033[92m
         SECP256K1 at Weierstrass Form  Copyright (C) 2021  Gustavo Madureira
         License GNU GPL-3.0-or-later <https://gnu.org/licenses/gpl.html>
         This program comes with ABSOLUTELY NO WARRANTY.
         This is free software, and you are welcome to redistribute it
         under certain conditions.

         \33[1;7m *** FOR EDUCATIONAL PURPOSES ONLY *** \033[0m


            \033[92mPoint Number: {data[0]}
             Private Key: {data[1]}
 Uncompressed Public Key: 04{data[2]}
                            {data[3]}
   Compressed Public Key: {data[4]}{data[2]}

            Elapsed time: {ELAPSED:.02f} seconds.
\033[0m""")
