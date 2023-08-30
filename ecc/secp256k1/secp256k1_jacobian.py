"""
[ LICENSE ]  ( GNU GPLv3 )  :

SECP256K1 at Jacobian Coordinates  -  Elliptic Curve Cryptography (ECC).
Copyright (C)  2021 - 2023  Gustavo Madureira

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
the Jacobian coordinates.

- Used to perform public key generation;
- Used to perform digital signature generation and verification with
  ECDSA and Schnorr;
- Used to perform data encryption and decryption.

(Curve used in cryptocurrencies such as Bitcoin, Ethereum, etc...)

For educational purposes only.

Works on Python 3.8 or higher.


    Source:

            https://github.com/gtmadureira/cryptography/blob/main/ecc/secp256k1/secp256k1_jacobian.py


    Author:

            • Gustavo Madureira (gtmadureira@gmail.com)
            • https://gtmadureira.github.io/


    Author Nickname:

            • Hattoshi Hanzōmoto (hanzomoto_hattoshi@protonmail.com)
"""


from types import NoneType
from typing import Final, Tuple

# Type Hints.
Point = Tuple[int, int]
JacobianCoordinate = Tuple[int, int, int]
ScalarBinary = Tuple[str, ...]


# *****************************************************************************
# *                                                                           *
# *       Manipulation of Binary, Hexadecimal and Integer values.             *
# *                                                                           *
# *****************************************************************************
# *                                                                           *
# *****************************************************************************
# ********************************* [ BEGIN ] *********************************
# *****************************************************************************


def is_bin_encoded(_bin: str) -> bool:
    """Checks if the input string is an encoded binary value."""
    assert isinstance(_bin, str)
    try:
        result = bool(((int(_bin, 2) == 0) or int(_bin, 2)) and
                      ((len(_bin) % 2) == 0))
    except ValueError:
        result = False
    assert isinstance(result, bool)
    return result


def is_hex_encoded(_hex: str) -> bool:
    """Checks if the input string is an encoded hexadecimal value."""
    assert isinstance(_hex, str)
    try:
        result = bool(((int(_hex, 16) == 0) or int(_hex, 16)) and
                      ((len(_hex) % 2) == 0))
    except ValueError:
        result = False
    assert isinstance(result, bool)
    return result


def hex_from_int(_x: int, output_length_bytes: int) -> str:
    """Converts an integer to encoded hexadecimal string."""
    assert isinstance(_x and output_length_bytes, int)
    result = hex(_x)[2:].zfill(output_length_bytes * 2)
    assert (isinstance(result, str) and
            is_hex_encoded(result))
    return result


# *****************************************************************************
# ********************************** [ END ] **********************************
# *****************************************************************************


# *****************************************************************************
# *                                                                           *
# *      Mathematical Domain Parameters of the Elliptic Curve SECP256K1,      *
# *      over Finite Field (Fp).                                              *
# *                                                                           *
# *      Source:   https://www.secg.org/sec2-v2.pdf                           *
# *                                                                           *
# *****************************************************************************
# *                                                                           *
# *****************************************************************************
# ********************************* [ BEGIN ] *********************************
# *****************************************************************************


# The Finite Field (Fp) is defined by:
FP_CURVE: Final[int] = \
    0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F

assert (isinstance(FP_CURVE, int) and
        is_hex_encoded(hex_from_int(FP_CURVE, 32)))


# The Elliptic Curve in short Weierstrass form ( y² = x³ + a⋅x + b )
# over (Fp) is defined by the coefficients:
A_CURVE: Final[int] = \
    0x0000000000000000000000000000000000000000000000000000000000000000

B_CURVE: Final[int] = \
    0x0000000000000000000000000000000000000000000000000000000000000007

assert (isinstance(A_CURVE and B_CURVE, int) and
        is_hex_encoded(hex_from_int(A_CURVE, 32)) and
        is_hex_encoded(hex_from_int(B_CURVE, 32)))


# The Generator point (G) is defined by:
X_COORD_GENERATOR_POINT_CURVE: Final[int] = \
    0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798

Y_COORD_GENERATOR_POINT_CURVE: Final[int] = \
    0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8

assert (isinstance(X_COORD_GENERATOR_POINT_CURVE and
                   Y_COORD_GENERATOR_POINT_CURVE, int) and
        is_hex_encoded(hex_from_int(X_COORD_GENERATOR_POINT_CURVE, 32)) and
        is_hex_encoded(hex_from_int(Y_COORD_GENERATOR_POINT_CURVE, 32)))

GENERATOR_POINT_CURVE: Final[Point] = (X_COORD_GENERATOR_POINT_CURVE,
                                       Y_COORD_GENERATOR_POINT_CURVE)

assert (isinstance(GENERATOR_POINT_CURVE, tuple) and
        is_hex_encoded(hex_from_int(GENERATOR_POINT_CURVE[0], 32)) and
        is_hex_encoded(hex_from_int(GENERATOR_POINT_CURVE[1], 32)))


# The order (n) of Generator point and the Cofactor (h) are defined by:
N_CURVE: Final[int] = \
    0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

H_CURVE: Final[int] = \
    0x0000000000000000000000000000000000000000000000000000000000000001

assert (isinstance(N_CURVE and H_CURVE, int) and
        is_hex_encoded(hex_from_int(N_CURVE, 32)) and
        is_hex_encoded(hex_from_int(H_CURVE, 32)))


# The point that points to infinity on the Elliptic Curve is defined by:
X_COORD_POINT_INFINITY_CURVE: Final[int] = 0

Y_COORD_POINT_INFINITY_CURVE: Final[int] = 0 if (B_CURVE != 0) else 1

assert isinstance(X_COORD_POINT_INFINITY_CURVE and
                  Y_COORD_POINT_INFINITY_CURVE, int)

POINT_INFINITY_CURVE: Final[Point] = (X_COORD_POINT_INFINITY_CURVE,
                                      Y_COORD_POINT_INFINITY_CURVE)

assert isinstance(POINT_INFINITY_CURVE, tuple)


# The point that points to infinity on the Elliptic Curve over Jacobian
# coordinates is defined by:
X_COORD_POINT_INFINITY_JACOBIAN: Final[int] = 1

Y_COORD_POINT_INFINITY_JACOBIAN: Final[int] = 1

Z_COORD_POINT_INFINITY_JACOBIAN: Final[int] = 0

assert isinstance(X_COORD_POINT_INFINITY_JACOBIAN and
                  Y_COORD_POINT_INFINITY_JACOBIAN and
                  Z_COORD_POINT_INFINITY_JACOBIAN, int)

POINT_INFINITY_JACOBIAN: Final[JacobianCoordinate] = \
    (X_COORD_POINT_INFINITY_JACOBIAN,
     Y_COORD_POINT_INFINITY_JACOBIAN,
     Z_COORD_POINT_INFINITY_JACOBIAN)

assert isinstance(POINT_INFINITY_JACOBIAN, tuple)


# *****************************************************************************
# ********************************** [ END ] **********************************
# *****************************************************************************


# *****************************************************************************
# *                                                                           *
# *        Elliptic Curve Arithmetic, using the Jacobian Coordinates,         *
# *        over a Finite Field (Fp).                                          *
# *                                                                           *
# *        Source:   https://www.hyperelliptic.org/EFD/                       *
# *                                                                           *
# *                  https://www.hyperelliptic.org/EFD/precomp.pdf            *
# *                                                                           *
# *                  https://www.hyperelliptic.org/HEHCC/                     *
# *                                                                           *
# *****************************************************************************
# *                                                                           *
# *****************************************************************************
# ********************************* [ BEGIN ] *********************************
# *****************************************************************************


def modular_inverse(_k: int, _p: int) -> int:
    """
    Using Fermat's little theorem on the Elliptic Curve to find the
    modular multiplicative inverse.

    Returns the modular multiplicative inverse of {_k mod _p}. Where the
    only integer {_x} is defined such that {(_k ⋅ _x) mod _p = 1}.

    {_p} must be a prime number.
    """
    assert isinstance(_k and _p, int)
    if (_k % _p) == 0:
        result = None
        assert isinstance(result, NoneType)
        return result  # type: ignore
    _x = pow(_k, - 1, _p)
    assert bool(((_k * _x) % _p) == 1)
    result = _x
    assert isinstance(result, int)
    return result


def is_infinite(point: Point) -> bool:
    """
    Returns True if the point is at infinity on the Elliptic Curve,
    otherwise it returns False.
    """
    assert isinstance(point, tuple)
    _xp, _yp = point
    result = bool((point == POINT_INFINITY_CURVE) or
                  ((_xp == 0) and (_yp != 0)))
    assert isinstance(result, bool)
    return result


def is_on_curve(point: Point) -> bool:
    """
    Returns True if the point lies on the Elliptic Curve, otherwise it
    returns False.
    """
    assert isinstance(point, tuple)
    if is_infinite(point):
        result = True
        assert isinstance(result, bool)
        return result
    _xp, _yp = point
    result = bool(((pow(_yp, 2, FP_CURVE) - pow(_xp, 3, FP_CURVE) - A_CURVE *
                    _xp - B_CURVE) % FP_CURVE) == 0)
    assert isinstance(result, bool)
    return result


def x_coordinate(point: Point) -> int:
    """
    Refers to {x} coordinate of point, assuming it is not at infinity,
    then returns {x}.
    """
    assert (isinstance(point, tuple) and
            is_on_curve(point) and not
            is_infinite(point))
    _xp, _ = point
    result = _xp
    assert isinstance(result, int)
    return result


def y_coordinate(point: Point) -> int:
    """
    Refers to {y} coordinate of point, assuming it is not at infinity,
    then returns {y}.
    """
    assert (isinstance(point, tuple) and
            is_on_curve(point) and not
            is_infinite(point))
    _, _yp = point
    result = _yp
    assert isinstance(result, int)
    return result


def has_even_y(point: Point) -> bool:
    """
    Where the point is not at infinity, it returns True so that the
    y-coordinate is an even value if {_yp mod 2 = 0}, otherwise it
    returns False.
    """
    assert (isinstance(point, tuple) and
            is_on_curve(point) and not
            is_infinite(point))
    _yp = y_coordinate(point)
    result = bool((_yp % 2) == 0)
    assert isinstance(result, bool)
    return result


def lift_x(_xp: int) -> Point:
    """
    Given an x-coordinate on the curve, return a corresponding affine
    point for which the y-coordinate is even:


    {x(point) = x}

    and

    {y(point) = y} if {y mod 2 = 0} else
    {y(point) = (Fp) - y} otherwise.
    """
    assert isinstance(_xp, int)
    if not 0 < _xp < FP_CURVE:
        result = POINT_INFINITY_CURVE
        assert isinstance(result, tuple)
        return result
    yp_2 = (pow(_xp, 3, FP_CURVE) + A_CURVE * _xp + B_CURVE) % FP_CURVE
    _yp = pow(yp_2, (FP_CURVE + 1) // 4, FP_CURVE)
    if pow(_yp, 2, FP_CURVE) != yp_2:
        result = POINT_INFINITY_CURVE
        assert isinstance(result, tuple)
        return result
    point = (_xp, _yp)
    result = (x_coordinate(point),
              y_coordinate(point) if has_even_y(point) else
              (FP_CURVE - y_coordinate(point)))
    assert isinstance(result, tuple)
    return result


def is_infinite_jacobian(jacobian: JacobianCoordinate) -> bool:
    """
    Returns True if the point is at infinity on the Elliptic Curve over
    Jacobian coordinates, otherwise it returns False.
    """
    assert isinstance(jacobian, tuple)
    _xp, _yp, _zp = jacobian
    result = bool((jacobian == POINT_INFINITY_JACOBIAN) or
                  ((_zp == 0) and (0 not in (_xp, _yp))))
    assert isinstance(result, bool)
    return result


def is_on_curve_jacobian(jacobian: JacobianCoordinate) -> bool:
    """
    Returns True if the point lies on the Elliptic Curve over Jacobian
    coordinate, otherwise it returns False.
    """
    assert isinstance(jacobian, tuple)
    if is_infinite_jacobian(jacobian):
        result = True
        assert isinstance(result, bool)
        return result
    _xp, _yp, _zp = jacobian
    zp_2 = pow(_zp, 2, FP_CURVE)
    zp_4 = pow(_zp, 4, FP_CURVE)
    result = bool(((pow(_yp, 2, FP_CURVE) - pow(_xp, 3, FP_CURVE) - A_CURVE *
                    _xp * zp_4 - B_CURVE * zp_2 * zp_4) % FP_CURVE) == 0)
    assert isinstance(result, bool)
    return result


def is_affine_jacobian(jacobian: JacobianCoordinate) -> bool:
    """
    Returns True if the point is affine form in Jacobian coordinates
    (x, y, 1).
    """
    assert (isinstance(jacobian, tuple) and
            is_on_curve_jacobian(jacobian) and not
            is_infinite_jacobian(jacobian))
    _, _, _zp = jacobian
    result = bool(_zp == 1)
    assert isinstance(result, bool)
    return result


def to_jacobian(point: Point) -> JacobianCoordinate:
    """
    Convert an affine point to Jacobian coordinates, or returns point at
    infinity.

    A Jacobian coordinates is represented as (x, y, z).
    """
    assert (isinstance(point, tuple) and
            is_on_curve(point))
    if is_infinite(point):
        jacobian = POINT_INFINITY_JACOBIAN
        result = jacobian
        assert (isinstance(result, tuple) and
                is_on_curve_jacobian(result))
        return result
    _xp, _yp = point
    _xp, _yp, _zp = (_xp, _yp, 1)  # pylint: disable=W0127
    jacobian = (_xp, _yp, _zp)
    result = jacobian
    assert (isinstance(result, tuple) and
            is_on_curve_jacobian(result))
    return result


def from_jacobian(jacobian: JacobianCoordinate) -> Point:
    """
    Convert a Jacobian coordinates to affine point, or returns point at
    infinity.

    An affine point is represented as (x, y).
    """
    assert (isinstance(jacobian, tuple) and
            is_on_curve_jacobian(jacobian))
    if is_infinite_jacobian(jacobian):
        point = POINT_INFINITY_CURVE
        result = point
        assert (isinstance(result, tuple) and
                is_on_curve(result))
        return result
    _xp, _yp, _zp = jacobian
    if is_affine_jacobian(jacobian):
        point = (_xp, _yp)
        result = point
        assert (isinstance(result, tuple) and
                is_on_curve(result))
        return result
    zp_inv = modular_inverse(_zp, FP_CURVE)
    zp_inv_2 = pow(zp_inv, 2, FP_CURVE)
    zp_inv_3 = pow(zp_inv, 3, FP_CURVE)
    point = ((zp_inv_2 * _xp) % FP_CURVE, (zp_inv_3 * _yp) % FP_CURVE)
    result = point
    assert (isinstance(result, tuple) and
            is_on_curve(result))
    return result


def jacobian_point_doubling(  # pylint: disable=R0914
        jacobian_p: JacobianCoordinate) -> JacobianCoordinate:
    """
    Point doubling on the Elliptic Curve over Jacobian coordinates
    (x, y, z).

    It doubles Point-P.

    - Point-P is defined as Jacobian-P.
    """
    assert (isinstance(jacobian_p, tuple) and
            is_on_curve_jacobian(jacobian_p))
    if is_infinite_jacobian(jacobian_p):
        jacobian_r = POINT_INFINITY_JACOBIAN
        result = jacobian_r
        assert (isinstance(result, tuple) and
                is_on_curve_jacobian(result))
        return result
    _xp, _yp, _zp = jacobian_p
    xp_2 = pow(_xp, 2, FP_CURVE)
    yp_2 = pow(_yp, 2, FP_CURVE)
    yp_4 = pow(yp_2, 2, FP_CURVE)
    zp_2 = pow(_zp, 2, FP_CURVE)
    _s = (2 * (pow(_xp + yp_2, 2, FP_CURVE) - xp_2 - yp_4)) % FP_CURVE
    if not A_CURVE:
        _m = (3 * xp_2) % FP_CURVE
        _t = pow(_m, 2, FP_CURVE)
        _xr = (_t - 2 * _s) % FP_CURVE
        _yr = (_m * (_s - _xr) - 8 * yp_4) % FP_CURVE
        _zr = (2 * _yp * _zp) % FP_CURVE
        jacobian_r = (_xr, _yr, _zr)
        result = jacobian_r
        assert (isinstance(result, tuple) and
                is_on_curve_jacobian(result))
        return result
    if (A_CURVE % FP_CURVE) == (FP_CURVE - 3):
        alpha = (_xp * yp_2) % FP_CURVE
        beta = (3 * (_xp - zp_2) * (_xp + zp_2)) % FP_CURVE
        _xr = (pow(beta, 2, FP_CURVE) - 8 * alpha) % FP_CURVE
        _yr = ((beta * (4 * alpha - _xr) - 8 *
               pow(yp_2, 2, FP_CURVE)) % FP_CURVE)
        _zr = (pow(_yp + _zp, 2, FP_CURVE) - yp_2 - zp_2) % FP_CURVE
        jacobian_r = (_xr, _yr, _zr)
        result = jacobian_r
        assert (isinstance(result, tuple) and
                is_on_curve_jacobian(result))
        return result
    if _zp == 1:
        _m = (3 * xp_2 + A_CURVE) % FP_CURVE
        _zr = (2 * _yp) % FP_CURVE
    else:
        _m = (3 * xp_2 + A_CURVE * pow(zp_2, 2, FP_CURVE)) % FP_CURVE
        _zr = (pow(_yp + _zp, 2, FP_CURVE) - yp_2 - zp_2) % FP_CURVE
    _t = (pow(_m, 2, FP_CURVE) - 2 * _s) % FP_CURVE
    _xr = _t
    _yr = (_m * (_s - _t) - 8 * yp_4) % FP_CURVE
    jacobian_r = (_xr, _yr, _zr)
    result = jacobian_r
    assert (isinstance(result, tuple) and
            is_on_curve_jacobian(result))
    return result


def jacobian_point_addition_affined_only(  # pylint: disable=R0914
        jacobian_p: JacobianCoordinate,
        jacobian_q: JacobianCoordinate) -> JacobianCoordinate:
    """
    Point addition on the Elliptic Curve over Jacobian coordinates
    (x, y, z), with both points with z = 1. That is , only with points
    in affine space (x, y, 1).

    It adds Point-P with Point-Q.

    - Point-P is defined as Jacobian-P and Point-Q is defined as
    Jacobian-Q.
    """
    assert (isinstance(jacobian_p and jacobian_q, tuple) and
            is_on_curve_jacobian(jacobian_p) and
            is_on_curve_jacobian(jacobian_q) and
            is_affine_jacobian(jacobian_p) and
            is_affine_jacobian(jacobian_q))
    _xp, _yp, _ = jacobian_p
    _xq, _yq, _ = jacobian_q
    if (_xp == _xq) and (_yp != _yq):
        jacobian_r = POINT_INFINITY_JACOBIAN
        result = jacobian_r
        assert (isinstance(result, tuple) and
                is_on_curve_jacobian(result))
        return result
    if (_xp == _xq) and (_yp == _yq):
        jacobian_r = jacobian_point_doubling(jacobian_p)
        result = jacobian_r
        assert (isinstance(result, tuple) and
                is_on_curve_jacobian(result))
        return result
    _h = (_xq - _xp) % FP_CURVE
    h_2 = pow(_h, 2, FP_CURVE)
    _i = (4 * h_2) % FP_CURVE
    _j = (_h * _i) % FP_CURVE
    _r = (2 * (_yq - _yp)) % FP_CURVE
    _v = (_xp * _i) % FP_CURVE
    _xr = (pow(_r, 2, FP_CURVE) - _j - 2 * _v) % FP_CURVE
    _yr = (_r * (_v - _xr) - 2 * _yp * _j) % FP_CURVE
    _zr = (2 * _h) % FP_CURVE
    jacobian_r = (_xr, _yr, _zr)
    result = jacobian_r
    assert (isinstance(result, tuple) and
            is_on_curve_jacobian(result))
    return result


def jacobian_point_addition_z_equals(  # pylint: disable=R0914
        jacobian_p: JacobianCoordinate,
        jacobian_q: JacobianCoordinate) -> JacobianCoordinate:
    """
    Point addition on the Elliptic Curve over Jacobian coordinates
    (x, y, z), with both points having z equal (Co-Z), but being != 1.
    That is, both are not in affine space but sharing the same
    z-coordinate (x, y, z != 1).

    It adds Point-P with Point-Q.

    - Point-P is defined as Jacobian-P and Point-Q is defined as
    Jacobian-Q.
    """
    assert (isinstance(jacobian_p and jacobian_q, tuple) and
            is_on_curve_jacobian(jacobian_p) and
            is_on_curve_jacobian(jacobian_q) and not
            is_affine_jacobian(jacobian_p) and not
            is_affine_jacobian(jacobian_q))
    _xp, _yp, _zp = jacobian_p
    _xq, _yq, _ = jacobian_q
    if (_xp == _xq) and (_yp != _yq):
        jacobian_r = POINT_INFINITY_JACOBIAN
        result = jacobian_r
        assert (isinstance(result, tuple) and
                is_on_curve_jacobian(result))
        return result
    if (_xp == _xq) and (_yp == _yq):
        jacobian_r = jacobian_point_doubling(jacobian_p)
        result = jacobian_r
        assert (isinstance(result, tuple) and
                is_on_curve_jacobian(result))
        return result
    _a = pow(_xq - _xp, 2, FP_CURVE)
    _b = (_xp * _a) % FP_CURVE
    _c = (_xq * _a) % FP_CURVE
    _d = pow(_yq - _yp, 2, FP_CURVE)
    _xr = (_d - _b - _c) % FP_CURVE
    _yr = ((_yq - _yp) * (_b - _xr) - _yp * (_c - _b)) % FP_CURVE
    _zr = (_zp * (_xq - _xp)) % FP_CURVE
    jacobian_r = (_xr, _yr, _zr)
    result = jacobian_r
    assert (isinstance(result, tuple) and
            is_on_curve_jacobian(result))
    return result


def jacobian_point_addition_mixed(  # pylint: disable=R0914
        jacobian_p: JacobianCoordinate,
        jacobian_q: JacobianCoordinate) -> JacobianCoordinate:
    """
    Point addition (mixed) on the Elliptic Curve over Jacobian
    coordinate (x, y, z) and affine form in Jacobian coordinates
    (x, y, 1).

    It adds Point-P with Point-Q.

    - Point-P is defined as Jacobian-P and Point-Q is defined as
    Jacobian-Q.
    """
    assert (isinstance(jacobian_p and jacobian_q, tuple) and
            is_on_curve_jacobian(jacobian_p) and
            is_on_curve_jacobian(jacobian_q) and not
            is_affine_jacobian(jacobian_p) and
            is_affine_jacobian(jacobian_q))
    _xp, _yp, _zp = jacobian_p
    _xq, _yq, _ = jacobian_q
    zp_2 = pow(_zp, 2, FP_CURVE)
    _uq = (_xq * zp_2) % FP_CURVE
    _sq = (_yq * _zp * zp_2) % FP_CURVE
    if (_xp == _uq) and (_yp != _sq):
        jacobian_r = POINT_INFINITY_JACOBIAN
        result = jacobian_r
        assert (isinstance(result, tuple) and
                is_on_curve_jacobian(result))
        return result
    if (_xp == _uq) and (_yp == _sq):
        jacobian_r = jacobian_point_doubling(jacobian_p)
        result = jacobian_r
        assert (isinstance(result, tuple) and
                is_on_curve_jacobian(result))
        return result
    _h = (_uq - _xp) % FP_CURVE
    h_2 = pow(_h, 2, FP_CURVE)
    _i = (4 * h_2) % FP_CURVE
    _j = (_h * _i) % FP_CURVE
    _r = (2 * (_sq - _yp)) % FP_CURVE
    _v = (_xp * _i) % FP_CURVE
    _xr = (pow(_r, 2, FP_CURVE) - _j - 2 * _v) % FP_CURVE
    _yr = (_r * (_v - _xr) - 2 * _yp * _j) % FP_CURVE
    _zr = (pow(_zp + _h, 2, FP_CURVE) - zp_2 - h_2) % FP_CURVE
    jacobian_r = (_xr, _yr, _zr)
    result = jacobian_r
    assert (isinstance(result, tuple) and
            is_on_curve_jacobian(result))
    return result


def jacobian_point_addition(  # pylint: disable=R0911, R0914, R0915
        jacobian_p: JacobianCoordinate,
        jacobian_q: JacobianCoordinate) -> JacobianCoordinate:
    """
    Point addition on the Elliptic Curve over Jacobian coordinates
    (x, y, z).

    It adds Point-P with Point-Q.

    - Point-P is defined as Jacobian-P and Point-Q is defined as
    Jacobian-Q.
    """
    assert (isinstance(jacobian_p and jacobian_q, tuple) and
            is_on_curve_jacobian(jacobian_p) and
            is_on_curve_jacobian(jacobian_q))
    if is_infinite_jacobian(jacobian_p) and is_infinite_jacobian(jacobian_q):
        jacobian_r = POINT_INFINITY_JACOBIAN
        result = jacobian_r
        assert (isinstance(result, tuple) and
                is_on_curve_jacobian(result))
        return result
    if is_infinite_jacobian(jacobian_p):
        jacobian_r = jacobian_q
        result = jacobian_r
        assert (isinstance(result, tuple) and
                is_on_curve_jacobian(result))
        return result
    if is_infinite_jacobian(jacobian_q):
        jacobian_r = jacobian_p
        result = jacobian_r
        assert (isinstance(result, tuple) and
                is_on_curve_jacobian(result))
        return result
    if is_affine_jacobian(jacobian_p) and is_affine_jacobian(jacobian_q):
        jacobian_r = jacobian_point_addition_affined_only(jacobian_p,
                                                          jacobian_q)
        result = jacobian_r
        assert (isinstance(result, tuple) and
                is_on_curve_jacobian(result))
        return result
    if is_affine_jacobian(jacobian_p):
        jacobian_r = jacobian_point_addition_mixed(  # pylint: disable=W1114
            jacobian_q, jacobian_p)
        result = jacobian_r
        assert (isinstance(result, tuple) and
                is_on_curve_jacobian(result))
        return result
    if is_affine_jacobian(jacobian_q):
        jacobian_r = jacobian_point_addition_mixed(jacobian_p, jacobian_q)
        result = jacobian_r
        assert (isinstance(result, tuple) and
                is_on_curve_jacobian(result))
        return result
    _xp, _yp, _zp = jacobian_p
    _xq, _yq, _zq = jacobian_q
    if _zp == _zq:
        jacobian_r = jacobian_point_addition_z_equals(jacobian_p, jacobian_q)
        result = jacobian_r
        assert (isinstance(result, tuple) and
                is_on_curve_jacobian(result))
        return result
    zp_2 = pow(_zp, 2, FP_CURVE)
    zq_2 = pow(_zq, 2, FP_CURVE)
    _up = (_xp * zq_2) % FP_CURVE
    _uq = (_xq * zp_2) % FP_CURVE
    _sp = (_yp * _zq * zq_2) % FP_CURVE
    _sq = (_yq * _zp * zp_2) % FP_CURVE
    if (_up == _uq) and (_sp != _sq):
        jacobian_r = POINT_INFINITY_JACOBIAN
        result = jacobian_r
        assert (isinstance(result, tuple) and
                is_on_curve_jacobian(result))
        return result
    if (_up == _uq) and (_sp == _sq):
        jacobian_r = jacobian_point_doubling(jacobian_p)
        result = jacobian_r
        assert (isinstance(result, tuple) and
                is_on_curve_jacobian(result))
        return result
    _h = (_uq - _up) % FP_CURVE
    _i = pow(2 * _h, 2, FP_CURVE)
    _j = (_h * _i) % FP_CURVE
    _r = (2 * (_sq - _sp)) % FP_CURVE
    _v = (_up * _i) % FP_CURVE
    _xr = (pow(_r, 2, FP_CURVE) - _j - 2 * _v) & FP_CURVE
    _yr = (_r * (_v - _xr) - 2 * _sp * _j) % FP_CURVE
    _zr = ((pow(_zp + _zq, 2, FP_CURVE) - zp_2 - zq_2) * _h) % FP_CURVE
    jacobian_r = (_xr, _yr, _zr)
    result = jacobian_r
    assert (isinstance(result, tuple) and
            is_on_curve_jacobian(result))
    return result


def fast_point_addition(point_p: Point, point_q: Point) -> Point:
    """
    Fast point addition on the Elliptic Curve over affine (x, y) to
    Jacobian coordinates (x, y, z).

    It adds Point-P with Point-Q.

    - Point-P is defined as Jacobian-P and Point-Q is defined as
    Jacobian-Q.
    """
    assert (isinstance(point_p and point_q, tuple) and
            is_on_curve(point_p) and
            is_on_curve(point_q))
    if is_infinite(point_p) and is_infinite(point_q):
        point_r = POINT_INFINITY_CURVE
        result = point_r
        assert (isinstance(result, tuple) and
                is_on_curve(result))
        return result
    if is_infinite(point_p):
        point_r = point_q
        result = point_r
        assert (isinstance(result, tuple) and
                is_on_curve(result))
        return result
    if is_infinite(point_q):
        point_r = point_p
        result = point_r
        assert (isinstance(result, tuple) and
                is_on_curve(result))
        return result
    jacobian_p = to_jacobian(point_p)
    jacobian_q = to_jacobian(point_q)
    assert (isinstance(jacobian_p and jacobian_q, tuple) and
            is_on_curve_jacobian(jacobian_p) and
            is_on_curve_jacobian(jacobian_q))
    jacobian_r = jacobian_point_addition_affined_only(jacobian_p, jacobian_q)
    assert (isinstance(jacobian_r, tuple) and
            is_on_curve_jacobian(jacobian_r))
    point_r = from_jacobian(jacobian_r)
    result = point_r
    assert (isinstance(result, tuple) and
            is_on_curve(result))
    return result


def fast_scalar_multiplication(scalar: int, point: Point) -> Point:
    """
    Fast scalar multiplication of point on the Elliptic Curve over
    affine (x, y) to Jacobian coordinates (x, y, z).

    It doubles Point-P and adds Point-P with Point-Q.

    - Point is defined as Jacobian and Current is defined as
    New-Point.
    """
    # *     [ SCALAR VERIFICATION ]

    # *  -> Initial verification of arguments received by function
    # *     parameters.
    assert (isinstance(scalar, int) and
            isinstance(point, tuple) and
            is_on_curve(point))
    if (scalar in (0, N_CURVE)) or is_infinite(point):
        new_point = POINT_INFINITY_CURVE
        result = new_point
        assert (isinstance(result, tuple) and
                is_on_curve(result))
        return result
    if scalar == 1:
        result = point
        assert (isinstance(result, tuple) and
                is_on_curve(result))
        return result
    if (scalar < 0) or (scalar > N_CURVE):
        new_point = fast_scalar_multiplication(scalar % N_CURVE, point)
        result = new_point
        assert (isinstance(result, tuple) and
                is_on_curve(result))
        return result

    # *     [ SCALAR PREPARATION ]

    # *  -> The next variable below, will be assigned from the convert
    # *     the scalar (represented as an integer) to string objects as
    # *     bits ("0" or "1"), and store them inside a
    # *     Tuple["1", "0", "0", "1", ...].
    # *  -> We guarantee, as a warning in the code, that it will be
    # *     declared as Final and cannot be reassigned.
    scalar_bits: Final[ScalarBinary] = tuple(bin(scalar)[2:])
    assert is_bin_encoded("".join(scalar_bits).zfill(256))

    # *     [ SCALAR MULTIPLICATION ]

    # *  -> For-Loop below, working from left-to-right.
    # *  -> "AND IGNORING" the Most Significant Bit (MSB) from incoming
    # *     scalar. This is usually the bit farthest to the left, or the
    # *     bit transmitted first in a sequence.
    jacobian = to_jacobian(point)
    current = jacobian
    assert (isinstance(jacobian and current, tuple) and
            is_on_curve_jacobian(jacobian) and
            is_on_curve_jacobian(current))
    for bit in scalar_bits[1:]:
        current = jacobian_point_doubling(current)
        if bit == "1":
            current = jacobian_point_addition(current, jacobian)
    new_point = from_jacobian(current)
    result = new_point
    assert (isinstance(result, tuple) and
            is_on_curve(result))
    return result


# *****************************************************************************
# ********************************** [ END ] **********************************
# *****************************************************************************


# *****************************************************************************
# *                                                                           *
# *                Elliptic Curve Scalar Multiplication Test.                 *
# *                                                                           *
# *****************************************************************************
# *                                                                           *
# *****************************************************************************
# ********************************* [ BEGIN ] *********************************
# *****************************************************************************


if __name__ == "__main__":

    try:
        # *  -> Non-Singularity test ( 4⋅a³ + 27⋅b² != 0 ) for the
        # *     Elliptic Curve, and also tests if the Generator point
        # *     lies on the Elliptic Curve.
        # *  -> The program will only start if it passes the tests.
        assert ((((- 16 * (4 * pow(A_CURVE, 3, FP_CURVE) + 27 *
                           pow(B_CURVE, 2, FP_CURVE))) % FP_CURVE) !=
                 (0 % FP_CURVE)) and is_on_curve(GENERATOR_POINT_CURVE))

        from os import name as system_type
        from os import system as run_command
        from random import randrange
        from time import perf_counter

        from colorama import just_fix_windows_console

        # Get ANSI escapes from color scheme to work on Windows
        # operating system.
        just_fix_windows_console()

        # Define clear_screen function.
        def clear_screen() -> None:
            """
            Tests the operating system type and sets the screen clear
            command.
            """
            # Screen clear command for Windows operating system.
            if system_type == "nt":
                run_command("cls")
            # Screen clear command for macOS/Linux operating system.
            elif system_type == "posix":
                run_command("clear")

        start = perf_counter()
        ELAPSED = 0.0
        COUNTER = 1
        while ELAPSED < 15.0:
            private_key = randrange(1, N_CURVE)
            public_key = fast_scalar_multiplication(
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
            SECP256K1 at Jacobian Coordinates
            Copyright (C) 2021 - 2023  Gustavo Madureira
            License GNU GPL-3.0-or-later <https://gnu.org/licenses/gpl.html>
            This program comes with ABSOLUTELY NO WARRANTY.
            This is free software, and you are welcome to redistribute it
            under certain conditions.

            \33[1;7m[ *** FOR EDUCATIONAL PURPOSES ONLY *** ]\033[0m


               \033[92mPoint Number: {data[0]}
                Private Key: {data[1]}
    Uncompressed Public Key: 04{data[2]}
                               {data[3]}
      Compressed Public Key: {data[4]}{data[2]}

               Elapsed Time: {ELAPSED:.02f} seconds.
   Number of Points Created: {COUNTER} points.
\033[0m""")
            COUNTER += 1
    except Exception:  # pylint: disable=W0703
        print("""\033[91m\33[1;7m
[-----------------------------------------------------------------------------]
[                                                                             ]
[     *** Something bad happened and the program had to be shut down ***      ]
[                                                                             ]
[  Possible Errors:                                                           ]
[                                                                             ]
[ [X] Mathematical domain parameters of the elliptic curve may be incorrect ! ]
[                                                                             ]
[-----------------------------------------------------------------------------]
\033[0m""")


# *****************************************************************************
# ********************************** [ END ] **********************************
# *****************************************************************************
