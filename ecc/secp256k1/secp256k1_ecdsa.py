"""
[ LICENSE ]  ( GNU GPLv3 )  :

SECP256K1 ECDSA-Signature  -  Elliptic Curve Cryptography (ECC).
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
Module to creates and verifies ECDSA-Signature with elliptic curve
SECP256K1, using the Jacobian form.

Works on Python 3.8 or higher.


    Source:

            https://github.com/gtmadureira/cryptography/blob/main/ecc/secp256k1/secp256k1_ecdsa.py


    Author:

            • Gustavo Madureira (gtmadureira@gmail.com)
            • https://gtmadureira.github.io/


    Author Nickname:

            • Hattoshi Hanzōmoto (hanzomoto_hattoshi@protonmail.com)
"""


from random import randrange
from typing import Final, Tuple
from colorama import just_fix_windows_console  # type: ignore
from secp256k1_jacobian import GENERATOR_POINT_CURVE, N_CURVE, \
    modular_inverse, fast_point_addition, fast_scalar_multiplication


# Type Hints.
Point = Tuple[int, int]
Signature = Tuple[int, int]


# Get ANSI escapes from color scheme to work on Windows.
just_fix_windows_console()


# Generator point and curve order.
G: Final[Point] = GENERATOR_POINT_CURVE
N: Final[int] = N_CURVE


def ecdsa_signature(private_key: int, message_hash: bytes) -> Signature:
    """This creates the message ECDSA-Signature."""
    data_hash = int.from_bytes(message_hash, byteorder="big")
    random_number = randrange(1, N)
    _da = private_key
    _z = data_hash
    _k = random_number
    _xp, _ = fast_scalar_multiplication(_k, G)
    _r = _xp % N
    _s = (modular_inverse(_k, N) * (_z + _r * _da)) % N
    result = (_r, _s)
    return result


def ecdsa_verification(
        public_key: Point, message_hash: bytes, signature: Signature) -> bool:
    """This verifies the message ECDSA-Signature."""
    result = False
    data_hash = int.from_bytes(message_hash, byteorder="big")
    _qa = public_key
    _z = data_hash
    _r, _s = signature
    if not 0 < _r < N or not 0 < _s < N:
        result = False
        return result
    _w = modular_inverse(_s, N)
    _sb = (_z * _w) % N
    _sd = (_r * _w) % N
    _p = fast_scalar_multiplication(_sb, G)
    _q = fast_scalar_multiplication(_sd, _qa)
    _xr, _ = fast_point_addition(_p, _q)
    if _r % N == _xr % N:
        result = True
        return result
    return result


if __name__ == "__main__":

    # ECDSA-Signature test.
    from hashlib import sha256
    from secp256k1_jacobian import x_coordinate, y_coordinate

    def hasher_double_sha256(data: bytes) -> bytes:
        """Get double hash through SHA-256."""
        result = sha256(sha256(data).digest()).digest()
        return result

    PRIVATE_KEY = \
        0xE05AF5BC208C749190567B921A0C28FE112CD8B54E9FF82F77FA58998B694D4C

    PUBLIC_KEY = fast_scalar_multiplication(PRIVATE_KEY, G)

    MESSAGE = b"My name is Gustavo Madureira. This is a ECDSA-Signature test."

    MESSAGE_HASH = hasher_double_sha256(MESSAGE)
    SIGNATURE = ecdsa_signature(PRIVATE_KEY, MESSAGE_HASH)
    SIGNATURE_VERIFICATION = ecdsa_verification(
        PUBLIC_KEY, MESSAGE_HASH, SIGNATURE)

    if SIGNATURE_VERIFICATION:
        SIGNATURE_RESULT = "\033[92m[✔] Good signature\033[0m"
    else:
        SIGNATURE_RESULT = "\033[95m[X] Bad signature\033[0m"

    DATA = (hex(PRIVATE_KEY)[2:].zfill(64).upper(),
            hex(x_coordinate(PUBLIC_KEY))[2:].zfill(64).upper(),
            hex(y_coordinate(PUBLIC_KEY))[2:].zfill(64).upper(),
            MESSAGE.decode(),
            MESSAGE_HASH.hex().upper(),
            hex(SIGNATURE[0])[2:].zfill(64).upper(),
            hex(SIGNATURE[1])[2:].zfill(64).upper(),
            SIGNATURE_RESULT)

    print(f"""\033[92m
        SECP256K1 ECDSA-Signature  Copyright (C) 2021  Gustavo Madureira
        License GNU GPL-3.0-or-later <https://gnu.org/licenses/gpl.html>
        This program comes with ABSOLUTELY NO WARRANTY.
        This is free software, and you are welcome to redistribute it
        under certain conditions.

        \33[1;7m *** FOR EDUCATIONAL PURPOSES ONLY *** \033[0m


           Private Key: {DATA[0]}

            Public Key: X = {DATA[1]}
                        Y = {DATA[2]}

               Message: {DATA[3]}

          Message Hash: {DATA[4]}

             Signature: R = {DATA[5]}
                        S = {DATA[6]}

Signature Verification: {DATA[7]}
""")
