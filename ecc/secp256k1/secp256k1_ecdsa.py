"""
[ LICENSE ] ( GNU GPLv3 ) :

SECP256K1 ECDSA-Signature - Elliptic Curve Cryptography (ECC).
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

Module to creates and verifies ECDSA-Signature in secp256k1 elliptic
curve, using Weierstrass form.


    Source:

            https://github.com/gtmadureira/cryptography/blob/main/ecc/secp256k1/secp256k1_ecdsa.py

    Author:

            • Gustavo Madureira (gtmadureira@gmail.com)
            • https://gtmadureira.github.io/
"""


from typing import Tuple
from random import randrange
from secp256k1_weierstrass import modular_inverse, _N_CURVE_
from secp256k1_weierstrass import ec_point_multiplication, ec_point_addition


# Type Hints.
Point = Tuple[int, int]


def ecdsa_signature(private_key: int, data_hash: bytes) -> tuple:
    """This creates the message ECDSA-Signature."""
    random_number = randrange(1, _N_CURVE_)
    xp, _ = ec_point_multiplication(random_number)
    r = xp % _N_CURVE_
    message_hash = int.from_bytes(data_hash, byteorder="big")
    s = ((message_hash + r * private_key) *
         (modular_inverse(random_number, _N_CURVE_))) % _N_CURVE_
    result = (r, s)
    return result


def ecdsa_verification(public_key: Point,
                       data_hash: bytes, signature: tuple) -> bool:
    """This verifies the message ECDSA-Signature."""
    result = False
    message_hash = int.from_bytes(data_hash, byteorder="big")
    r, s = signature
    w = modular_inverse(s, _N_CURVE_)
    up = ec_point_multiplication((message_hash * w) % _N_CURVE_)
    uq = ec_point_multiplication((r * w) % _N_CURVE_, public_key)
    x, _ = ec_point_addition(up, uq)  # type: ignore
    if r == x:
        result = True
    return result


if __name__ == "__main__":

    from hashlib import sha256
    from secp256k1_weierstrass import x, y

    def hasher_double_sha256(data: bytes) -> bytes:
        """Get double hash through SHA-256."""
        result = sha256(sha256(data).digest()).digest()
        return result

    # ECDSA-Signature test.
    private_key = \
        0xE05AF5BC208C749190567B921A0C28FE112CD8B54E9FF82F77FA58998B694D4C
    public_key = ec_point_multiplication(private_key)

    message = b"My name is Gustavo Madureira. This is a ECDSA-Signature test."

    message_hash = hasher_double_sha256(message)
    signature = ecdsa_signature(private_key, message_hash)
    signature_verification = ecdsa_verification(
        public_key, message_hash, signature)

    if signature_verification:
        signature_result = "\033[92m[✔] Good signature\033[0m"
    else:
        signature_result = "\033[95m[X] Bad signature\033[0m"

    print(f"""\033[92m
        SECP256K1 ECDSA-Signature  Copyright (C) 2021  Gustavo Madureira
        This program comes with ABSOLUTELY NO WARRANTY.
        This is free software, and you are welcome to redistribute it
        under certain conditions.\033[0m


           Private Key: {hex(private_key)[2:].zfill(64).upper()}

            Public Key: X = {hex(x(public_key))[2:].zfill(64).upper()}
                        Y = {hex(y(public_key))[2:].zfill(64).upper()}

               Message: {message.decode()}

          Message Hash: {message_hash.hex().upper()}

             Signature: R = {hex(signature[0])[2:].zfill(64).upper()}
                        S = {hex(signature[1])[2:].zfill(64).upper()}

Signature Verification: {signature_result}
""")
