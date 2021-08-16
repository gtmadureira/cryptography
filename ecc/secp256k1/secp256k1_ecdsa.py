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
from secp256k1_weierstrass import _GENERATOR_POINT_CURVE_, _N_CURVE_, \
    modular_inverse, ec_point_addition, ec_point_multiplication


# Type Hints.
Point = Tuple[int, int]
Signature = Tuple[int, int]


G = _GENERATOR_POINT_CURVE_
N = _N_CURVE_


def ecdsa_signature(private_key: int, message_hash: bytes) -> Signature:
    """This creates the message ECDSA-Signature."""
    data_hash = int.from_bytes(message_hash, byteorder="big")
    random_number = randrange(1, N)
    dA = private_key
    z = data_hash
    k = random_number
    xp, _ = ec_point_multiplication(k, G)
    r = xp % N
    s = (modular_inverse(k, N) * (z + r * dA)) % N
    result = (r, s)
    return result


def ecdsa_verification(public_key: Point,
                       message_hash: bytes, signature: Signature) -> bool:
    """This verifies the message ECDSA-Signature."""
    result = False
    data_hash = int.from_bytes(message_hash, byteorder="big")
    QA = public_key
    z = data_hash
    r, s = signature
    if not 0 < r < N or not 0 < s < N:
        result = False
        return result
    w = modular_inverse(s, N)
    up = (z * w) % N
    uq = (r * w) % N
    xp = ec_point_multiplication(up, G)
    xq = ec_point_multiplication(uq, QA)
    xr, _ = ec_point_addition(xp, xq)  # type: ignore
    if r % N == xr % N:
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
    public_key = \
        (0x9E0833545033AFE57A6C605D1D91B808CF495DD525217754AF9ED826E13A8AD6,
         0x3F2A1F3371300A6A0CFFE391361D9F838C67EB55B5DC6B701B0CA03668A7B8CF)

    message = b"My name is Gustavo Madureira. This is a ECDSA-Signature test."

    message_hash = hasher_double_sha256(message)
    signature = ecdsa_signature(private_key, message_hash)
    signature_verification = ecdsa_verification(
        public_key, message_hash, signature)

    if signature_verification:
        signature_result = "\033[92m[✔] Good signature\033[0m"
    else:
        signature_result = "\033[95m[X] Bad signature\033[0m"

    data = (hex(private_key)[2:].zfill(64).upper(),
            hex(x(public_key))[2:].zfill(64).upper(),
            hex(y(public_key))[2:].zfill(64).upper(),
            message.decode(),
            message_hash.hex().upper(),
            hex(signature[0])[2:].zfill(64).upper(),
            hex(signature[1])[2:].zfill(64).upper(),
            signature_result)

    print(f"""\033[92m
        SECP256K1 ECDSA-Signature  Copyright (C) 2021  Gustavo Madureira
        This program comes with ABSOLUTELY NO WARRANTY.
        This is free software, and you are welcome to redistribute it
        under certain conditions.\033[0m


           Private Key: {data[0]}

            Public Key: X = {data[1]}
                        Y = {data[2]}

               Message: {data[3]}

          Message Hash: {data[4]}

             Signature: R = {data[5]}
                        S = {data[6]}

Signature Verification: {data[7]}
""")
