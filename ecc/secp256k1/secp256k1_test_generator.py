# pylint: disable=C0114


# Works on Python 3.8 or higher.


from typing import Final, Tuple
from colorama import just_fix_windows_console  # type: ignore
from secp256k1_jacobian import fast_scalar_multiplication
from secp256k1_weierstrass import GENERATOR_POINT_CURVE, has_even_y, \
    ec_point_multiplication


Point = Tuple[int, int]


just_fix_windows_console()


G: Final[Point] = GENERATOR_POINT_CURVE


PRIVATE_KEY = \
    0xE05AF5BC208C749190567B921A0C28FE112CD8B54E9FF82F77FA58998B694D4C


public_key_w = ec_point_multiplication(PRIVATE_KEY, G)
if has_even_y(public_key_w):
    PREFIX_W = "02"
else:
    PREFIX_W = "03"

public_key_j = fast_scalar_multiplication(PRIVATE_KEY, G)
if has_even_y(public_key_j):
    PREFIX_J = "02"
else:
    PREFIX_J = "03"

data = (hex(PRIVATE_KEY)[2:].zfill(64).upper(),
        hex(public_key_w[0])[2:].zfill(64).upper(),
        hex(public_key_w[1])[2:].zfill(64).upper(),
        PREFIX_W,
        hex(public_key_j[0])[2:].zfill(64).upper(),
        hex(public_key_j[1])[2:].zfill(64).upper(),
        PREFIX_J)

print(f"""\033[92m
     [ Weierstrass Form ]

             Private Key: {data[0]}
 Uncompressed Public Key: 04{data[1]}{data[2]}
   Compressed Public Key: {data[3]}{data[1]}


        [ Jacobian Form ]

             Private Key: {data[0]}
 Uncompressed Public Key: 04{data[4]}{data[5]}
   Compressed Public Key: {data[6]}{data[4]}
\033[0m""")
