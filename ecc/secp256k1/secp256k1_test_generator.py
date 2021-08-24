from secp256k1_weierstrass import _GENERATOR_POINT_CURVE_, int_from_hex, \
    has_even_y, ec_point_multiplication
from secp256k1_jacobian import fast_jacobian_point_multiplication


G = _GENERATOR_POINT_CURVE_


private_key = int_from_hex(
    "E05AF5BC 208C7491 90567B92 1A0C28FE 112CD8B5 4E9FF82F 77FA5899 8B694D4C")


public_key_w = ec_point_multiplication(private_key, G)
if has_even_y(public_key_w):
    prefix_w = "02"
else:
    prefix_w = "03"

public_key_j = fast_jacobian_point_multiplication(private_key, G)
if has_even_y(public_key_j):
    prefix_j = "02"
else:
    prefix_j = "03"

data = (hex(private_key)[2:].zfill(64).upper(),
        hex(public_key_w[0])[2:].zfill(64).upper(),
        hex(public_key_w[1])[2:].zfill(64).upper(),
        prefix_w,
        hex(public_key_j[0])[2:].zfill(64).upper(),
        hex(public_key_j[1])[2:].zfill(64).upper(),
        prefix_j)

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
