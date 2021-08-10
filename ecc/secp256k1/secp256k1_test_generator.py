import secp256k1_jacobian
import secp256k1_weierstrass


private_key = "634737958D20C72F558B35AC9ED9AC9F" \
              "530990EF4E10E9BB5456A650EB439E9C"


public_key_w = secp256k1_weierstrass.ec_point_multiplication(
    int("0x" + private_key, 16))
if public_key_w[1] % 2 == 1:
    prefix_w = "03"
else:
    prefix_w = "02"

public_key_j = secp256k1_jacobian.ec_point_multiplication(
    int("0x" + private_key, 16))
if public_key_j[1] % 2 == 1:
    prefix_j = "03"
else:
    prefix_j = "02"

data = (private_key,
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
