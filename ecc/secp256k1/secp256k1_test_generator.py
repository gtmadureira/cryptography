import secp256k1_weierstrass
import secp256k1_jacobian


private_key = "634737958D20C72F558B35AC9ED9AC9F" \
              "530990EF4E10E9BB5456A650EB439E9C"


public_key_w = secp256k1_weierstrass.ec_point_multiplication(
    int("0x" + private_key, 16))
if public_key_w[1] % 2 == 1:
    i = "03"
else:
    i = "02"

print("\n\033[92m\n" +
      "[ Weierstrass Form ]" + "\n"
      "            Private Key: " + private_key + "\n"
      "Uncompressed Public Key: " + "04" +
      hex(public_key_w[0])[2:].zfill(64).upper() +
      hex(public_key_w[1])[2:].zfill(64).upper() + "\n"
      "  Compressed Public Key: " + i +
      hex(public_key_w[0])[2:].zfill(64).upper() + "\n" +
      "\033[0m")


public_key_j = secp256k1_jacobian.ec_point_multiplication(
    int("0x" + private_key, 16))
if public_key_j[1] % 2 == 1:
    i = "03"
else:
    i = "02"

print("\n\033[92m\n" +
      "   [ Jacobian Form ]" + "\n"
      "            Private Key: " + private_key + "\n"
      "Uncompressed Public Key: " + "04" +
      hex(public_key_j[0])[2:].zfill(64).upper() +
      hex(public_key_j[1])[2:].zfill(64).upper() + "\n"
      "  Compressed Public Key: " + i +
      hex(public_key_j[0])[2:].zfill(64).upper() + "\n" +
      "\033[0m")
