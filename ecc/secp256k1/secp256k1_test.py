from secp256k1 import ec_point_multiplication

private_key = "C8AF8B035B655F8B854FFD1EA8D9E439" \
              "4BF9EEBE07ED1B8D9F03AFFFD9DBD7F5"

print()
print("              Private Key: " + private_key)
public_key = ec_point_multiplication(int("0x" + private_key, 16))
if public_key[1] % 2 == 1:
    i = "03"
else:
    i = "02"
print("  Uncompressed Public Key: " + "04" +
      hex(public_key[0])[2:].zfill(64).upper() +
      hex(public_key[1])[2:].zfill(64).upper())
print("    Compressed Public Key: " + i +
      hex(public_key[0])[2:].zfill(64).upper())
print()
