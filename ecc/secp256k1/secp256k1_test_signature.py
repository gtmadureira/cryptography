from hashlib import sha256
from secrets import randbits
from secp256k1_weierstrass import modular_inverse, _N_CURVE_
from secp256k1_weierstrass import ec_point_addition, ec_point_multiplication


private_key = int(
    "0x" +
    "634737958D20C72F558B35AC9ED9AC9F530990EF4E10E9BB5456A650EB439E9C",
    16)


def random_secure_number():
    random_number = randbits(256)
    if not 0 < random_number < _N_CURVE_:
        random_secure_number()
    else:
        return random_number


random_number = random_secure_number()
public_key = ec_point_multiplication(private_key)
message = "My name is Gustavo Madureira and this is a test message."
message_hash_hex = sha256(message.encode("utf-8")).hexdigest()
message_hash_int = int("0x" + message_hash_hex, 16)

# This creates the message signature.
xs, ys = ec_point_multiplication(random_number)
r = xs % _N_CURVE_
s = ((message_hash_int + r * private_key) *
     (modular_inverse(random_number, _N_CURVE_))) % _N_CURVE_

# This verifies the message signature.
w = modular_inverse(s, _N_CURVE_)
u1 = ec_point_multiplication((message_hash_int * w) % _N_CURVE_)
u2 = ec_point_multiplication((r * w) % _N_CURVE_, public_key)
x, y = ec_point_addition(u1, u2)  # type: ignore
if r == x:
    result = "\t\t\033[92m[✔] Good signature\033[0m"
else:
    result = "\t\t\033[95m[X] Bad signature\033[0m"

data = (hex(private_key)[2:].zfill(64).upper(),
        hex(public_key[0])[2:].zfill(64).upper(),
        hex(public_key[1])[2:].zfill(64).upper(),
        message,
        message_hash_hex.upper(),
        hex(r)[2:].zfill(64).upper(),
        hex(s)[2:].zfill(64).upper(),
        result)

print(f"""
 Private Key: {data[0]}

  Public Key: X = {data[1]}
              Y = {data[2]}

     Message: {data[3]}

Message Hash: {data[4]}

   Signature: R = {data[5]}
              S = {data[6]}

              {result}
""")
