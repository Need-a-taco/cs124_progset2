import random

f = open("inputfile.txt", "a")
n = 257
big_n = 2 * (257**2)
for _ in range(big_n):
    x = random.randint(0, 10)
    f.write(f"{x}\n")