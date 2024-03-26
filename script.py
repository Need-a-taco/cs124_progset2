import numpy as np

def include_edge(p):     
    result = np.random.binomial(1,p)
    return result

n = 1024
matrix_dim = n**2
p = 0.01

matrices = [124] * (matrix_dim * 2)

with open('inputfile.txt', 'w+') as matrix_file:


    for _ in range(matrix_dim):
        x = include_edge(p)
        matrix_file.write(f"{x}\n")