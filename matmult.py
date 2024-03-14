import numpy as np

# class Matrix:
#     def __init__(self):

A = [[1, 2], [3, 4]]
B = [[5, 6], [7, 8]]
u = [3, 4, 5]
v = [3, 4, 5]
s = [[3], [4], [5]]

def dotproduct(u, v):
    dotproduct = 0
    for i in range(len(u)):
        dotproduct += u[i] * v[i]
    return dotproduct

# print(dotproduct(u, v))
        
def dotproduct_row_column (u, v):
    dotproduct = 0
    for i in range(len(u)):
        dotproduct += u[i] * v[i][0]
    return dotproduct

# print(dotproduct_row_column(v, s))

def matmult(mat1, mat2):
    new_mat = []
    m = len(mat1)
    n = len(mat1[0])
    for i in range(n):
        new_row = []
        for j in range(n):
            sum = 0
            for k in range(n):
                sum += mat1[i][k] * mat2[k][j]       
            new_row.append(sum)    
        new_mat.append(new_row)
    return new_mat

print(matmult(A, B))