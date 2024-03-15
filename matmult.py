import numpy as np

# class Matrix:
#     def __init__(self):

A = [[1, 2], [3, 4]]
B = [[5, 6], [7, 8]]
C = [[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]]
D = [[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]]
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


def strassen_matmult(mat1, mat2):
    # Meant for square matrices
    n = len(mat1)
    upleft_block1 = []
    upright_block1 = []
    downleft_block1 = []
    downright_block1 = []
    upleft_block2 = []
    upright_block2 = []
    downleft_block2 = []
    downright_block2 = []
    # This is wrong, needs to make new rows...
    for i in range(n):
        for j in range(n):
            if (i < n/2) and (j < n/2):
                upleft_block1.append(mat1[i][j])
            elif (i < n/2) and (j > n/2):
                upright_block1.append(mat1[i][j])
            elif (i > n/2) and (j < n/2):
                downleft_block1.append(mat1[i][j])
            else:
                downright_block1.append(mat1[i][j])
    for i in range(n):
        for j in range(n):
            if (i < n/2) and (j < n/2):
                upleft_block2.append(mat2[i][j])
            elif (i < n/2) and (j > n/2):
                upright_block2.append(mat2[i][j])
            elif (i > n/2) and (j < n/2):
                downleft_block2.append(mat2[i][j])
            else:
                downright_block2.append(mat2[i][j])
    print(upleft_block1, upleft_block2)
    
print(strassen_matmult(C, D))