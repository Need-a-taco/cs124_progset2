import numpy as np
import sys
import time
import random

###########################################################
#                                                         #
# These are all helper functions for strassen's algorithm #
#                                                         #
###########################################################

def get_matrix_block(mat, x0, x1, y0, y1):
        matrix_block = []
        for i in range(x0, x1):
            new_row = []
            for j in range(y0, y1):    
                new_row.append(mat[i][j])
            matrix_block.append(new_row)
        return matrix_block

def mergeblocks(upleft, upright, downleft, downright):
    new_mat = []
    for i in range(len(upleft)):
        new_mat.append(upleft[i] + upright[i])
    
    for i in range(len(downleft)):
        new_mat.append(downleft[i] + downright[i])
        
    return new_mat

def pad_matrix(mat):
    for row in mat:
        row.append(0)
    new_row_length = len(mat[0])  # This now includes the newly added column of zeros.
    mat.append([0] * new_row_length)
    return mat

def matrix_addition(mat1, mat2):
    new_mat = []
    n = len(mat1)
    for i in range(n):
        new_row = []
        for j in range(n):
            new_row.append(mat1[i][j] + mat2[i][j])
        new_mat.append(new_row)
    return new_mat

def matrix_subtraction(mat1, mat2):
    new_mat = []
    n = len(mat1)
    for i in range(n):
        new_row = []
        for j in range(n):
            new_row.append(mat1[i][j] - mat2[i][j])
        new_mat.append(new_row)
    return new_mat

#############################################################
#                                                           #
# Now onto the actual matrix multiplication implementations #
#                                                           #
#############################################################

def conventional_matmult(mat1, mat2):
    new_mat = []
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

def strassen_matmult(mat1, mat2, threshold):
    init_len = len(mat1)
    if (init_len <= threshold):
        mat_prod = conventional_matmult(mat1, mat2)
        return mat_prod
    if init_len % 2 == 1:
        mat1 = pad_matrix(mat1)
        mat2 = pad_matrix(mat2)
    # Establish an arbitrary base case. To be adjusted

        
    
    # --- Divide the matrices into blocks --- #
    n = len(mat1)
    x = y = (n / 2)
    x = int(x)
    y = int(y)
    
    a = get_matrix_block(mat1, 0, x, 0, y)
    e = get_matrix_block(mat2, 0, x, 0, y) # Upper right blocks
    
    b = get_matrix_block(mat1, 0, x, y, n)
    f = get_matrix_block(mat2, 0, x, y, n) # Upper left blocks
    
    c = get_matrix_block(mat1, x, n, 0, y)
    g = get_matrix_block(mat2, x, n, 0, y) # Lower left blocks
    
    d = get_matrix_block(mat1, x, n, y, n)
    h = get_matrix_block(mat2, x, n, y, n) # Lower right blocks

    
    # --- Calculate multiplication with only 7 variables --- # 
    p1 = strassen_matmult(a, matrix_subtraction(f, h), threshold)
    p2 = strassen_matmult(matrix_addition(a, b), h, threshold)
    p3 = strassen_matmult(matrix_addition(c, d), e, threshold)
    p4 = strassen_matmult(d, matrix_subtraction(g, e), threshold)
    p5 = strassen_matmult(matrix_addition(a, d), matrix_addition(e, h), threshold)
    p6 = strassen_matmult(matrix_subtraction(b, d), matrix_addition(g, h), threshold)
    p7 = strassen_matmult(matrix_subtraction(a, c), matrix_addition(e, f), threshold)

    # --- Calculate the new blocks with these 7 variables --- #
    
    upleft_block = matrix_addition(matrix_subtraction(matrix_addition(p5, p4), p2), p6)
    upright_block = matrix_addition(p1, p2)
    downleft_block = matrix_addition(p3, p4)
    downright_block = matrix_subtraction(matrix_subtraction(matrix_addition(p1, p5), p3), p7)
    
    # Merge the new blocks
    matmult = mergeblocks(upleft_block, upright_block, 
                          downleft_block, downright_block)
    
    return matmult

#######################################
#                                     #
# Implement the Programming Set Tasks #
#                                     #
#######################################

# Task 1
def testing_threshold():
    dimension = 256
    threshold = 8
    stras_list = []
    conv_list = []
    for dim in range(2, dimension):
        matrix =[[random.randint(0, 20) for _ in range(dim)] for _ in range(dim)]
        #strassens
        stras_avg_runtime = 0
        for _ in range(5):
            start = time.time()
            matmult = strassen_matmult(matrix, matrix, dim-1)
            matmult = [row[:dim] for row in matmult[:dim]]
            end = time.time()
            stras_avg_runtime += (end - start)
        stras_avg_runtime /= 5
        stras_list.append(stras_avg_runtime)
        
        #conventional
        conv_avg_runtime = 0
        for _ in range(5):
            start = time.time()
            matmult = conventional_matmult(matrix, matrix)
            end = time.time()
            conv_avg_runtime += (end - start)
        conv_avg_runtime /= 5
        conv_list.append(conv_avg_runtime)
        
    
# Task 2: What we actually submit to Gradescope
def matrix_multplication():
    # Parse input file
    dimension = int(sys.argv[2])
    
    inputfile = open(sys.argv[3], "r")
    entries = []
    for line in (inputfile):
        entries.append(int(line))
    entries_len = len(entries)
    mat1 = []
    mat2 = []
    
    for i in range(0, (entries_len // 2), dimension):
        mat1.append(entries[i : i + dimension])
    for i in range((entries_len // 2), entries_len, dimension):
        mat2.append(entries[i : i + dimension])
    matmult = strassen_matmult(mat1, mat2, ((dimension // 2) +1))
    
    if dimension % 2 == 1:
        matmult = [row[:dimension] for row in matmult[:dimension]]
    for i in range(dimension):
        print(matmult[i][i])


# Task 3: Finding triangles in graphs, also purely testing
def count_triangles(n, p): 
    adj_mat = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(len(adj_mat)):
        for j in range(i,len(adj_mat)):
            val = np.random.binomial(1,p)
            adj_mat[i][j] = val
            if i != j:
                adj_mat[j][i] = val
    # Find the cube of adjacency matrix
    adj_mat_cubed = strassen_matmult(adj_mat, strassen_matmult(adj_mat, adj_mat, 50), 50)
    print("\n")
    print(adj_mat_cubed)
    
    num_triangles = 0
    for i in range(n):
        num_triangles += adj_mat_cubed[i][i]
    num_triangles //= 6
    print(num_triangles)
    
        
def main():
    matrix_multplication()

if __name__ == "__main__":
    main()
