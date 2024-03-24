import numpy as np
import sys

C = [[1, 1, 1, 1, 2, 2, 2, 2],
     [1, 1, 1, 1, 2, 2, 2, 2],
     [1, 1, 1, 1, 2, 2, 2, 2],
     [1, 1, 1, 1, 2, 2, 2, 2],
     [3, 3, 3, 3, 4, 4, 4, 4],
     [3, 3, 3, 3, 4, 4, 4, 4],
     [3, 3, 3, 3, 4, 4, 4, 4],
     [3, 3, 3, 3, 4, 4, 4, 4],]
D = [[1, 1, 1, 1, 1, 2, 2, 2, 2, 2],
     [1, 1, 1, 1, 1, 2, 2, 2, 2, 2],
     [1, 1, 1, 1, 1, 2, 2, 2, 2, 2],
     [1, 1, 1, 1, 1, 2, 2, 2, 2, 2],
     [1, 1, 1, 1, 1, 2, 2, 2, 2, 2],
     [3, 3, 3, 3, 3, 4, 4, 4, 4, 4],
     [3, 3, 3, 3, 3, 4, 4, 4, 4, 4],
     [3, 3, 3, 3, 3, 4, 4, 4, 4, 4],
     [3, 3, 3, 3, 3, 4, 4, 4, 4, 4],
     [3, 3, 3, 3, 3, 4, 4, 4, 4, 4]]

E = [[1,2,3],
     [1,2,3],
     [1,2,3]]
D = [[1,2],
     [1,2]]

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
    original_size = len(mat) 
    next_power_of_2 = 2**np.ceil(np.log2(original_size)).astype(int)
    pad_len = next_power_of_2 - original_size
    pad_np_mat = np.pad(mat, ((0, pad_len), (0, pad_len)), 'constant', constant_values=0)
    pad_mat = pad_np_mat.tolist()
    return pad_mat

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

    # Make dimensions a power of 2
    if len(mat1) > 0 and np.log2(len(mat1)) % 1 != 0:
        mat1 = pad_matrix(mat1)
        mat2 = pad_matrix(mat2)
    # Establish an arbitrary base case. To be adjusted
    if (len(mat1) == 2):
        mat_prod = conventional_matmult(mat1, mat2)
        return mat_prod
        
    
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
    p1 = strassen_matmult(a, matrix_subtraction(f, h))
    p2 = strassen_matmult(matrix_addition(a, b), h)
    p3 = strassen_matmult(matrix_addition(c, d), e)
    p4 = strassen_matmult(d, matrix_subtraction(g, e))
    p5 = strassen_matmult(matrix_addition(a, d), matrix_addition(e, h))
    p6 = strassen_matmult(matrix_subtraction(b, d), matrix_addition(g, h))
    p7 = strassen_matmult(matrix_subtraction(a, c), matrix_addition(e, f))

    # --- Calculate the new blocks with these 7 variables --- #
    
    upleft_block = matrix_addition(matrix_subtraction(matrix_addition(p5, p4), p2), p6)
    upright_block = matrix_addition(p1, p2)
    downleft_block = matrix_addition(p3, p4)
    downright_block = matrix_subtraction(matrix_subtraction(matrix_addition(p1, p5), p3), p7)

    
    # Merge the new blocks
    matmult = mergeblocks(upleft_block, upright_block, 
                          downleft_block, downright_block)
    
    return matmult

def main():
    # Parse input file
    dim = int(sys.argv[2])
    
    inputfile = open(sys.argv[3], "r")
    entries = []
    for line in (inputfile):
        entries.append(int(line))
    entries_len = len(entries)
    mat1 = []
    mat2 = []
    
    for i in range(0, (entries_len // 2), dim):
        mat1.append(entries[i : i + dim])
    for i in range((entries_len // 2), entries_len, dim):
        mat2.append(entries[i : i + dim])
    
    matmult = strassen_matmult(mat1, mat2)
    if len(matmult) != dim:
        matmult = [row[:dim] for row in matmult[:dim]]
    res = []
    for i in range(dim):
        res.append(matmult[i][i])
    print(res)


if __name__ == "__main__":
    main()


