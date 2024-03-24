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
    def mat_to_lst(mat):
        new_lst = []
        for row in mat:
            for element in row:
                new_lst.append(element)
        return new_lst

    def add_from_block(lst, mat):
        if mat:
            lst.extend(mat[0])

    ul = mat_to_lst(upleft)
    ur = mat_to_lst(upright)
    dl = mat_to_lst(downleft)
    dr = mat_to_lst(downright)

    new_mat = []
    n = len(upleft)
    cutoff = n // 2
    for i in range(n):
        new_row = []
        for j in range(n):
            if (i < cutoff) and (j < cutoff):
                add_from_block(new_row, ul)
            elif (i < cutoff) and (j >= cutoff):
                add_from_block(new_row, ur)
            elif (i >= cutoff) and (j < cutoff):
                add_from_block(new_row, dl)
            else:
                add_from_block(new_row, dr)
        new_mat.append(new_row)

    return new_mat

def pad_matrix(mat):
    np_mat = np.array(mat)
    pad_np_mat = np.pad(np_mat, ((0, 1), (0, 1)), mode='constant', constant_values=0)
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

    # Make dimensions are even, not odd.
    if (len(mat1) % 2 != 0):
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
    
    upleft_block = matrix_subtraction(matrix_addition(p5, p4), 
                                matrix_addition(p2, p6))
    upright_block = matrix_addition(p1, p2)
    downleft_block = matrix_addition(p3, p4)
    downright_block = matrix_subtraction(matrix_addition(p1, p5), 
                                         matrix_subtraction(p3, p7))
    
    # Merge the new blocks
    matmult = mergeblocks(upleft_block, upright_block, 
                          downleft_block, downright_block)
    
    return matmult

def main():
    # Parse input file
    dim = int(sys.argv[2])
    
    inputfile = open("inputfile.txt", "r")
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
    print(matmult)

if __name__ == "__main__":
    main()


