# cs124_progset2
Pedro and Josh, the dream team again

# Matrix Multiplication
To compile the code... use the command

<g++ -o matmult matmult.cpp>

To talk through the code real quick, we use three for
loops to implement conventional grade school matrix
multiplication. This makes sense, because this is an
O(n^3) algorithm. 

Josh notes: to implement matrix multiplication, it helped
to think through a dot product implementation. 

def dotproduct(u, v):
    dotproduct = 0
    for i in range(len(u)):
        dotproduct += u[i] * v[i]
    return dotproduct
        
def dotproduct_row_column (u, v):
    dotproduct = 0
    for i in range(len(u)):
        dotproduct += u[i] * v[i][0]
    return dotproduct

def include_edge(p):     
    result = np.random.binomial(1,p)
    return result


adj_mat = []
for i in range(n):
    new_row = [124] * n
    adj_mat.append(new_row)
for i in range(n):
    for j in range(n):
        if (adj_mat[i][j] == 124):
            adj_mat[i][j] = adj_mat[j][i] = include_edge(p)