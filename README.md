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