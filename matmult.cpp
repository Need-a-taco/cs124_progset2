#include <vector>
#include <iostream>

using namespace std;

vector<vector<int>> 
matmult(const vector<vector<int>>& mat1,
        const vector<vector<int>>& mat2) {
    vector<vector<int>> new_mat;
    int m = mat1.size();  // Number of rows in mat1
    int n = mat2[0].size();  // Number of columns in mat2
    int p = mat1[0].size();  // Number of columns in mat1 or rows in mat2

    for (int i = 0; i < m; i++) 
    {
        vector<int> new_row;
        for (int j = 0; j < n; j++) 
        {
            int sum = 0;
            for (int k = 0; k < p; k++) 
            {
                sum += mat1[i][k] * mat2[k][j];
            }
            new_row.push_back(sum);
        }
        new_mat.push_back(new_row);
    }
    return new_mat;
}


int main() {
    vector<vector<int>> matrix1 = {{1, 2}, {3, 4}};
    vector<vector<int>> matrix2 = {{5, 6}, {7, 8}};

    vector<vector<int>> resultMatrix = matmult(matrix1, matrix2);

    // Output the result matrix
    for (const auto& row : resultMatrix) {
        for (int val : row) {
            cout << val << " ";
        }
        cout << endl;
    }

    return 0;
}

