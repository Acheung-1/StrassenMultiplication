#include <stdio.h> 
#include <mpi.h>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>
using namespace std;


vector<vector<vector<int>>> quadrants(const vector<vector<int>> A) {
    // Function splits a matrix into four quadrants, starting from top left and going clockwise
    // Input: a matrix
    // Output: vector of matrixes {A11, A12, A21, A22}

    vector<int> newRowA;
    vector<vector<int>> one, two, three, four;
    int rowA = A.size();
    
    
    // For input matrix
    for (int i = 0; i < A.size(); i++) {
        newRowA = {};
        for (int j = 0; j < A[0].size(); j++) {
            if (j == rowA/2) {
                if (i < rowA/2) {
                    one.push_back(newRowA);
                }
                if (i >= rowA/2) {
                    three.push_back(newRowA);
                }
                newRowA = {};
            }
            newRowA.push_back(A[i][j]);
        }

        if (i < rowA/2) {
            two.push_back(newRowA);
        }
        if (i >= rowA/2) {
            four.push_back(newRowA);
        }
    }

    return {one, two, three, four};
}

vector<vector<int>> combine(const vector<vector<vector<int>>> seven) {
    // Function combines the 7 matrix multiplications from strassen and outputs the complete matrix
    // Input: a vector of 7 matrixes {M1, M2, M3, M4 ,M5 ,M6 ,M7}
    // Output: combined matrix

    vector<int> newRow;
    vector<vector<int>> multi1, multi2, multi3, multi4, multi5, multi6, multi7;
    vector<vector<int>> result;


    multi1 = seven[0];
    multi2 = seven[1];
    multi3 = seven[2];
    multi4 = seven[3];
    multi5 = seven[4];
    multi6 = seven[5];
    multi7 = seven[6];

    int h = multi1.size();

    for (int i = 0; i < h; ++i) {
        newRow = {};
        for (int j = 0; j < 2*h; ++j) {
            if (j < h) {
                newRow.push_back(multi1[i][j]+multi4[i][j]-multi5[i][j]+multi7[i][j]);
            }
            else if (j >= h){
                newRow.push_back(multi3[i][j-h]+multi5[i][j-h]);
            }
        }
        result.push_back(newRow);
    }

    for (int i = 0; i < h; ++i) {
        newRow = {};
        for (int j = 0; j < 2*h; ++j) {
            if (j < h) {
                newRow.push_back(multi2[i][j]+multi4[i][j]);
            }
            else if (j >= h){
                newRow.push_back(multi1[i][j-h]-multi2[i][j-h]+multi3[i][j-h]+multi6[i][j-h]);
            }
        }
        result.push_back(newRow);
    }

    return result;

}

vector<vector<vector<int>>> sevenSplit(const vector<vector<int>> A, const vector<vector<int>> B) {
    // Function takes in two matrixes and returns 14 submatrixes that need to be multiplied
    // Input: 2 matrixes
    // Output: a vector of 14 matrixes {Add1,Add2,Add3,B11,A11,Add4,A22,Add5,Add6,B22,Add7,Add8,Add9,Add10}
    // Add1*Add2, Add3*B11, A11*Add4, etc... from Strassen
    

    vector<vector<vector<int>>> seven;
    
    vector<vector<int>> A11, A12, A21, A22, B11, B12, B21, B22;

    vector<vector<vector<int>>> resultsA, resultsB;
    resultsA = quadrants(A);
    resultsB = quadrants(B);
    A11 = resultsA[0];
    A12 = resultsA[1];
    A21 = resultsA[2];
    A22 = resultsA[3];
    B11 = resultsB[0];
    B12 = resultsB[1];
    B21 = resultsB[2];
    B22 = resultsB[3];

    vector<vector<int>> Add1(A11.size(), vector<int>(A11[0].size()));
    vector<vector<int>> Add2(A11.size(), vector<int>(A11[0].size()));
    vector<vector<int>> Add3(A11.size(), vector<int>(A11[0].size()));
    vector<vector<int>> Add4(A11.size(), vector<int>(A11[0].size()));
    vector<vector<int>> Add5(A11.size(), vector<int>(A11[0].size()));
    vector<vector<int>> Add6(A11.size(), vector<int>(A11[0].size()));
    vector<vector<int>> Add7(A11.size(), vector<int>(A11[0].size()));
    vector<vector<int>> Add8(A11.size(), vector<int>(A11[0].size()));
    vector<vector<int>> Add9(A11.size(), vector<int>(A11[0].size()));
    vector<vector<int>> Add10(A11.size(), vector<int>(A11[0].size()));

    for (int i = 0; i < A11.size(); ++i) {
        for (int j = 0; j < A11[0].size(); ++j) {
            Add1[i][j] = A11[i][j]+A22[i][j];
            Add2[i][j] = B11[i][j]+B22[i][j];
            Add3[i][j] = A21[i][j]+A22[i][j];
            Add4[i][j] = B12[i][j]-B22[i][j];
            Add5[i][j] = B21[i][j]-B11[i][j];
            Add6[i][j] = A11[i][j]+A12[i][j];
            Add7[i][j] = A21[i][j]-A11[i][j];
            Add8[i][j] = B11[i][j]+B12[i][j];
            Add9[i][j] = A12[i][j]-A22[i][j];
            Add10[i][j] = B21[i][j]+B22[i][j];
        }
    }


    seven = {Add1,Add2,Add3,B11,A11,Add4,A22,Add5,Add6,B22,Add7,Add8,Add9,Add10};


    return seven;
}

vector<vector<int>> strassen(const vector<vector<int>> A, const vector<vector<int>> B) {
    // Recursive strassen function
    // Input: two matrixes
    // Output: expected matrix after multiplication of inputs

    int rowA, colA, rowB, colB;
    vector<vector<int>> answer;
    

    rowA = A.size();
    colA = A[0].size();
    rowB = B.size();
    colB = B[0].size();

    if (rowA == 2) {
        int M1, M2, M3, M4, M5, M6, M7;
        M1 = (A[0][0]+A[1][1])*(B[0][0]+B[1][1]);
        M2 = (A[1][0]+A[1][1])*B[0][0];
        M3 = A[0][0]*(B[0][1]-B[1][1]);
        M4 = A[1][1]*(B[1][0]-B[0][0]);
        M5 = (A[0][0]+A[0][1])*B[1][1];
        M6 = (A[1][0]-A[0][0])*(B[0][0]+B[0][1]);
        M7 = (A[0][1]-A[1][1])*(B[1][0]+B[1][1]);

        answer.push_back({M1+M4-M5+M7, M3+M5});
        answer.push_back({M2+M4, M1-M2+M3+M6});


        return answer;
    }

    vector<vector<vector<int>>> split;

    split = sevenSplit(A, B);

    vector<vector<int>> M1, M2, M3, M4, M5, M6, M7;

    M1 = strassen(split[0],split[1]);
    M2 = strassen(split[2],split[3]);
    M3 = strassen(split[4],split[5]);
    M4 = strassen(split[6],split[7]);
    M5 = strassen(split[8],split[9]);
    M6 = strassen(split[10],split[11]);
    M7 = strassen(split[12],split[13]);


    answer = combine({M1, M2, M3, M4, M5, M6, M7});


    return answer; // Return the result matrix
}



int main (int argc, char *argv[]) { 
    int rank, size, rowA, colA, rowB, colB, square, matrixSize;
    
    vector<int> scatterArray;
    int chunks;
    vector<vector<int>> result;

    // Starting time
    auto start_time = std::chrono::high_resolution_clock::now();
    

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank); 
    MPI_Comm_size (MPI_COMM_WORLD, &size); 

    srand(time(0) + rank);

    if (rank == 0) {

        // Creating the two matrixes, a matrixSize x matrixSize matrix 
        // Can control size by changing matrixSize
        // All elements within matrix will be randomly generated numbers from 0-10
        vector<vector<int>> matrixA, matrixB;
        vector<int> genA, genB;

        // edit this to change size of matrix

        matrixSize = 512;
        for (int i = 0; i < matrixSize; i++) {
            genA = {};
            genB = {};

            for (int j = 0; j < matrixSize; j++) {
                int ele = rand() % 10;
                genA.push_back(ele);
                ele = rand() % 10;
                genB.push_back(ele);
            }
            matrixA.push_back(genA);
            matrixB.push_back(genB);
        }

        rowA = matrixA.size();
        colA = matrixA[0].size();
        rowB = matrixB.size();
        colB = matrixB[0].size();

        square = max({rowA, colA, rowB, colB});

        // If matrix isn't a power of two, increase the size of it so it becomes the next power of two

        int i = 0;
        while (pow(2, i) < square) {
            i++;
        }

        vector<vector<int>> newA(pow(2, i), vector<int>(pow(2, i),0));
        vector<vector<int>> newB(pow(2, i), vector<int>(pow(2, i),0));


        // Copying values from matrixA into newA and matrixB into newB
        for (int i = 0; i < matrixA.size(); ++i) {
            for (int j = 0; j < matrixA[i].size(); ++j) {
                newA[i][j] = matrixA[i][j];
            }
        }
        for (int i = 0; i < matrixB.size(); ++i) {
            for (int j = 0; j < matrixB[i].size(); ++j) {
                newB[i][j] = matrixB[i][j];
            }
        }


        matrixA = newA;
        matrixB = newB;

        // printf("MatrixA : \n");
        // for (int i = 0; i < matrixA.size(); i++) {
        //     for (int j = 0; j < matrixA[0].size(); j++) {
        //         printf("%d ", matrixA[i][j]);
        //     }
        //     printf("\n");
        // }
        // printf("MatrixB : \n");
        // for (int i = 0; i < matrixB.size(); i++) {
        //     for (int j = 0; j < matrixB[0].size(); j++) {
        //         printf("%d ", matrixB[i][j]);
        //     }
        //     printf("\n");
        // }

        // Setting up strassen BASED on the number of nodes

        // 7 nodes

        // Splitting the matrix into 14 pieces (every 2 matrixes need to be multiplied together)


        vector<vector<vector<int>>> split1;
        split1 = sevenSplit(matrixA, matrixB);
        
        int rows = matrixA.size()/2;
        chunks = 2*rows*rows;


        // Turn into vector format to scatter

        
        for (int i = 0; i < split1.size(); i++) {
            for (int j = 0; j < rows; j++) {
                for (int k = 0; k < rows; k++) {
                    scatterArray.push_back(split1[i][j][k]);
                }
            }
        }
    }
    else {
        int rowA, colA, rowB, colB;
        vector<vector<int>> res;
        vector<int> newRowA, newRowB;
        vector<vector<vector<int>>> resultsA, resultsB;
        vector<vector<int>> A11, A12, A21, A22, B11, B12, B21, B22;
        vector<vector<int>> Add1, Add2, Add3, Add4, Add5, Add6, Add7, Add8, Add9, Add10;
    }

    // Identify how many chunks each node gets
    MPI_Bcast(&chunks, 1, MPI_INT, 0, MPI_COMM_WORLD);

    vector<int> matrixes(chunks);
    MPI_Scatter(scatterArray.data(), chunks, MPI_INT, matrixes.data(), chunks, MPI_INT, 0, MPI_COMM_WORLD);
    

    // Recreating the matrixes that need to be multiplied together
    vector<vector<int>> subA, subB;

    int subMatrixSize = matrixes.size()/2;
    vector<int> rows;
    int limit = static_cast<int>(sqrt(subMatrixSize));


    for (int i = 0; i < matrixes.size(); i++) {
        rows.push_back(matrixes[i]);
        if ((i+1) % limit == 0) {
            if (i < subMatrixSize) {
                subA.push_back(rows);
            }
            else {
                subB.push_back(rows);
            }
            rows = {};
        }
    }


    
    // Use recursive strassen to multiply the two submatrixes from the scatter
    vector<vector<int>> multi(subA.size(), vector<int>(subA[0].size()));

    multi = strassen(subA,subB);

    vector<int> gatherArray(multi.size()*multi[0].size());


    for (int i = 0; i < multi.size(); i++) {
        for (int j = 0; j < multi[0].size(); j++) {
            gatherArray[i*multi[0].size()+j] = multi[i][j];
        }
    } 



    vector<int> answer(size*gatherArray.size());

    MPI_Gather(gatherArray.data(), gatherArray.size(), MPI_INT, answer.data(), gatherArray.size(), MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {

        // Putting back into matrix form
        vector<vector<int>> allRows, mul1, mul2, mul3, mul4, mul5, mul6, mul7;
        int length = answer.size()/7;
        int rowSize = static_cast<int>(sqrt(length));
        vector<vector<vector<int>>> matrixSep;
        vector<vector<int>> subMat;
        rows = {};

        for (int i = 0; i < answer.size(); i++) {
            rows.push_back(answer[i]);
            if (rows.size() == rowSize) {
                subMat.push_back(rows);
                rows.clear();
            }
            if (subMat.size() == rowSize) {
                matrixSep.push_back(subMat);
                subMat.clear();
            }
        }

        
        mul1 = matrixSep[0];
        mul2 = matrixSep[1];
        mul3 = matrixSep[2];
        mul4 = matrixSep[3];
        mul5 = matrixSep[4];
        mul6 = matrixSep[5];
        mul7 = matrixSep[6];

        vector<vector<int>> result;


        result = combine({mul1, mul2, mul3, mul4, mul5, mul6, mul7});

        printf("answer has %d rows and %d columns \n", result.size(), result[0].size());


        // printf("answer is : \n");

        // for (int i = 0; i < rowA; i++) {
        //     for (int j = 0; j < colB; j++) {
        //         printf("%d, ", result[i][j]);
        //     }
        //     printf("\n");
        // }

        auto end_time = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

        // Print the duration
        printf("The calculation takes %d seconds : \n", duration/1000000);
    }

    MPI_Finalize();
    return 0;
}