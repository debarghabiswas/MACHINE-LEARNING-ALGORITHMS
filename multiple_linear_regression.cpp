#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

using namespace std;

class MultipleLinearRegression{
public:
    const double EPS = 1e-9;
    vector<double> b; //STORING THE COEFFICIENTS

    vector<int> dimmatrix(const vector<vector<double>>& x){
        vector<int> n(2);
        n[0] = x.size();
        n[1] = x[0].size();
        return n;
    }

    vector<vector<double>> invmatrix(const vector<vector<double>>& x){ //USING GAUSS JORDEN ELIMINATION METHOD FOR FINDING THE INVERSE MATRIX
        vector<int> n = dimmatrix(x);
        
        if(n[0] == 0 || n[1] != n[0]){
            throw runtime_error("Matrix must be non-empty and sqaure.");
        }

        vector<vector<double>> aug(n[0], vector<double>(2*n[0], 0.0));//Augmented matrix

        for(int i = 0; i < n[0]; i++){
            if(x[i].size() != n[0]){
                throw runtime_error("Matrix must be sqaure.");
            }
            for (int j=0; j<n[0]; j++){
                aug[i][j] = x[i][j];
            }
            aug[i][n[0]+i] = 1.0;
        }

        //GAUSS-JORDON ELIMINATION
        for(int i=0; i<n[0]; i++){
            int pivot = i;
            for(int row = i+1; row<n[0]; row++){
                if(fabs(aug[row][i]) > fabs(aug[pivot][i])){
                    pivot = row;
                }
            }

            if (fabs(aug[pivot][i]) < EPS){
                throw runtime_error("Maxtrix is singular and cannot be inverted.");
            }

            swap(aug[i], aug[pivot]);

            double pivotVal = aug[i][i];

            for(int j=0; j<2*n[0]; j++){
                aug[i][j] /= pivotVal;
            }

            for(int row = 0; row < n[0]; row++){
                if(row!=i){
                    double factor = aug[row][i];
                    for(int j=0; j<2*n[0]; j++){
                        aug[row][j] -= factor * aug[i][j];
                    }
                }
            }
        }

        vector<vector<double>> inv(n[0], vector<double>(n[0],0.0));
        for(int i=0; i<n[0]; i++){
            for(int j=0; j<n[0]; j++){
                inv[i][j] = aug[i][n[0]+j];
            }
        }

        return inv;
    }

    vector<vector<double>> transposematrix(vector<vector<double>> x){
        vector<int> n = dimmatrix(x);
        vector<vector<double>> xT(n[1], vector<double>(n[0], 0.0));

        for(int i=0; i<n[0]; i++){
            for(int j=0; j<n[1]; j++){
                xT[j][i] = x[i][j];
            }
        }

        return xT;
    }

    vector<vector<double>> matrixMul2_2d(vector<vector<double>> A, vector<vector<double>> B){
        int m = A.size();
        int n = A[0].size();
        int p = B[0].size();

        if (B.size() != n){
            throw invalid_argument("Number of columns in A must equal number of rows in B.");
        }

        vector<vector<double>> result(m, vector<double>(p, 0.0));

        for(int i=0; i<m; ++i){
            for(int j=0; j<p; ++j){
                for(int k=0; k<n; ++k){
                    result[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return result;
    }

    vector<double> matrixMul2_1d(vector<vector<double>> A, vector<double> B){
        int m = A.size();
        int n = A[0].size();
        if(B.size() != n){
            throw invalid_argument("Matrix column count must match vector size.");
        }
        vector<double> result(m,0.0);
        for(int i=0; i<m; ++i){
            for(int j=0; j<n; ++j){
                result[i] += A[i][j]*B[j];
            }
        }
        return result;
    }

    void fit(vector<vector<double>> dependentVar, vector<double> independentVar){
        // vector<int> n = dimmatrix(dependentVar);
        //INSERTING A COLUMN 1s TO DEPENDENT VARIABLE
        vector<vector<double>> X = dependentVar;
        vector<double> Y = independentVar;
        double newval = 1;
        for(auto& row: X){
            row.insert(row.begin(), newval);
        }

        //FOLLOWING IS FOR DEBUGGING PURPOSE
        // cout << "X matrix after adding intercept column:" << endl;
        // for(int i=0; i<X.size(); i++) {
        //     for(int j=0; j<X[0].size(); j++) {
        //         cout << X[i][j] << " ";
        //     }
        //     cout << endl;
        // }
        
        //Y = Xb + e
        //b = (Xt * X)-1 * (Xt * Y)
        vector<vector<double>> Xt = transposematrix(X);
        vector<vector<double>> XtX = matrixMul2_2d(Xt, X);

        //FOLLOWING IS FOR DEBUGGING PURPOSE
        // cout << "XtX matrix:" << endl;
        // for(int i=0; i<XtX.size(); i++){
        //     for(int j=0; j<XtX[0].size(); j++) {
        //         cout << XtX[i][j] << " ";
        //     }
        //     cout << endl;
        // }

        vector<double> XtY = matrixMul2_1d(Xt, Y);
        vector<vector<double>> A_inv = invmatrix(XtX);
        b = matrixMul2_1d(A_inv, XtY);
    }

    double predict(const vector<double>& x){
        if(b.empty()) {
            throw runtime_error("Model not fitted yet.");
        }
        double sum = b[0]; // Intercept
        for(int i=0; i<x.size(); i++){
            sum += x[i] * b[i+1];
        }
        return sum;
    }
};

int main(){
    vector<vector<double>> x = { //DEPENDENT VARIABLE
        {2,3,5},
        {1,0,2},
        {3,2,4},
        {5,1,3},
        {4,4,1},
        {6,2,0},
        {3,3,3},
        {2,5,1},
        {4,1,4},
        {5,3,2}
    };

    vector<double> y = {28,11,25,29,30,31,27,25,31,33}; //INDEPENDENT VARIABLE

    MultipleLinearRegression model;
    model.fit(x,y);

    vector<double> pred = {5,3,2};

    cout << "Predicted Value: " << model.predict(pred);

    return 0;
}