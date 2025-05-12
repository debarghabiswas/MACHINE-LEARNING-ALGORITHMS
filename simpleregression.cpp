#include <iostream>

using namespace std;

class SimpleRegression{ //y = b0 + (b1*x)
public:
    double b0,b1;
    double sumofprod(double a[], double b[], int size){
        int n = size;
        double sum=0;
        for(int i=0; i<n; i++){
            sum += (a[i]*b[i]);
        }
        return sum;
    }
    double sumofx(double a[], int size){
        int n = size;
        double sum=0;
        for(int i=0; i<n; i++){
            sum += a[i];
        }
        return sum;
    }
    double sumofy(double a[], int size){
        int n = size;
        double sum=0;
        for(int i=0; i<n; i++){
            sum+=a[i];
        }
        return sum;
    }
    double sumofx2(double a[], int size){
        int n = size;
        double sum=0;
        for(int i=0; i<n; i++){
            sum += (a[i]*a[i]);
        }
        return sum;
    }
    double mean(double a[], int size){
        int n = size;
        double sum=0;
        for(int i=0; i<n; i++){
            sum += a[i];
        }
        return (sum/n);
    }
    void fit(double dependentVar[], double independentVar[], int size){
        int n = size; //number of observations
        double x[n];
        for (int i = 0; i < n; i++) {
            x[i] = dependentVar[i];
        }
        double y[n];
        for(int i=0; i<n; i++){
            y[i] = independentVar[i];
        }

        //for finding the regression coefficient b1
        double numerator_b1, denomenator_b1;
        numerator_b1 = (n*(sumofprod(x,y,n))) - ((sumofx(x,n))*(sumofy(y,n)));
        denomenator_b1 = (n*(sumofx2(x,n))) - (sumofx(x,n) * sumofx(x,n));
        b1 = numerator_b1/denomenator_b1;

        //for finding the intercept b0
        b0 = mean(y,n) - (b1*mean(x,n));
    }

    double predict(double x){
        return (b0+(b1*x));
    }
};

int main(){
    double x[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};//dependent variable
    double y[10] = {5, 7, 9, 11, 13, 15, 17, 19, 21, 23}; //independent variable

    SimpleRegression ob;
    ob.fit(x,y,10);

    cout << "Predicted Value: " << ob.predict(9.0);

    return 0;
}