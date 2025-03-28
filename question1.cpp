#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

long long factorial(int n) {
    long long fact = 1;
    for (int i = 2; i <= n; ++i) {
        fact *= i;
    }
    return fact;
}
double lagrangeInterpolation(const vector<double>& x_vals, const vector<double>& y_vals, double x, int degree) {
    double result = 0.0;
    
    for (int i = 0; i <= degree; ++i) {
        double term = y_vals[i];
        for (int j = 0; j <= degree; ++j) {
            if (j != i) {
                term *= (x - x_vals[j]) / (x_vals[i] - x_vals[j]);
            }
        }
        result += term;
    }
    return result;
}

double errorBound(const vector<double>& x_vals, int degree, double x, double max_derivative) {
    double product = 1.0;
    for (int i = 0; i <= degree; ++i) {
        product *= (x - x_vals[i]);
    }
    return max_derivative * fabs(product) / factorial(degree + 1);
}

int main() {
    vector<double> x_vals = {0.698, 0.733, 0.768, 0.803};
    vector<double> y_vals = {0.7661, 0.7432, 0.7193, 0.6946};
    double target = 0.750;
    double max_derivative = sin(0.750); 
    // Degree 1
    double y_1 = lagrangeInterpolation(x_vals, y_vals, target, 1);
    cout << "Degree 1 Approximation: " << y_1 << endl;
    double error_1 = errorBound(x_vals, 1, target, max_derivative);
    cout << "Degree 1 Error Bound: " << error_1 << endl;
    // Degree 2
    double y_2 = lagrangeInterpolation(x_vals, y_vals, target, 2);
    cout << "Degree 2 Approximation: " << y_2 << endl;
    double error_2 = errorBound(x_vals, 2, target, max_derivative);
    cout << "Degree 2 Error Bound: " << error_2 << endl;
    // Degree 3 
    double y_3 = lagrangeInterpolation(x_vals, y_vals, target, 3);
    cout << "Degree 3 Approximation: " << y_3 << endl;
    double error_3 = errorBound(x_vals, 3, target, max_derivative);
    cout << "Degree 3 Error Bound: " << error_3 << endl;
    // Degree 4 
    x_vals.push_back(0.750);
    y_vals.push_back(cos(0.750));  
    double y_4 = lagrangeInterpolation(x_vals, y_vals, target, 4);
    cout << "Degree 4 Approximation: " << y_4 << endl;
    double error_4 = errorBound(x_vals, 4, target, max_derivative);
    cout << "Degree 4 Error Bound: " << error_4 << endl;

    return 0;
}
