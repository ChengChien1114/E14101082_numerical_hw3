// E14101082_陳政謙
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>  

using namespace std;

double lagrangeInterpolation(const vector<double>& X, const vector<double>& E_X, double Guess) {
    double ANS = 0;
    int n = X.size();  
    for (int i = 0; i < n; ++i) {
        double L = 1;
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                L = L * (Guess - X[j]) / (X[i] - X[j]);
            }
        }
        ANS += E_X[i] * L;
    }
    return ANS;
}

int main() {
    vector<double> X = {0.3, 0.4, 0.5, 0.6};
    vector<double> E_X = {0.740818, 0.670320, 0.606531, 0.548812};
    double Error = 0.00001;  
    double Guess = 0.6;  
    while (true) {
        double ANS = lagrangeInterpolation(X, E_X, Guess);
        if (fabs(ANS - Guess) < Error) {
            break;
        } else {
            Guess = ANS;
        }
    }

    cout << "x: " << fixed << setprecision(8) << Guess << endl;

    return 0;
}
