// E14101082_陳政謙
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <iomanip>
using namespace std;

void buildHermite(const vector<double>& x,
                  const vector<double>& f,
                  const vector<double>& fp,
                  vector<double>& X,
                  vector<double>& Y)
{
    int n = (int)x.size();
    int m = 2 * n;
    X.resize(m);
    Y.resize(m);
    vector<vector<double>> dd(m, vector<double>(m, 0.0));

    for(int i = 0; i < n; i++){
        X[2*i]     = x[i];
        X[2*i + 1] = x[i];
        dd[2*i][0]     = f[i];
        dd[2*i + 1][0] = f[i];
    }

    for(int i = 0; i < n; i++){
        dd[2*i][1] = fp[i];
        if(i < n - 1){
            dd[2*i + 1][1] = (f[i+1] - f[i]) / (x[i+1] - x[i]);
        }
    }

    for(int j = 2; j < m; j++){
        for(int i = 0; i < m - j; i++){
            double denom = (X[i + j] - X[i]);
            dd[i][j] = (dd[i + 1][j - 1] - dd[i][j - 1]) / denom;
        }
    }

    for(int k = 0; k < m; k++){
        Y[k] = dd[0][k];
    }
}

double evalHermite(double t, const vector<double>& X, const vector<double>& Y)
{
    int m = (int)X.size();
    double p = Y[m - 1];
    for(int i = m - 2; i >= 0; i--){
        p = Y[i] + (t - X[i]) * p;
    }
    return p;
}

double evalHermiteDeriv(double t, const vector<double>& X, const vector<double>& Y)
{
    int m = (int)X.size();
    double dP = 0.0;
    for(int k = 1; k < m; k++){
        double coeff = Y[k];
        double sumOfProducts = 0.0;
        for(int j = 0; j < k; j++){
            double prod = 1.0;
            for(int r = 0; r < k; r++){
                if(r != j) {
                    prod *= (t - X[r]);
                }
            }
            sumOfProducts += prod;
        }
        dP += coeff * sumOfProducts;
    }
    return dP;
}

double findTimeExceedSpeed(double speedLimit,
                           double left,
                           double right,
                           const vector<double>& X,
                           const vector<double>& Y)
{
    const int N = 100000;
    double h = (right - left) / N;
    double tOld = left;
    double vOld = evalHermiteDeriv(tOld, X, Y);

    for(int i = 1; i <= N; i++){
        double tNew = left + i * h;
        double vNew = evalHermiteDeriv(tNew, X, Y);
        if(vOld <= speedLimit && vNew > speedLimit){
            double a = tOld, b = tNew;
            for(int j = 0; j < 50; j++){
                double mid = 0.5 * (a + b);
                double vm = evalHermiteDeriv(mid, X, Y);
                if(vm > speedLimit) {
                    b = mid;
                } else {
                    a = mid;
                }
            }
            return 0.5 * (a + b);
        }
        tOld = tNew;
        vOld = vNew;
    }
    return numeric_limits<double>::quiet_NaN();
}

double findMaxSpeed(double left,
                    double right,
                    const vector<double>& X,
                    const vector<double>& Y,
                    double& maxT)
{
    const int N = 200000;
    double h = (right - left) / N;
    maxT = left;
    double maxV = -1e9;

    for(int i = 0; i <= N; i++){
        double curT = left + i * h;
        double curV = evalHermiteDeriv(curT, X, Y);
        if(curV > maxV){
            maxV = curV;
            maxT = curT;
        }
    }
    return maxV;
}

int main()
{
    vector<double> tData = {0, 3, 5, 8, 13};
    vector<double> sData = {0, 200, 375, 620, 990};
    vector<double> vData = {75, 77, 80, 74, 72};

    vector<double> X, Y;
    buildHermite(tData, sData, vData, X, Y);

    cout << fixed << setprecision(6);

    double pos10 = evalHermite(10.0, X, Y);
    double vel10 = evalHermiteDeriv(10.0, X, Y);
    cout << "(a) At t=10 seconds:" << endl;
    cout << "    Position H(10) = " << pos10 << " ft" << endl;
    cout << "    Speed   H'(10) = " << vel10 << " ft/s" << endl << endl;

    double speedLimit = 80.6;
    double exceedTime = findTimeExceedSpeed(speedLimit, 0.0, 13.0, X, Y);
    if(std::isnan(exceedTime)){
        cout << "(b) Speed never exceeds 80.6 ft/s in [0, 13]." << endl << endl;
    } else {
        cout << "(b) Speed first exceeds 80.6 ft/s at t = " << exceedTime << " seconds." << endl << endl;
    }

    double maxT;
    double maxV = findMaxSpeed(0.0, 13.0, X, Y, maxT);
    double maxV_mph = maxV * 0.681818;
    cout << "(c) Max speed = " << maxV << " ft/s"
         << " at t = " << maxT << " s"
         << " ( = " << maxV_mph << " mph)" << endl;

    return 0;
}
