#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

// Function to check if a matrix is diagonally dominant
bool isDiagonallyDominant(vector<vector<double>> &matrix) {
    for (int i = 0; i < 3; i++) {
        double sum = 0;
        for (int j = 0; j < 3; j++) {
            if (i != j) sum += fabs(matrix[i][j]);
        }
        if (fabs(matrix[i][i]) <= sum) return false;
    }
    return true;
}

// Jacobi method
void jacobiMethod(vector<vector<double>> &A, vector<double> &b, vector<double> &x, int maxIter = 100, double tol = 1e-6) {
    vector<double> x_new(3);
    int iter = 0;

    cout << "\nJacobi Iterations:\n";
    while (iter < maxIter) {
        for (int i = 0; i < 3; i++) {
            double sum = 0;
            for (int j = 0; j < 3; j++) {
                if (i != j) sum += A[i][j] * x[j];
            }
            x_new[i] = (b[i] - sum) / A[i][i];
        }
        
        double error = 0;
        for (int i = 0; i < 3; i++) {
            error += fabs(x_new[i] - x[i]);
            x[i] = x_new[i];
        }
        
        cout << "Iteration " << iter + 1 << ": ";
        for (int i = 0; i < 3; i++) cout << "x" << i + 1 << " = " << x[i] << " ";
        cout << "Error: " << error << endl;
        
        if (error < tol) break;
        iter++;
    }
}

// Gauss-Seidel method
void gaussSeidelMethod(vector<vector<double>> &A, vector<double> &b, vector<double> &x, int maxIter = 100, double tol = 1e-6) {
    int iter = 0;
    
    cout << "\nGauss-Seidel Iterations:\n";
    while (iter < maxIter) {
        double error = 0;
        for (int i = 0; i < 3; i++) {
            double sum = 0;
            for (int j = 0; j < 3; j++) {
                if (i != j) sum += A[i][j] * x[j];
            }
            double x_new = (b[i] - sum) / A[i][i];
            error += fabs(x_new - x[i]);
            x[i] = x_new;
        }
        
        cout << "Iteration " << iter + 1 << ": ";
        for (int i = 0; i < 3; i++) cout << "x" << i + 1 << " = " << x[i] << " ";
        cout << "Error: " << error << endl;
        
        if (error < tol) break;
        iter++;
    }
}

int main() {
    vector<vector<double>> A(3, vector<double>(3));
    vector<double> b(3), x(3, 0.0);

    cout << "Enter the coefficients of the 3 equations (A matrix):\n";
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cout << "A[" << i + 1 << "][" << j + 1 << "]: ";
            cin >> A[i][j];
        }
    }

    cout << "Enter the constants (b vector):\n";
    for (int i = 0; i < 3; i++) {
        cout << "b[" << i + 1 << "]: ";
        cin >> b[i];
    }

    if (!isDiagonallyDominant(A)) {
        cout << "The matrix is not diagonally dominant. The Jacobi and Gauss-Seidel methods may not converge.\n";
        return 1;
    }

    cout << "Initial guesses for x1, x2, x3: ";
    for (int i = 0; i < 3; i++) cin >> x[i];

    // Solve using Jacobi Method
    jacobiMethod(A, b, x);

    // Reset initial guess for Gauss-Seidel
    fill(x.begin(), x.end(), 0.0);

    cout << "\nEnter initial guesses for Gauss-Seidel method (x1, x2, x3): ";
    for (int i = 0; i < 3; i++) cin >> x[i];

    // Solve using Gauss-Seidel Method
    gaussSeidelMethod(A, b, x);

    // Wait for user to press Enter before closing
    cout << "Press Enter to exit...";
    cin.ignore();  // Ignore any leftover newline character
    cin.get();

    return 0;
}
