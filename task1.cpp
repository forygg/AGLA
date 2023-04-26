#include <iostream>
#include <vector>
#include <iomanip>


#ifdef WIN32
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"
#else
#define GNUPLOT_NAME "gnuplot -persist"
#endif

using namespace std;

class ColumnVector{
    vector <double> vec;
public:
    int size;
    explicit ColumnVector(int n){
        vec.resize(n);
        size = n;
        fill(vec.begin(), vec.end(), 0);
    }
    friend istream &operator>>(istream &in, ColumnVector &v) {
        for (int i = 0; i < v.vec.size(); i++) {
            in >> v.vec[i];
        }
        return in;
    }

    friend ostream &operator<<(ostream &out, ColumnVector &v) {
        for (int i = 0; i < v.vec.size(); i++) {
            if (abs(v[i]) < 1e-10) out << fixed << setprecision(2) << 0.00 << "\n";
            else out << fixed << setprecision(4) << v[i] << "\n";
        }
        return out;
    }

    double &operator[](int i) {
        return vec[i];
    }

    ColumnVector operator*(double a) {
        ColumnVector result(size);
        for (int i = 0; i < size; i++) {
            result[i] = vec[i] * a;
        }
        return result;
    }

    ColumnVector operator*(ColumnVector v) {
        ColumnVector result(size);
        for (int i = 0; i < size; i++) {
            result[i] = vec[i] * v[i];
        }
        return result;
    }

    ColumnVector operator+(ColumnVector v) {
        ColumnVector result(size);
        for (int i = 0; i < size; i++) {
            result[i] = vec[i] + v[i];
        }
        return result;
    }

    ColumnVector operator-(ColumnVector v) {
        ColumnVector result(size);
        for (int i = 0; i < size; i++) {
            result[i] = vec[i] - v[i];
        }
        return result;
    }

    ColumnVector operator/(ColumnVector v) {
        ColumnVector result(size);
        for (int i = 0; i < size; i++) {
            result[i] = vec[i] / v[i];
        }
        return result;
    }

    ColumnVector operator/(double a) {
        ColumnVector result(size);
        for (int i = 0; i < size; i++) {
            result[i] = vec[i] / a;
        }
        return result;
    }
    double norm(){
        double res = 0;
        for (int i = 0; i < this->size; i++){
            res += this->vec[i] * this->vec[i];
        }
        return res;
    }
};

class Matrix {
protected:
    vector<vector<double>> v;
public:
    int n{};
    int m{};
    Matrix(int n, int m) {
        this->n = n;
        this->m = m;
        v = vector<vector<double>>(n, vector<double>(m, 0));
    }

    Matrix() = default;

    vector<double> &operator[](int i){
        return this->v[i];
    }
    friend istream &operator>>(istream &in, Matrix &matrix) {
        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.m; j++) {
                in >> matrix[i][j];
            }
        }
        return in;
    }

    friend ostream &operator<<(ostream &out, Matrix &matrix) {
        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.m; j++) {
                if (abs(matrix[i][j]) < 1e-10) out << 0.00 << " ";
                else out << fixed << setprecision(4) << matrix[i][j] << " ";
            }
            out << "\n";
        }
        return out;
    }

    Matrix &operator=(Matrix const &matrix) = default;

    Matrix operator+(Matrix &matrix){
        if (this->n != matrix.n || this->m != matrix.m){
            cout << "Error: the dimensional problem occurred\n";
            return {0,0};
        }
        Matrix result(matrix.n, matrix.m);
        for (int i = 0; i < result.n; i++){
            for (int j = 0; j < result.m; j++){
                result[i][j] = this->v[i][j] + matrix[i][j];
            }
        }
        return result;
    }

    Matrix operator-(Matrix &matrix){
        if (this->n != matrix.n || this->m != matrix.m){
            cout << "Error: the dimensional problem occurred\n";
            return {0,0};
        }
        Matrix result(matrix.n, matrix.m);
        for (int i = 0; i < result.v.size(); i++){
            for (int j = 0; j < result[i].size(); j++){
                result[i][j] = this->v[i][j] - matrix[i][j];
            }
        }
        return result;
    }

    Matrix operator*(Matrix &matrix){
        if (this->m != matrix.n){
            cout << "Error: the dimensional problem occurred\n";
            return {0,0};
        }
        Matrix result(this->n, matrix.m);
        for (int i = 0; i < this->n; i++){
            for (int j = 0; j < matrix.m; j++){
                double sum = 0;
                for (int k = 0; k < this->m; k++){
                    sum += this->v[i][k] * matrix.v[k][j];
                }
                result[i][j] = sum;
            }
        }
        return result;
    }
    Matrix transposed(){
        Matrix result(this->m, this->n);
        for (int i = 0; i < result.n; i++){
            for (int j = 0; j < result.m; j++){
                result[i][j] = this->v[j][i];
            }
        }
        return result;
    }


};

class SquareMatrix : public Matrix{
public:
    SquareMatrix() = default;
    explicit SquareMatrix(int n) : Matrix(n, n){};
    SquareMatrix operator+(SquareMatrix &matrix){
        if (this->n != matrix.n){
            cout << "Error: the dimensional problem occurred\n";
            return SquareMatrix(0);
        }
        Matrix *m1 = this;
        Matrix *m2 = &matrix;
        Matrix result = (*m1 + *m2);
        return *((SquareMatrix *) &result);
    }
    SquareMatrix operator-(SquareMatrix &matrix){
        if (this->n != matrix.n){
            cout << "Error: the dimensional problem occurred\n";
            return SquareMatrix(0);
        }
        Matrix *m1 = this;
        Matrix *m2 = &matrix;
        Matrix result = (*m1 - *m2);
        return *((SquareMatrix *) &result);
    }
    SquareMatrix operator*(SquareMatrix &matrix){
        if (this->n != matrix.n){
            cout << "Error: the dimensional problem occurred\n";
            return SquareMatrix(0);
        }
        Matrix *m1 = this;
        Matrix *m2 = &matrix;
        Matrix result = (*m1 * *m2);
        return *((SquareMatrix *) &result);
    }
    SquareMatrix transposed(){
        Matrix *m = this;
        Matrix result = m->transposed();
        return *((SquareMatrix *) &result);
    }
};
class IdentityMatrix : public SquareMatrix{
public:
    explicit IdentityMatrix(int n) : SquareMatrix(n){
        for (int i = 0; i < n; i++){
            this->v[i][i] = 1;
        }
    }
};
class EliminationMatrix : public SquareMatrix{
public:
    int i, j;
    EliminationMatrix(Matrix &matrix, int i, int j){
        this->i = i--;
        this->j = j--;
        IdentityMatrix im(matrix.m);
        double coefficient = matrix[i][j] / matrix[j][j];
        for (int k = 0; k < im.n; k++){
            im[i][k] -= coefficient * im[j][k];
        }
        *this = *(EliminationMatrix*)(&im);
    }
};
class PermutationMatrix : public SquareMatrix{
public:
    int i, j;
    PermutationMatrix(Matrix &matrix, int i, int j){
        this->i = i--;
        this->j = j--;
        IdentityMatrix im(matrix.m);
        swap(im[i], im[j]);
        *this = *(PermutationMatrix*)(&im);
    }
};

class AugmentedMatrix{
public:
    SquareMatrix sm;
    IdentityMatrix im = IdentityMatrix(0);
    explicit AugmentedMatrix(SquareMatrix &A){
        this->sm = A;
        im = IdentityMatrix(A.n);
    }
    SquareMatrix inverse(){
        SquareMatrix A = this->sm;
        SquareMatrix I = *(SquareMatrix*)(&this->im);
        int n = A.n;
        int step = 0;
        for (int i = 0; i < n-1; i++){
            int maxPivotIndex = i;
            for (int j = i + 1; j < n; j++){
                if (abs(A[j][i]) > abs(A[maxPivotIndex][i])){
                    maxPivotIndex = j;
                }
            }
            if (maxPivotIndex != i) {
                PermutationMatrix pm(A, i+1, maxPivotIndex + 1);
                A = pm * A;
                PermutationMatrix pmI(I, i+1, maxPivotIndex + 1);
                I = pmI * I;
            }
            for (int j = i + 1; j < n; j++){
                if (A[j][i] != 0){
                    for (int k = 0; k < I.n; k++)
                        I[j][k] -= (I[i][k] * A[j][i]/A[i][i]);
                    EliminationMatrix em(A, j+1, i+1);
                    A = em * A;
                }
            }
        }
        for (int i = n-1; i >=0; i--){
            for (int j = i-1; j >= 0; j--){
                if (A[j][i] != 0){
                    for (int k = 0; k < I.n; k++)
                        I[j][k] -= (I[i][k] * A[j][i]/A[i][i]);
                    EliminationMatrix em(A, j+1, i+1);
                    A = em * A;
                }
            }
        }
        for (int i = 0; i < n; i++) {
            for (int k = 0; k < I.n; k++)
                I[i][k] /= A[i][i];
            A[i][i] = 1;
        }
        return I;
    }
};
ColumnVector solveSystem (SquareMatrix A, ColumnVector b){
    int n = A.n;
    int step = 0;
    cout << "step #" << step << ":\n";
    cout << A << b;
    for (int i = 0; i < n-1; i++){
        int maxPivotIndex = i;
        for (int j = i + 1; j < n; j++){
            if (abs(A[j][i]) > abs(A[maxPivotIndex][i])){
                maxPivotIndex = j;
            }
        }
        if (maxPivotIndex != i) {
            PermutationMatrix pm(A, i+1, maxPivotIndex + 1);
            A = pm * A;
            swap(b[i], b[maxPivotIndex]);
            cout << "step #" << ++step << ": permutation\n";
            cout << A << b;
        }
        for (int j = i + 1; j < n; j++){
            if (A[j][i] != 0){
                b[j] -= b[i] * A[j][i]/A[i][i];
                EliminationMatrix em(A, j+1, i+1);
                A = em * A;
                cout << "step #" << ++step << ": elimination\n";
                cout << A << b;
            }
        }
    }
    for (int i = n-1; i >=0; i--){
        for (int j = i-1; j >= 0; j--){
            if (A[j][i] != 0){
                b[j] -= b[i] * A[j][i]/A[i][i];
                EliminationMatrix em(A, j+1, i+1);
                A = em * A;
                cout << "step #" << ++step << ": elimination\n";
                cout << A << b;
            }
        }
    }
    cout << "Diagonal normalization:\n";
    for (int i = 0; i < n; i++) {
        b[i] /= A[i][i];
        A[i][i] = 1;
    }
    cout << A << b;
    cout << "result:\n";
    return b;
}

double findDeterminant(SquareMatrix &A){
    int n = A.n;
    int step = 0;
    for (int i = 0; i < n; i++){
        int maxPivotIndex = i;
        for (int j = i + 1; j < n; j++){
            if (abs(A[j][i]) > abs(A[maxPivotIndex][i])){
                maxPivotIndex = j;
            }
        }
        if (maxPivotIndex != i) {
            PermutationMatrix pm(A, i+1, maxPivotIndex + 1);
            A = pm * A;
            cout << "step #" << ++step << ": permutation\n" << A;
        }
        for (int j = i + 1; j < n; j++){
            if (A[j][i] != 0){
                EliminationMatrix em(A, j+1, i+1);
                A = em * A;
                cout << "step #" << ++step << ": elimination\n" << A;
            }
        }
    }
    double det = 1;
    for (int i = 0; i < n; i++){
        det *= A[i][i];
    }
    cout << "result:\n";
    return det;
}

Matrix LeastSquareApproximation(Matrix A, ColumnVector b, int degree){
    cout << "A:\n" << A;
    Matrix A_T = A.transposed();
    Matrix temp = A_T*A;
    cout << "A_T*A:\n" << temp;
    SquareMatrix temp2 = *(SquareMatrix*)(&temp);
    AugmentedMatrix am(temp2);
    temp2 = am.inverse();
    cout << "(A_T*A)^-1:\n" << temp2;
    Matrix B(b.size, 1);
    for (int i = 0; i < b.size; i++)
        B[i][0] = b[i];
    Matrix temp3 = A_T * B;
    cout << "A_T*b:\n" << temp3;
    Matrix ans = (Matrix)temp2 * temp3;
    cout << "x~:\n" << ans;
};

int main(){
#ifdef WIN32
    FILE *pipe = _popen(GNUPLOT_NAME, "w");
#else
    FILE* pipe = popen(GNUPLOT_NAME, "w");
#endif
    if (pipe != NULL) {
    int m;
    cin >> m;
    ColumnVector t(m);
    ColumnVector b(m);
    
    fprintf(pipe, "%s\n", "plot '-' title 'Data' with points, '-' title 'Least " 
            "Squares Approximation' with lines");
        
    for (int i = 0; i < m; i++){
        cin >> t[i];
        cin >> b[i];
    }
    int degree;
    cin >> degree;
    Matrix A(m, degree + 1);
    for (int i = 0; i < m; i++){
        A[i][0] = 1;
        A[i][1] = t[i];
        for (int j = 2; j < degree + 1; j++){
            A[i][j] = A[i][j-1] * t[i];
        }
    }
    Matrix Answer = LeastSquareApproximation(A, b, degree);
    for (int i = 0; i < m; i++) {
         int x, y;
         x = t[i];
         y = b[i];
         fprintf(pipe, "%d\t%d\n", x, y);
    }
    fprintf(pipe, "%s\n", "e");
    for (double x = -10; x < 10; x+=0.1) {
        double y = 0;
        for (int i = 0; i < Answer.rows; i++)
            y += Answer.myMatrix[i][0] * pow(x, i);
        fprintf(pipe, "%f\t%f\n", x, y);
    }
        fprintf(pipe, "%s\n", "e");
        fflush(pipe);
        pclose(pipe);
    }
}

