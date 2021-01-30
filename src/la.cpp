#include "la.h"

#include <iostream>
#include "fact.h"

namespace LA {

    Vector::Vector(){
        ptr = nullptr;
        len = 0;
    }

    Vector::Vector(size_t len){
        this->len = len;
        ptr = new double[len];
    }

    Vector::Vector(const Vector& v){
        this->len = (&v)->len;
        this->ptr = new double[len];
        memcpy(this->ptr, (&v)->ptr, len * sizeof(double));
    }
    
    Vector::Vector(double* list, size_t len){
        this->len = len;
        this->ptr = new double[len];

        for(size_t i = 0; i < len; i++)
            (*this)(i) = list[i];
    }
    
    Vector::~Vector(){
        delete ptr;
    }

    void Vector::print(){
        if(len == 0){
            printf("Vector is empty!\n");
            return;
        }

        printf("Vector[%lu] : \n{", len);
        for(size_t i=0; i < len; i++){
            printf("%10.4f", (*this)(i));
            if (i < len - 1) printf("\n ");
        }
        printf("  }\n");
    }

    size_t Vector::length(){
        return len;
    }

    inline    
    double& Vector::operator() (size_t pos) {
        if (pos >= len)
            throw IndexOutOfRange("Vector ubscript out of bounds");
        return ptr[pos];
    }

    inline 
    double Vector::operator() (size_t pos) const {
        if (pos >= len)
            throw IndexOutOfRange("Vector subscript out of bounds");
        return ptr[pos];
    }

    double Vector::operator * (Vector m){
        return dot(m);
    }

    Vector Vector::operator + (Vector m){
        if(len != m.len)
            throw IncompatibleArguments("Cannot add vectors of different dimensions!");

        Vector tmp = Vector(len);

        for(size_t i = 0; i < len; i++)
                tmp(i) = (*this)(i) + m(i);

        return tmp;
    }

    Vector Vector::operator += (Vector m){
        if(m.len != len)
            throw IncompatibleArguments("Cannot add vectors of different dimensions!");

        for(size_t i = 0; i < m.len; i++)
                (*this)(i) += m(i);

        return *this;
    }

    Vector Vector::operator - (Vector m){
        if(len != m.len)
            throw IncompatibleArguments("Cannot subtract vectors of different dimensions!");

        Vector tmp = Vector(len);

        for(size_t i = 0; i < len; i++)
                tmp(i) = (*this)(i) - m(i);

        return tmp;
    }

    Vector Vector::operator -= (Vector m){
        if(m.len != len)
            throw IncompatibleArguments("Cannot subtract vectors of different dimensions!");

        for(size_t i = 0; i < m.len; i++)
                (*this)(i) -= m(i);

        return *this;
    }

    Vector Vector::operator *= (double m){
        for(size_t i = 0; i < len; i++)
            (*this)(i) *= m;
        return *this;
    }

    Vector Vector::operator /= (double m){
        for(size_t i = 0; i < len; i++)
            (*this)(i) /= m;
        return *this;
    }

    Vector Vector::operator -= (double m){
        for(size_t i = 0; i < len; i++)
            (*this)(i) -= m;
        return *this;
    }

    Vector Vector::operator += (double m){
        for(size_t i = 0; i < len; i++)
            (*this)(i) += m;
        return *this;
    }

    Vector Vector::operator + (double m){
        Vector tmp = Vector(len);

        for(size_t i = 0; i < len; i++)
                tmp(i) = (*this)(i) + m;

        return tmp;
    }

    Vector Vector::operator - (double m){
        Vector tmp = Vector(len);

        for(size_t i = 0; i < len; i++)
                tmp(i) = (*this)(i) - m;

        return tmp;
    }

    Vector Vector::operator * (double m){
        Vector tmp = Vector(len);

        for(size_t i = 0; i < len; i++)
                tmp(i) = (*this)(i) * m;

        return tmp;
    }

    Vector Vector::operator / (double m){
        Vector tmp = Vector(len);

        for(size_t i = 0; i < len; i++)
                tmp(i) = (*this)(i) / m;

        return tmp;
    }

    double Vector::dot(Vector a){
        if(len != a.len){
            printf("Dimensions aren't compatible for dot operation (Vectors)!\n");
            return -MAXFLOAT;
        }

        double s = 0;

        for(size_t i = 0; i < len; i++){
            s += a(i) * (*this)(i);
        }

        return s;
    }

    double Vector::norm(size_t n){
        if(n <= 0)
            throw IndexOutOfRange("Norm must be greater than 1");
        
        if(n == 1){
            double s = 0;
            for(size_t i = 0; i < len; i++)
                s += _abs((*this)(i));
            return s;
        }

        if(n == INFTY){
            double s = _abs((*this)(0));
            for(size_t i = 1; i < len; i++){
                s = _max(_abs((*this)(i)), s);
            }
            return s;
        }

        double s = 0;
        for(size_t i = 0; i < len; i++)
            s += pow((*this)(i), n);
        return pow(s, 1/(double)n);
    }

    double Vector::sum(){
        double s = 0;
        for(size_t i = 0; i < len; i++)
            s += (*this)(i);
        return s;
    }

    /******************************************
     *              Matrix Space              *
     ******************************************/

    Matrix::Matrix(){
        ptr = nullptr;
        row = 0;
        col = 0;
    }

    Matrix::Matrix(size_t s){
        this->row = s;
        this->col = s;
        ptr = new double[s * s];
    }

    Matrix::Matrix(size_t row, size_t col){
        this->row = row;
        this->col = col;
        ptr = new double[row * col];
    }

    Matrix::Matrix(const Matrix& m){
        this->row = (&m)->row;
        this->col = (&m)->col;
        this->ptr = new double[row * col];
        memcpy(this->ptr, (&m)->ptr, row * col * sizeof(double));
    }

    Matrix& Matrix::operator=(const Matrix& m){
        row = (&m)->row;
        col = (&m)->col;

        delete ptr;
        ptr = new double[row * col];
        memcpy(ptr, (&m)->ptr, row * col * sizeof(double));

        return *this;
    }

    Matrix::Matrix(Vector v){
        row = v.length();
        col = 1;
        this->ptr = new double[row * col];

        for(size_t i = 0; i < row; i++)
            (*this)(i, 0) = v(i);
    }
    
    Matrix::~Matrix(){
        delete ptr;
    }

    void Matrix::random(double min, double max, bool sym){
        if(!is_square() && sym)
            throw IncompatibleArguments("Trying to generate a symmetric matrix for a non square matrix!");

        for(size_t i = 0; i < this->row; i++)
            for(size_t j = 0; j < this->col; j++)
                if(sym && i > j)
                    this->ptr[i*this->col+j] = this->ptr[j*this->col+i];
                else
                    this->ptr[i*this->col+j] = _rand(min,max);
    }

    void Matrix::from_list(double* list, size_t size){
        if(size != row * col){
            printf("Not matching size\n");
            return;
        }

        for(size_t i = 0; i < this->row * this->col; i++)
            this->ptr[i] = list[i];
    }

    void Matrix::print(){
        if(row == 0 && col == 0){
            printf("Matrix is empty!\n");
            return;
        }

        printf("Matrix[%lu,%lu] :\n", row, col);
        
        printf("{{");

        for(size_t i=0; i < row; i++){
            for(size_t d = 0; d < col; d++){
                printf("%10.4f ", (*this)(i, d));
            }
            if(i < row - 1)
                printf("}\n {");
            else printf("}}\n\n");
        }
    }

    // __DEPRECATED
    double Matrix::get_elem(size_t row, size_t col){
        if(row >= this->row || col >= this->col){
            printf("Indices out of range\n");
            exit(1);
        }

        return this->ptr[row * this->col + col];
    }

    void Matrix::set_elem(size_t row, size_t col, double val){
        if(row >= this->row || col >= this->col){
            printf("Indices out of range\n");
            exit(1);
        }

        this->ptr[row * this->col + col] = val;
    }

    inline      // better use inline for small most used functions
    double& Matrix::operator() (size_t row, size_t col) {
        if (row >= this->row || col >= this->col)
            throw IndexOutOfRange("Matrix subscript out of bounds");
        return ptr[row * this->col + col];
    }

    inline 
    double Matrix::operator() (size_t row, size_t col) const {
        if (row >= this->row || col >= this->col)
            throw IndexOutOfRange("Matrix subscript out of bounds");
        return ptr[row * this->col + col];
    }

    Vector Matrix::get_col(size_t c){
        if(c >= col)
            throw IndexOutOfRange("Column index out of range");
        
        Vector v = Vector(row);

        for (size_t i = 0; i < row; i++)
            v(i) = (*this)(i, c);
        
        return v;
    }

    Vector Matrix::get_row(size_t r){
        if(r >= row)
            throw IndexOutOfRange("Row index out of range");
        
        Vector v = Vector(col);

        for (size_t i = 0; i < col; i++)
            v(i) = (*this)(r, i);
        
        return v;
    }

    Vector Matrix::get_diag(){
        if(row != col)
            throw Exception("It's not wise getting the diagonal of non square matrices!");
        
        Vector v = Vector(col);

        for (size_t i = 0; i < col; i++)
            v(i) = (*this)(i, i);
        
        return v;
    }

    void Matrix::set_col(size_t c, Vector v){
        if (v.length() != row)
            throw IncompatibleArguments("vector length doesn't match the number of rows!");
        if(c >= col)
            throw IndexOutOfRange("Column index out of range");

        for (size_t i = 0; i < row; i++)
            (*this)(i, c) = v(i);
    }

    void Matrix::set_row(size_t r, Vector v){
        if (v.length() != row)
            throw IncompatibleArguments("vector length doesn't match the number of cols!");
        if(r >= row)
            throw IndexOutOfRange("Row index out of range");

        for (size_t i = 0; i < col; i++)
            (*this)(r, i) = v(i);
    }

    size_t Matrix::rows(){
        return this->row;
    }

    size_t Matrix::cols(){
        return this->col;
    }

    Matrix Matrix::identity(size_t n){
        Matrix I = Matrix(n, n);
        for(size_t i = 0; i < n; i++)
            I(i, i) = 1;
        return I;
    }

    Matrix Matrix::ones(size_t n, size_t m){
        Matrix I = Matrix(n, m);
        for(size_t i = 0; i < n; i++)
        for(size_t j = 0; j < m; j++)
            I(i, j) = 1;
        return I;
    }

    Matrix Matrix::diag(double* b, size_t n){
        Matrix I = Matrix(n, n);
        for(size_t i = 0; i < n; i++)
            I(i, i) =  b[i];
        return I;
    }

    // O(n^3)
    size_t Matrix::rank(){
        Matrix m = *this;
        size_t rk = _min(row, col);

        for(size_t r = 0; r < rk; r++){
            if(m(r, r)){
                for(size_t c = 0; c < m.rows(); c++){
                    if(c != r){
                        double mult = (*this)(c, r) / m(r, r);
                        for(size_t i = 0; i < rk; i++)
                            m(c, i) = m(c, i) - mult * (*this)(r, i);
                    }
                }
            } else {
                int red = 1;
                for(size_t i = r + 1; i < m.rows(); i++){
                    if(m(i, r)){
                        // swap
                        for(size_t t = 0; t < rk; t++)
                            _swap(m(r, t), m(i, t))
                        red = 0;
                        break;
                    }
                }

                if(red){
                    rk--;
                    for(size_t i = 0; i < m.rows(); i++)
                        m(i, r) =  m(i, rk);
                }

                r--;
            }
        }
        
        return rk;
    }
    // O(n^2)
    Matrix Matrix::get_minor(size_t r, size_t c){
        int k = 0;
        int l = 0;
        Matrix m = Matrix(row - 1, row - 1);
        for (size_t i = 0; i < row; i++) {
            for (size_t j = 0; j < col; j++) {
                if (i != r && j != c) {
                    m(k, l) = (*this)(i, j);
                    l++;
                }
            }
            if (i != r) 
                k++;
            l = 0;
        }

        return m;
    }
    
    double Matrix::det(){
        if(!is_square()){
            printf("Cannot calculate determinant of no-squared matrices!\n");
            return -INT64_MAX;
        }
        
        if (row == 2)
            return ((*this)(0, 0) * (*this)(1, 1)) - ((*this)(0, 1) * (*this)(1, 0));
        double d = 0;
        for (size_t j = 0; j < row; j++) {
            Matrix Aa = get_minor(0, j);
            d +=  pow(-1, j) * (*this)(0, j) * Aa.det();
        }
        return d;
    }

    Matrix Matrix::inv(){
        if(!is_square()){
            printf("Cannot find inverse of no-squared matrices!\n");
            return Matrix();
        }
	
	    Matrix A = Matrix(row, row);

        if(row == 1){
            A(0, 0) = 1/(*this)(0, 0);
            return A;
        }

        float d = det();
        if (d == 0){
            printf("Matrix isn't invertible\n");
            return Matrix();
        }

        if (row == 2) {
            A(0, 0) = (*this)(1, 1) / d;
            A(0, 1) = -(*this)(0, 1) / d;
            A(1, 0) = -(*this)(1, 0) / d;
            A(1, 1) = (*this)(0, 0) / d;
            return A;
        }

        for (size_t i = 0; i < row; i++) {
            for (size_t j = 0; j < row; j++) {
                Matrix Aa = get_minor(i, j);
                A(i, j) = pow((-1), i+j) * Aa.det() / d;
            }
        }

        return A.T();	
    }

    double Matrix::norm(){
        double nr = -1;
        for(size_t j = 0; j < col; j++){
            double s = 0;
            for(size_t i = 0; i < row; i++)
                s += _abs((*this)(i, j));
            nr = _max(nr, s);
        }

        return nr;
    }

    double Matrix::norm_infty(){
        double nr = -1;
        for(size_t i = 0; i < row; i++){
            double s = 0;
            for(size_t j = 0; j < col; j++)
                s += _abs((*this)(i, j));
            nr = _max(nr, s);
        }

        return nr;
    }

    double Matrix::norm_frobinus(){
        double s = 0;
        for(size_t i = 0; i < row; i++)
            for(size_t j = 0; j < col; j++)
                s += pow((*this)(i, j), 2);

        return sqrt(s);
    }

    Matrix Matrix::T(){
        Matrix at = Matrix(this->col, this->row);

        for(size_t r = 0; r < this->row; r++)
            for(size_t c = 0; c < this->col; c++){
                at(c, r) = (*this)(r, c);
            }

        return at;
    }

    Matrix Matrix::operator + (Matrix m){
        if(m.rows() != row || m.cols() != col){
            printf("Cannot add matrices of different dimensions!\n");
            return Matrix();
        }

        Matrix tmp = Matrix(row, col);

        for(size_t i = 0; i < m.row; i++)
            for(size_t j = 0; j < m.col; j++)
                tmp(i, j) = (*this)(i, j) + m(i, j);

        return tmp;
    }

    Matrix Matrix::operator - (Matrix m){
        if(m.rows() != row || m.cols() != col){
            printf("Cannot subtract matrices of different dimensions!\n");
            return Matrix();
        }

        Matrix tmp = Matrix(row, col);

        for(size_t i = 0; i < m.row; i++)
            for(size_t j = 0; j < m.col; j++)
                tmp(i, j) = (*this)(i, j) - m(i, j);

        return tmp;
    }

    Matrix Matrix::operator * (Matrix m){
        return dot(m);
    }

    Vector Matrix::operator*(Vector v){
        if(col != v.length())
            throw IncompatibleArguments("Incmpatible dimensions for Matrix-Vector multiplication!");

        Vector o = Vector(row);

        for(size_t i = 0; i < row; i++)
            o(i) = get_row(i) * v;

        return o;
    }

    Matrix Matrix::operator / (Matrix m){
        return m.inv().dot(*this);
    }

    Matrix Matrix::operator + (double s){
        Matrix tmp = Matrix(row, col);

        for(size_t i = 0; i < row; i++)
            for(size_t j = 0; j < col; j++)
                tmp(i, j) = (*this)(i, j) + s;

        return tmp;
    }

    Matrix Matrix::operator - (double s){
        Matrix tmp = Matrix(row, col);

        for(size_t i = 0; i < row; i++)
            for(size_t j = 0; j < col; j++)
                tmp(i, j) = (*this)(i, j) - s;

        return tmp;
    }

    Matrix Matrix::operator * (double s){
        Matrix tmp = Matrix(row, col);

        for(size_t i = 0; i < row; i++)
            for(size_t j = 0; j < col; j++)
                tmp(i, j) = (*this)(i, j) * s;

        return tmp;
    }

    Matrix Matrix::operator / (double s){
        Matrix tmp = Matrix(row, col);

        for(size_t i = 0; i < row; i++)
            for(size_t j = 0; j < col; j++)
                tmp(i, j) = (*this)(i, j) / s;

        return tmp;
    }

    Matrix Matrix::operator += (Matrix m){
        if(m.rows() != row || m.cols() != col){
            printf("Cannot add matrices of different dimensions!\n");
            return *this;
        }

        for(size_t i = 0; i < m.row; i++)
            for(size_t j = 0; j < m.col; j++)
                (*this)(i, j) += m(i, j);

        return *this;
    }

    Matrix Matrix::operator -= (Matrix m){
        if(m.rows() != row || m.cols() != col){
            printf("Cannot subtract matrices of different dimensions!\n");
            return *this;
        }

        for(size_t i = 0; i < row; i++)
            for(size_t j = 0; j < col; j++)
                (*this)(i, j) -= m(i, j);

        return *this;
    }

    Matrix Matrix::operator += (double s){
        for(size_t i = 0; i < row; i++)
            for(size_t j = 0; j < col; j++)
                (*this)(i, j) += s;

        return *this;
    }

    Matrix Matrix::operator -= (double s){
        for(size_t i = 0; i < row; i++)
            for(size_t j = 0; j < col; j++)
                (*this)(i, j) -= s;

        return *this;
    }

    Matrix Matrix::operator *= (double s){
        for(size_t i = 0; i < row; i++)
            for(size_t j = 0; j < col; j++)
                (*this)(i, j) *= s;

        return *this;
    }

    Matrix Matrix::operator /= (double s){
        for(size_t i = 0; i < row; i++)
            for(size_t j = 0; j < col; j++)
                (*this)(i, j) /= s;

        return *this;
    }

    Matrix Matrix::dot(Matrix a){
        if(cols() != a.rows()){
            printf("Dimensions aren't compatible for multiplication!\n");
            return *this;
        }

        Matrix tmp = Matrix(row, a.cols());

        for(size_t i = 0; i < tmp.row; i++){
            for(size_t j = 0; j < tmp.col; j++){
                double s = 0;
                for(size_t k = 0; k < col; k++){
                    s += (*this)(i, k) * a(k, j);
                }
                tmp(i, j) = s;
            }
        }
        return tmp;
    }

    bool Matrix::compare(Matrix o, double thresh){
        if (row != o.row || col != o.col)
            throw IncompatibleArguments("Matrices dimensions aren't compatible!");

        for (size_t i = 0; i < row; i++) {
            for (size_t j = 0; j < col; j++) {
            double difference = _abs((*this)(i,j) - o(i,j));
            if (difference > thresh) {
                return false;
            }
            }
        }
        return true;
    }

    bool Matrix::is_square(){
        if (row == col)
            return true;
        return false;
    }


    bool Matrix::is_symmetric(){
        if(!is_square()) 
            return false;
            
        for(size_t i = 0; i < row; i++)
            for(size_t j = 0; j < col; j++)
                if(i != j && (*this)(i, j) != (*this)(j, i))
                    return false;
        
        return true;
    }

    bool Matrix::is_diag(){
        if(!is_square()) 
            return false;
            
        for(size_t i = 0; i < row; i++)
            for(size_t j = 0; j < col; j++)
                if(i != j && (*this)(i, j) != 0)
                    return false;
        
        return true;
    }

    bool Matrix::is_null(){   
        for(size_t i = 0; i < row; i++)
            for(size_t j = 0; j < col; j++)
                if((*this)(i, j) != 0)
                    return false;
        
        return true;
    }

    bool Matrix::is_pd(){
        Chol c = Chol(*this);
        if(c.is_spd() && rank() == row)
            return true;
        return false;
    }

} // namespace LA
