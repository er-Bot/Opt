#include "fact.h"
#include <iostream>  // remove
namespace LA {

    // Not correct algorithm
    QR::QR(Matrix A){
        if(!A.is_symmetric())
            throw IncompatibleArguments("Cannot apply QR on non-symetic matrices!");

        _q = Matrix(A.rows(), A.cols());
        _r = Matrix(A.rows(), A.cols());

        decompose(A, _q, _r);
    }

    void QR::decompose(Matrix& A, Matrix& Q, Matrix& R){
        for (size_t i = 0; i < A.rows(); i++) {
            Vector u_i = A.get_col(i);
            Vector q_i = u_i;
            for (size_t j = 0; j < i; j++) {
                Vector q_j = Q.get_row(j);
                R(j, i) = q_j * u_i;
                Vector temp = q_j * R(j, i);
                q_i -= temp;
            }
            R(i, i) = q_i.norm();
            q_i /= R(i, i);
            Q.set_row(i, q_i);
        }
        Q = Q.T();
    }

    QR::~QR(){}

    Matrix QR::q(){
        return _q;
    }

    Matrix QR::r(){
        return _r;
    }

    Chol::Chol(Matrix A){
        if(!A.is_square()){
            printf("Can't perfor cholesky decomposition of a non-square matrix!\n");
            return;
        }
    
        _l = Matrix(A.rows(), A.cols());
        spd = true;

        for (size_t j = 0; j < A.rows(); j++) {
            double s = A(j,j);
            for (size_t k = 0; k < j; k++)
                s -= _l(j,k) * _l(j, k);
            if (s  <= 0) {
                printf("ERROR: non-positive definite matrix!\n");
                spd = false;
                return;
            }
            _l(j,j) = sqrt(s);
            
            for (size_t i = j+1; i < A.rows(); i++) {
                s = A(i,j);
                for (size_t k = 0; k < i; k++)
                    s -= _l(i, k) * _l(j, k);
                _l(i,j) = s / _l(j,j);
            }
        }

    }

    Chol::~Chol(){}

    Matrix Chol::l(){
        return _l;
    }

    bool Chol::is_spd(){
        return spd;
    }

    // Reference: Golub and Van Loan, "Matrix Computations", second edition, p 137.  
    Mod_chol::Mod_chol(Matrix A){
        if(!A.is_square()){
            printf("Can't perform cholesky decomposition of a non-square matrix!\n");
            return;
        }
    
        _l = Matrix(A.rows(), A.cols());
        _d = Matrix(A.rows(), A.cols());
        indef = false;

        for (size_t j = 0; j < A.rows(); j++) {
            double sum = A(j,j);
            for (size_t s = 0; s < j; s++)
                sum -= _d(s,s) * _l(j,s) * _l(j,s);
            if (sum == 0) {
                printf("ERROR: indefinite matrix!\n");
                indef = true;
                return;
            }
            _d(j,j) = sum;
            
            for (size_t i = j+1; i < A.rows(); i++) {
                sum = A(i,j);
                for (size_t s = 0; s < j; s++)
                    sum -= _d(s,s) * _l(i, s) * _l(j, s);
                _l(i,j) = sum / _d(j,j);
            }
        }

        _l += Matrix::identity(_l.rows());
    }

    Mod_chol::~Mod_chol(){}

    Matrix Mod_chol::l(){
        return _l;
    }

    Matrix Mod_chol::d(){
        return _d;
    }

    bool Mod_chol::is_def(){
        return !indef;
    }

    Eig::Eig(Matrix A, size_t max_iterations, double thresh){
        if(!A.is_square())
            throw IncompatibleArguments("Can't perform eigenvalue decomposition of a non-square matrix!\n");
    
        qr = QR(A);
        Matrix R = qr.r();
        Matrix Q = qr.q();

        _d = Matrix(A.rows(), A.cols());
        _q = Q;

        Matrix prev;

        size_t iterations = 0;
        do {
            prev = _d;
            _d = R * Q;
            qr.decompose(_d, Q, R);
            _q = _q * Q;
            iterations++;
        } while (iterations < max_iterations && !prev.compare(_d, thresh));
    }

    Eig::~Eig(){}

    Matrix Eig::q(){
        return _q;
    }

    Matrix Eig::d(){
        return _d;
    }

    Vector Eig::eigen_values(){
        return _d.get_diag();
    }

    Vector Eig::eigen_vector(size_t i){
        return _q.get_col(i);
    }

} // namespace LA
