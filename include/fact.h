#pragma once

#include "la.h"

namespace LA {
    
    class QR {              // A = Q R
    private:
        Matrix _q;
        Matrix _r;
    public:
        QR(){};
        QR(Matrix);
        ~QR();

        Matrix q();
        Matrix r();
        void decompose(Matrix&, Matrix&, Matrix&);
    };  

    class Chol {            // A = L L'
    private:
        Matrix _l;
        bool spd;
    public:
        Chol(Matrix);
        ~Chol();

        Matrix l();
        bool is_spd();
    };   

    class Mod_chol {        // A = L D L' (L : lower triangular)
    private:
        Matrix _l;
        Matrix _d;
        bool indef;
    public:
        Mod_chol(Matrix);
        ~Mod_chol();

        Matrix l();
        Matrix d();
        bool is_def();
    }; 

    class Eig {        // A = Q D Q' (Q : orthogonal)
    private:
        QR qr;
        Matrix _q;
        Matrix _d;
    public:
        Eig(Matrix, size_t max_iterations=100, double thresh=1e-6);
        ~Eig();

        Matrix q();
        Matrix d();
        Vector eigen_values();
        Vector eigen_vector(size_t);
    }; 

} // namespace LA
