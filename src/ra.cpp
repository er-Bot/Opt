#include "ra.h"

namespace RA {
    
    Function::Function(size_t n){
        order = n;
    }
    
    Function::~Function(){}

    double Function::operator() (Vector x) const{
        return f(x);
    }

    void Function::set_func(realfunc func){
        f = func;
    }

    double Function::partial(Vector x, size_t comp){
        Vector nx = x;
        Vector px = x;

        nx(comp) += EPSILON;
        px(comp) -= EPSILON;

        return (f(nx) - f(px)) / (2 * EPSILON);
    }

    Vector Function::grad(Vector x){
        Vector g = Vector(x.length());

        for(size_t i = 0; i < x.length(); i++)
            g(i) = partial(x, i);

        return g;
    }

    Matrix Function::hess(Vector x){
        size_t n = x.length();
        Matrix h = Matrix(n, n);
        double fx = f(x);
        double ep2 = EPSILON * EPSILON;
        double ep3 = 4 * EPSILON * EPSILON;

        for(size_t i = 0; i < n; i++){
            Vector x1 = x;
            Vector x2 = x;
            x1(i) -= EPSILON;
            x2(i) += EPSILON;

            h(i,i) = (f(x2) - 2 * fx + f(x1)) / ep2;

            size_t j = i + 1;
            while(j < n){
                x1(j) -= EPSILON;
                x2(j) += EPSILON;
                double v4 = f(x1);
                double v1 = f(x2);

                x1(j) += 2 * EPSILON;
                x2(j) -= 2 * EPSILON;
                double v2 = f(x1);
                double v3 = f(x2);

                h(i,j) = (v1 + v4 - v2 - v3)/ep3;
                h(j,i) = h(i,j);

                x1(j) = x(j);
                x2(j) = x(j);
                j++;
            }
        }
        
        return h;
    }

} // namespace RA
