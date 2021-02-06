#pragma once

#include "ra.h"

using namespace RA;
using namespace LA;

namespace NLO {
    
    /***  One Dimentional Line Searcher  ***/
    class Line_Search{
    public:
        Line_Search(Function _f) : f(_f){
            constant(1);
        }
        Line_Search() {
            constant(1);
        }
        ~Line_Search(){}

        typedef enum {
            Constant,
            Scaled_Constant,
            Armijo,
            Wolfe,
            Strong_Wolfe,
            Lipschitz,
            Fibonacci
        } Algorithm;

        void constant(double alpha);
        void scaled_constant(double alpha);
        void armijo(double c, double rho);
        void wolfe(double c1, double c2, double rho);
        void strong_wolfe(double c1, double c2, double rho);
        void lipschitz(double eps, double rho);
        void fibonacci(double lb, double ub, double precision);

        Algorithm type(){return _type;}

        double get();
        double get(Vector p);
        double get(Vector x, Vector p);
    private:
        Function f;
        double c1;
        double c2;
        double rho;
        Algorithm _type;

        inline static double get_armijo(Function f, Vector x, Vector p, double c1 = .5, double rho=0.95); 
        inline static double get_wolfe(Function f, Vector x, Vector p, double c1 = .25, double c2 = .75, double rho=0.95); 
        inline static double get_strong_wolfe(Function f, Vector x, Vector p, double c1 = .25, double c2 = .75, double rho=0.95); 
        inline static double get_lipschitz(Function f, Vector x, Vector p, double eps, double rho);
        inline static double get_fibonacci(Function f, Vector x, Vector p, double lb, double ub, double prec);
    };
    
    /***  Non Linear Optimizers  ***/
    class Optimizer{
    protected:
        //double step;
        double precision;
        size_t max_iter = INT64_MAX;
        size_t iter;
    public:
        Optimizer(){}
        ~Optimizer(){}

        void set_iterations(size_t);
        size_t nbr_iter();
        virtual Vector optimize(Vector x, double prec=1e-6) = 0;
    };

    class GD_Optimizer : public Optimizer {
    private:
        Function f;
        Line_Search step;
    public:
        GD_Optimizer(Function _f) : f(_f), step(Line_Search(_f)){}
        GD_Optimizer(Function _f, Line_Search ls) : f(_f), step(ls){}
        ~GD_Optimizer(){}

        Vector optimize(Vector, double prec=1e-6) override;
    };

    class Newton_Rahpson_Optimizer : public Optimizer {
    private:
        Function f;
    public:
        Newton_Rahpson_Optimizer(Function _f) : f(_f){}
        ~Newton_Rahpson_Optimizer(){}

        Vector optimize(Vector, double prec=1e-6) override;
    };

    class Lin_CG_Optimizer : public Optimizer {
    private:
        Matrix A;
        Vector b;
    public:
        Lin_CG_Optimizer(Matrix _A, Vector _b) : A(_A), b(_b) {}
        ~Lin_CG_Optimizer(){}

        Vector optimize(Vector, double prec=1e-6) override;
    };
    
    class Fletcher_Reeves_Optimizer : public Optimizer {
    private:
        Function f;
        Line_Search step;
    public:
        Fletcher_Reeves_Optimizer(Function _f) : f(_f), step(Line_Search(_f)){}
        Fletcher_Reeves_Optimizer(Function _f, Line_Search ls) : f(_f), step(ls){}
        ~Fletcher_Reeves_Optimizer(){}

        Vector optimize(Vector, double prec=1e-6) override;
    };

} // namespace NLO
