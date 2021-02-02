#include "nlo.h"
#include <iostream>  // remove

namespace NLO {
    
    void Line_Search::constant(double alpha){
        c1 = alpha;
        _type = Constant;
    }

    void Line_Search::scaled_constant(double alpha){
        c1 = alpha;
        _type = Scaled_Constant;
    }

    void Line_Search::armijo(double c1 = .5, double rho=0.95){
        this->c1 = c1;
        this->rho = rho;
        _type = Armijo;
    }

    void Line_Search::wolfe(double c1 = .5, double c2 = .1, double rho=0.95){
        this->c1 = c1;
        this->c1 = c2;
        this->rho = rho;
        _type = Wolfe;
    }

    void Line_Search::strong_wolfe(double c1 = .5, double c2 = .1, double rho=0.95){
        this->c1 = c1;
        this->c1 = c2;
        this->rho = rho;
        _type = Wolfe;
    }

    void Line_Search::lipschitz(double eps, double rho=0.95){
        this->c1 = eps;
        this->rho = rho;
        _type = Lipschitz;
    }

    double Line_Search::get(){
        if(_type != Constant)
            throw NotImplementedError("The used line search isn't of type constant!");
        return c1;
    }

    double Line_Search::get(Vector p){
        if(_type != Scaled_Constant)
            throw NotImplementedError("The used line search isn't of type scaled constant!");
        return c1 / p.norm();
    }

    double Line_Search::get(Vector x, Vector p){
        switch (_type){
        case Armijo:
            return get_armijo(f, x, p, c1, rho);
        case Wolfe:
            return get_wolfe(f, x, p, c1, c2, rho);
        case Strong_Wolfe:
            return get_wolfe(f, x, p, c1, c2, rho);
        case Lipschitz:
            return get_lipschitz(f, x, p, c1, rho);
        default:
            throw NotImplementedError("Try to specify what step size line search by setting it's hyperparameters!");
        }
    }

    double Line_Search::get_armijo(Function f, Vector x, Vector p, double c, double rho){
        if(rho <= 0 || rho >= 1){
            printf("Backtracking: 'rho' must be in (0,1)!, reset to 0.95\n");
            rho = .95;
        }
        if(c <= 0 || c >= 1){
            printf("Backtracking: 'c' must be in (0,1)!, reset to 0.5\n");
            c = .5;
        }

        double alpha = 1;
        // armijo
        while(f(x + alpha * p) > f(x) + c * alpha * f.grad(x).dot(p)){
            alpha = rho * alpha;
        }
        if (alpha < 1e-16) throw ValueError("Step size too small!");
        
        return alpha;
    }

    double Line_Search::get_wolfe(Function f, Vector x, Vector p, double c1, double c2, double rho){
        double lb = 0;
        double ub = __DBL_MAX__;

        if (rho >= 1 || rho <= 0)
            throw IncompatibleArguments("rho must be in (0,1)!");

        if(c2 >= 1 || c1 <= 0 || c2 >= c1)
            throw IncompatibleArguments("wolfe is applicable only for 0 < c1 < c2 < 1!");

        Vector g = f.grad(x);

        double alpha = 1;
        while (true){
            if(f(x + alpha * p) > f(x) + c1 * alpha * (g * p)){
                ub = alpha;
                alpha = .5 * (lb + ub);
            } else if(p * f.grad(x + alpha * p) < c2 * p * g){
                lb = alpha;
                if (_is_inf(ub))    alpha = 2 * lb;
                else                alpha = .5 * (lb + ub);
            } else break;
        }
        if (alpha < 1e-16) throw ValueError("Step size too small!");

        return alpha;
    }

    double Line_Search::get_strong_wolfe(Function f, Vector x, Vector p, double c1, double c2, double rho){
        double lb = 0;
        double ub = __DBL_MAX__;

        if (rho >= 1 || rho <= 0)
            throw IncompatibleArguments("rho must be in (0,1)!");

        if(c2 >= 1 || c1 <= 0 || c2 >= c1)
            throw IncompatibleArguments("wolfe is applicable only for 0 < c1 < c2 < 1!");

        Vector g = f.grad(x);

        double alpha = 1;
        while (true){
            if(f(x + alpha * p) > f(x) + c1 * alpha * (g * p)){
                ub = alpha;
                alpha = .5 * (lb + ub);
            } else if(_abs(p * f.grad(x + alpha * p)) > c2 * p * g){
                lb = alpha;
                if (_is_inf(ub))    alpha = 2 * lb;
                else                alpha = .5 * (lb + ub);
            } else break;
        }
        if (alpha < 1e-16) throw ValueError("Step size too small!");

        return alpha;
    }

    double Line_Search::get_lipschitz(Function f, Vector x, Vector p, double eps, double rho){
        if (rho >= 1 || rho <= 0)
            throw IncompatibleArguments("rho must be in (0,1)!");

        Vector g = f.grad(x);
        double alpha = 1;
        while (true){
            std::cout << "Lip: " << f(x + alpha * p) - f(x) << " > " << alpha * g * p << " + " << pow(alpha * p.norm(), 2) / (2 * alpha) << " + " << eps << "(= " << alpha * g * p + pow(alpha * p.norm(), 2) / (2 * alpha) + eps << ")\n";
            if(f(x + alpha * p) - f(x) > alpha * g * p + pow(alpha * p.norm(), 2) / (2 * alpha) + eps){
                alpha *= rho;
            } else break;
        }
        if (alpha < 1e-16) throw ValueError("Step size too small!");

        return alpha;
    }



    void Optimizer::set_iterations(size_t n){
        max_iter = n;
    }

    size_t Optimizer::nbr_iter(){
        return iter;
    }

    Vector GD_Optimizer::optimize(Vector x, double prec) {
        
        precision = prec;
        double alpha = 1;

        Vector g = f.grad(x);
        x.print();
        for(iter = 1; iter  < max_iter; iter++){

            if(g.norm() <= precision)
                break;

            Vector d = - g;

            if (step.type() == Line_Search::Constant)
                alpha = step.get();
            else
                alpha = step.get(x, d);
                

            printf("Iter %lu : grad_norm %.15f : step size %.18f\n", iter, g.norm(), alpha);
            f.grad(x).print();

            x += alpha * d;
            g = f.grad(x);
        }

        return x;
    }

    Vector Newton_Rahpson_Optimizer::optimize(Vector x, double prec) {
        
        precision = prec;

        if(f.grad(x).norm() < precision) 
            return x;
        
        for(iter = 1; iter  < max_iter; iter++){
            Vector d = f.hess(x).inv() * f.grad(x);
            x -= d;
            
            if(d.norm() < precision){ 
                break;
            }
        }

        return x;
    }

    Vector Lin_CG_Optimizer::optimize(Vector x, double prec) {
        precision = prec;

        if (!A.is_pd())
            printf("Matrix isn't positive definite, CG algorithm will diverge or a division by zero may occures.\n");
        
        Vector r = b - A * x;
        Vector p = r;

        double r_old = (transpose(r) * r).ravel();

        for(iter = 1; iter < max_iter; iter++){
            // alpha = (r' r) / (p' A p)
            double alpha = r_old / (transpose(p) * A * p).ravel();

            // x = x + alpha p
            x += alpha * p;
            
            // r_n = r_o + alpha A p
            r -= alpha * (A * p);
            
            // p = r + (r_n' r_n) / (r_o' r_o) * p
            double r_new = (transpose(r) * r).ravel();
            if(sqrt(r_new) < precision)
                break;
            
            p = r + p * (r_new / r_old);

            r_old = r_new;
        }
        return x;
    }

    /*** Revise: seems like there is an error   ***/
    Vector Fletcher_Reeves_Optimizer::optimize(Vector x, double prec) {
        precision = prec;

        Vector g_old = f.grad(x);
        Vector p = - g_old;
        
        double alpha;
        for(iter = 1; iter  < max_iter; iter++){
            if(g_old.norm() < precision)
                break;
            
            if (step.type() == Line_Search::Constant)
                alpha = step.get();
            else
                alpha = step.get(x, p);
            
            printf("Iter %lu : grad_norm %.15f : step size %.18f\n", iter, g_old.norm(), alpha);
            f.grad(x).print();
            Vector x_n = x + alpha * p;
            
            Vector g_new = f.grad(x);
            double beta = (transpose(g_new) * g_new).ravel() / (transpose(g_old) * g_old).ravel(); 
            p = -g_new + beta * p;
            g_old = g_new;
        }

        return x;
    }


} // namespace NLO
