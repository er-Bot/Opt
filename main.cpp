#include <iostream>

#include "nlo.h"

LA::Matrix a;
LA::Vector b;
double quad(Vector x) {
    return (transpose(x) * a * x + transpose(x) * b).ravel();
}

double rosenbrock(Vector x){
    return pow(1 - x(1), 2) + 100 * pow(x(0) - x(1)*x(1), 2);
}

double brent(Vector x){
    return pow(x(0) + 10, 2) + pow(x(1) + 10, 2) + exp(-x(0)*x(0) - x(1) * x(1));
}

double foo(Vector x){
    return exp(-x(0)+1) * pow(x(1) - 4, 2) + pow(x(0) - 3, 2) - 1;
};

void test_NRO(){
    try{
        RA::Function f = RA::Function(2);
        f.set_func(foo);

        double v[] = {5, 5};
        LA::Vector x = LA::Vector(v, 2);

        NLO::Newton_Rahpson_Optimizer opt = NLO::Newton_Rahpson_Optimizer(f);
        opt.set_iterations(10);
        Vector xo = opt.optimize(x, 1e-5);

        std::cout << "Minimum of f is at point : " << "\n";
        xo.print();
        std::cout << "In " << opt.nbr_iter() << " iterations." << "\n";

    } catch(LA::Exception& e) {
        std::cerr << e.get_message() << '\n';
    }
}

void test_LCG(){
    try{
        Matrix A = Matrix(3, 3);
        double Av[] = {25, 15, -5, 15, 18, 0, -5, 0, 11};
        A.from_list(Av, 9);
        double bv[] = {5, 2, 3};
        Vector b = Vector(bv, 3);

        NLO::Lin_CG_Optimizer opt = NLO::Lin_CG_Optimizer(A, b);
        double xv[] = {1, 1, 1};
        Vector x0 = Vector(xv, 3);

        opt.set_iterations(10);

        Vector xopt = opt.optimize(x0, 1e-3);
        xopt.print();
    } catch(LA::Exception& e) {
        std::cerr << e.get_message() << '\n';
    }
}

void test_FR(){
    try{
        RA::Function f = RA::Function(2);
        f.set_func(foo);

        double v[] = {5, 5};
        LA::Vector x = LA::Vector(v, 2);

        (-x).print();

        NLO::Line_Search ls = NLO::Line_Search(f);
        ls.armijo(.5, .95);

        NLO::Fletcher_Reeves_Optimizer opt = NLO::Fletcher_Reeves_Optimizer(f, ls);
        opt.set_iterations(100000);
        Vector xo = opt.optimize(x, 1e-6);

        std::cout << "Minimum of f is at point : " << "\n";
        xo.print();

        std::cout << "In " << opt.nbr_iter() << " iterations." << "\n";
        std::cout << "f.grad(xopt) = " << "\n";
        f.grad(xo).print();

    } catch(LA::Exception& e) {
        std::cerr << e.get_message() << '\n';
    }
}

void test_GD(){
    try{
        double av[] = {1,0,0,1};
        a = Matrix(2);
        a.from_list(av, 4);
        double bv[] = {0,0};
        b = Vector(bv, 2);

        double xv[] = {2,3};
        LA::Vector x = LA::Vector(xv, 2);

        RA::Function f = RA::Function(2);
        f.set_func(brent);

        NLO::Line_Search ls = NLO::Line_Search(f);
        ls.lipschitz(1e-16, .95);

        NLO::GD_Optimizer opt = NLO::GD_Optimizer(f, ls);
        opt.set_iterations(1000);
        Vector xo = opt.optimize(x, 1e-6);

        std::cout << "Minimum of f is at point : " << "\n";
        xo.print();

        std::cout << "In " << opt.nbr_iter() << " iterations." << "\n";
        std::cout << "f.grad(xopt) = " << "\n";
        f.grad(xo).print();

    } catch(LA::Exception& e) {
        std::cerr << e.get_message() << '\n';
    }
}


int main(){
    test_FR();
}