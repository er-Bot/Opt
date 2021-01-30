#include <iostream>

#include "nlo.h"

double foo(Vector x){
    return exp(-x(0)+1) * pow(x(1) - 4, 2) + pow(x(0) - 3, 2) - 1;
};

int main(){

    try {
        RA::Function f = RA::Function(2);
        f.set_func(foo);
        double v[] = {5, 5};

        LA::Vector x = LA::Vector(v, 2);

        NLO::Newton_Rahpson_Optimizer opt = NLO::Newton_Rahpson_Optimizer(f, 1e-5);
        opt.set_iterations(10);
        Vector xo = opt.optimize(x);

        std::cout << "Minimum of Fun is at "<< std::endl;
        xo.print();
        
    } catch(LA::Exception& e){
        std::cout << e.get_message() << std::endl;
    }
    
}