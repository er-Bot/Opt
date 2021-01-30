#include "nlo.h"
#include <iostream>  // remove

namespace NLO {
    
    SD_optimizer::SD_optimizer(/* args */){

    }
    
    SD_optimizer::~SD_optimizer(){}

    Newton_Rahpson_Optimizer::Newton_Rahpson_Optimizer(Function f, double precision){
        this->f = f;
        this->precision = precision;
    }
    
    Newton_Rahpson_Optimizer::~Newton_Rahpson_Optimizer(){}

    void Newton_Rahpson_Optimizer::set_iterations(size_t n){
        max_iter = n;
    }

    Vector Newton_Rahpson_Optimizer::optimize(Vector x){
        
        if(f.grad(x).norm() < precision) 
            return x;
        
        for(size_t i = 0; i  < max_iter; i++){
            Vector d = f.hess(x).inv() * f.grad(x);
            x -= d;
            
            if(d.norm() < precision){ 
                break;
            }
        }

        return x;
    }

} // namespace NLO
