#pragma once

#include "ra.h"

using namespace RA;
using namespace LA;

namespace NLO {
    
    class SD_optimizer {
    private:
        
    public:
        SD_optimizer(/* args */);
        ~SD_optimizer();
    };

    class Newton_Rahpson_Optimizer{
    private:
        Function f;
        double precision;
        size_t max_iter;
    public:
        Newton_Rahpson_Optimizer(Function, double prec=1e-6);
        ~Newton_Rahpson_Optimizer();

        void set_iterations(size_t);
        Vector optimize(Vector);
    };

    class CG_Optimizer{
    private:
        Function f;
        double precision;
        size_t max_iter;
    public:
        CG_Optimizer(Function, double prec=1e-6);
        ~CG_Optimizer();

        void set_iterations(size_t);
        Vector optimize(Vector);
    };
    
    

} // namespace NLO
