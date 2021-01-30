#pragma once

#include "la.h"

using namespace LA;

namespace RA {
    
    typedef double (*realfunc)(Vector);

    class Function{
    private:
        size_t order;
        realfunc f;
    public:
        Function(){}
        Function(size_t);
        ~Function();

        double operator() (Vector) const;
        void set_func(realfunc);

        double partial(Vector, size_t);
        Vector grad(Vector);
        Matrix hess(Vector);
    };
    
    
    

} // namespace RA
