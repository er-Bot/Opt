#pragma once

#include "m.h"

#include <cstdlib>
#include <cstring>
#include <string>
#include <exception>

namespace LA {
    
    class Vector{
    private:
        size_t len;
        double* ptr;
    public:
        Vector();
        Vector(size_t);
        Vector(double*, size_t);
        Vector(const Vector&);
        ~Vector();

        void print();
        size_t length();

        double& operator() (size_t);
        double operator() (size_t) const;
        double operator*(Vector);
        Vector operator-(Vector);
        Vector operator-=(Vector);
        Vector operator+(Vector);
        Vector operator+=(Vector);
        Vector operator*(double);
        Vector operator/(double);
        Vector operator+(double);
        Vector operator-(double);
        Vector operator*=(double);
        Vector operator/=(double);
        Vector operator+=(double);
        Vector operator-=(double);
        double dot(Vector);

        double norm(size_t n = 2);
        double sum();
    };

    class Matrix{
    private:
        size_t row;
        size_t col;
        double* ptr;
    public:
        Matrix();
        Matrix(size_t);
        Matrix(size_t, size_t);
        Matrix(const Matrix&);
        Matrix& operator=(const Matrix&);
        Matrix(const Vector v);
        ~Matrix();

        void random(double min=0, double max=1, bool sym=false);
        void from_list(double*, size_t);
        void print();

        double get_elem(size_t, size_t);
        void set_elem(size_t, size_t, double);
        double& operator() (size_t, size_t);
        double operator() (size_t, size_t) const;
        Vector get_col(size_t);
        Vector get_row(size_t);
        Vector get_diag();
        void set_col(size_t, Vector);
        void set_row(size_t, Vector);
        size_t rows();
        size_t cols();

        static Matrix identity(size_t);
        static Matrix ones(size_t, size_t);
        static Matrix diag(double*, size_t);
    
        size_t rank();
        Matrix get_minor(size_t, size_t);
        double det();
        Matrix inv();
        double norm();
        double norm_infty();
        double norm_frobinus();

        Matrix T();
        Matrix operator+(Matrix);
        Matrix operator-(Matrix);
        Matrix operator*(Matrix);
        Vector operator*(Vector);
        Matrix operator/(Matrix);
        Matrix operator+(double);
        Matrix operator-(double);
        Matrix operator*(double);
        Matrix operator/(double);
        Matrix operator+=(Matrix);
        Matrix operator-=(Matrix);
        Matrix operator+=(double);
        Matrix operator-=(double);
        Matrix operator*=(double);
        Matrix operator/=(double);
        Matrix dot(Matrix);

        bool compare(Matrix, double);
        bool is_square();
        bool is_symmetric();
        bool is_diag();
        bool is_null();
        bool is_pd();
    };

    class Exception : public std::exception {
    private:
        std::string _msg;
    public:
        Exception(const std::string& msg) : _msg(msg) {}
        ~Exception(){}
        const char* what() const noexcept {return _msg.c_str();}
        std::string get_message() const { return _msg;}
    };
    class IndexOutOfRange : public Exception {
    public:
        IndexOutOfRange() : Exception("IndexOutOfRange: Index out of range exception!") {}
        IndexOutOfRange(const std::string& msg) : Exception("IndexOutOfRangeIndexOutOfRange: " + msg) {}
    };
    class IncompatibleArguments : public Exception {
    public:
        IncompatibleArguments() : Exception("IncompatibleArguments: Incompatible arguments exception!") {}
        IncompatibleArguments(const std::string& msg) : Exception("IncompatibleArguments: " + msg) {}
    };

    
} // namespace LA