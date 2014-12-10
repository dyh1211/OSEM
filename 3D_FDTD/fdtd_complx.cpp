#include <iostream>
#include <cmath>
#include "fdtd_complx.hpp"
using namespace std;


double real(complx in)
{
    return (in.real);
}


double imag(complx in)
{
    return (in.imag);
}


double abs(complx in)
{
    return sqrt(in.real*in.real + in.imag*in.imag);
}

double norm(complx in)
{
    return pow(sqrt(in.real*in.real + in.imag*in.imag),2);
}

complx conj(complx in)
{
    in.imag = in.imag * -1;
    return (in);
}

// define constructor
complx::complx( double r, double i )
{
    real = r;
    imag = i;
}

// define overloaded += (plus) operator
void complx::operator+= (const complx& c)
{
    real = (this->real + c.real);
    imag = (this->imag + c.imag);
}

// define overloaded -= (plus) operator
void complx::operator-= (const complx& c)
{
    real = (this->real - c.real);
    imag = (this->imag - c.imag);
}


// define overloaded *= (plus) operator
void complx::operator*= (const complx& c)
{
    real = (this->real * c.real);
    imag = (this->imag * c.imag);
}


// define overloaded + (plus) operator
complx complx::operator+ (const complx& c) const
{
    complx result;
    result.real = (this->real + c.real);
    result.imag = (this->imag + c.imag);
    return result;
}


// define overloaded + (plus) operator
complx complx::operator+ (const double& c) const
{
    complx result;
    result.real = (this->real + c);
    return result;
}

// define overloaded - (unary minus) operator
complx complx::operator- (void)
{
    complx result;
    result.real = (this->real * -1);
    result.imag = (this->imag * -1);
    return result;
}


// define overloaded - (binary minus) operator
complx complx::operator- (const complx& c) const
{
    complx result;
    result.real = (this->real - c.real);
    result.imag = (this->imag - c.imag);
    return result;
}


// define overloaded / (divide) operator
complx complx::operator/ (const complx& c) const
{
    complx result;
    result.real = (this->real * c.real + this->imag * c.imag) / (c.real * c.real + c.imag * c.imag);
    result.imag = (this->imag * c.real - this->real * c.imag) / (c.real * c.real + c.imag * c.imag);
    return result;
}

// define overloaded * (multiply) operator
complx complx::operator* (const complx& c) const
{
    complx result;
    result.real = (this->real * c.real - this->imag * c.imag);
    result.imag = (this->real * c.imag + this->imag * c.real);
    return result;
}


// define overloaded * (*) operator
complx complx::operator* (const double& c) const
{
    complx result;
    result.real = (this->real * c);
    result.imag = (this->imag * c);
    return result;
}



// define overloaded = (minus) operator
void complx::operator= (const complx& c)
{
    real = c.real;
    imag = c.imag;
}


