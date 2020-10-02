/*
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2019
 */

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include "interpolation.hpp"


#ifndef LINEINFO
    #define LINEINFO __FILE__ << ":" << __LINE__
#endif

typedef unsigned int uint;
using std::isinf;
using std::isnan;

/*
 * Intialize interpolation
 * Returns -1 in case of error, 0 otherwise
 */
int Interpolator::Initialize()
{
    int status=0;
    out_of_range_errors = true;
    if (ready)
    {
        // Interpolator is already initialized: clear it (GSL part)
        // and re-initialize
        switch(method)
        {
            case INTERPOLATE_SPLINE:
                if (spline != NULL)
                {
                    gsl_spline_free(spline);
                    spline=NULL;
                }
                if (acc != NULL)
                {
                    gsl_interp_accel_free(acc);
                    acc=NULL;
                }
                break;
            default:
                cerr << "Undefined interpolation method!" << endl;
                exit(1);
        }
        ready=false;
        
    }
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            acc = gsl_interp_accel_alloc();
            spline = gsl_spline_alloc(gsl_interp_cspline, points);
            status = gsl_spline_init(spline, xdata, ydata, points);
            break;
 
        default:
            exit(1);
    }
    ready=true;
    if (status)
    {
        cerr << "Interpolator initialization failed at " << LINEINFO << endl;
        return -1;
    }
    return 0;   //ok, there is no error handling at the moment...
}


double Interpolator::Evaluate(double x)
{
    if (isnan(x) or isinf(x))
    {
        cerr << "Trying to evaluate interpolator with x=" << x << " at " << LINEINFO << endl;
        exit(1);
    }
    
    if (!ready)
    {
        cerr << "Interpolator is not ready! Did you forget to call Interpolator::Initialize()?" << endl;
        return 0;
    }

    if (x<minx or x>maxx)
    {
		if (freeze)
		{
			if (x<minx) return freeze_underflow;
			else return freeze_overflow;
		}
		if (x < 0.9999*minx or x > 1.00001*maxx)	// if not true, no need to display error
        {
            if (out_of_range_errors)
                cerr << "x=" << x << " is not within limits [" << minx << ", " << maxx << "], forcing "
                    << "it in that interval! " << LINEINFO << endl;
        }
        if (x<minx) x=minx*1.00001;
        if (x>maxx) x=maxx*0.999999;
    }
    
    double res, yerr; int status;
    res=0;
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            status = gsl_spline_eval_e(spline, x, acc, &res);
            if (status)
            {
                cerr << "Interpolation failed at " << LINEINFO << ", error " << gsl_strerror(status)
                 << " (" << status << "), x=" << x << ", minx=" << xdata[0]
                 << ", maxx=" << xdata[points-1] << ", result=" << res << endl;
                 exit(1);
            }
            break;
        default:
            cerr << "Interpolation method is invalid! " << LINEINFO << endl;
            exit(1);
    }

    if (isnan(res) or isinf(res))
    {
        cerr << "Interpolation at x=" << x << " gives " << res << endl;
		return 0;
        exit(1);
    }

    
    return res;   
}

double Interpolator::Derivative(double x)
{
    double res=0; int status=0;
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            status = gsl_spline_eval_deriv_e(spline, x, acc, &res);
            break;
    }
    if (status)
        cerr << "An error occurred while evaluating the derivative at x=" << x
        << " result " << res << " " << LINEINFO << endl;

    return res;
}

double Interpolator::Derivative2(double x)
{
    double res; int status=0;
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            status = gsl_spline_eval_deriv2_e(spline, x, acc, &res);
            break;
    }

    if (status)
    {
        cerr << "2nd derivative interpolation failed at x=" << x <<", result "
        << res << " " << LINEINFO << endl;
    }
    return res;

}

Interpolator::Interpolator(double *x, double *y, int p)
{
    points=p;
    xdata=x;
    ydata=y;
    minx=x[0];
    maxx=x[p-1];
    method = INTERPOLATE_SPLINE;
    allocated_data=false;
    ready=false;
    freeze=false;
    freeze_underflow = y[0];
    freeze_overflow = y[p-1];

    for (int i=0; i<p; i++)
    {
        // Check that x values are monotonically increasing
        if (i>0)
        {
            if (xdata[i-1]>=xdata[i])
            {
                cerr << "Grid points are not monotonically increasing! grid["
                    << i-1 <<"]=" << xdata[i-1] <<", grid["<<i<<"]="<< xdata[i]
                    << " " << LINEINFO << endl;
            }
        }
    }

    Initialize();
}

Interpolator::Interpolator(std::vector<double> &x, std::vector<double> &y)
{
    points = x.size();
    xdata = new double[points];
    ydata = new double[points];
    allocated_data=true;

    for (uint i=0; i<x.size(); i++)
    {
        xdata[i]=x[i];
        ydata[i]=y[i];

        // Check that x values are monotonically increasing
        if (i>0)
        {
            if (xdata[i-1]>=xdata[i])
            {
                cerr << "Grid points are not monotonically increasing! grid["
                    << i-1 <<"]=" << xdata[i-1] <<", grid["<<i<<"]="<< xdata[i]
                    << " " << LINEINFO << endl;
                exit(1);
            }
        }
    }
    minx=xdata[0]; maxx=xdata[x.size()-1];
    method = INTERPOLATE_SPLINE;
    ready=false;
    freeze=false;
    freeze_overflow = y[y.size()-1];
    freeze_underflow = y[0];

    Initialize();
}

void Interpolator::SetMethod(INTERPOLATION_METHOD m)
{
    method = m;
}

void Interpolator::Clear()
{
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            if (spline != NULL)
            {
                gsl_spline_free(spline);
                spline=NULL;
            }
            if (acc != NULL)
            {
                gsl_interp_accel_free(acc);
                acc=NULL;
            }
            break;
    }

    if (allocated_data)
    {
        delete[] xdata;
        delete[] ydata;
        allocated_data=false;
    }
}

Interpolator::~Interpolator()
{
    Clear();

}

double* Interpolator::GetXData()
{
    return xdata;
}
double* Interpolator::GetYData()
{
    return ydata;
}
int Interpolator::GetNumOfPoints() const
{
    return points;
}
INTERPOLATION_METHOD Interpolator::GetMethod() const
{
    return method;
}

// Copy data from given class and initialize this, as this is
// the copy constructor
Interpolator::Interpolator(const Interpolator& inter)
{
    points=inter.GetNumOfPoints();
    xdata = new double[points];
    ydata = new double[points];
    allocated_data=true;


    gsl_spline *tmpspline = inter.GetGslSpline();
    for (int i=0; i<points; i++)
    {
        xdata[i] = tmpspline->x[i];
        ydata[i] = tmpspline->y[i];
    }
    minx = xdata[0]; maxx=xdata[points-1];
    method = inter.GetMethod();
    ready=false;
    Initialize();
}

gsl_spline* Interpolator::GetGslSpline() const
{
    return spline;
}

double Interpolator::MinX()
{
	return minx;
}

double Interpolator::MaxX()
{
	return maxx;
}


bool Interpolator::Freeze()
{
	return freeze;
}
void Interpolator::SetFreeze(bool f)
{
	freeze=f;
}
void Interpolator::SetUnderflow(double min)
{
	freeze_underflow=min;
}
 void Interpolator::SetOverflow(double max)
 {
	 freeze_overflow=max;
 }
double Interpolator::UnderFlow()
{
	 return freeze_underflow;
}
double Interpolator::OverFlow()
{
	return freeze_overflow;
}

void Interpolator::SetMaxX(double x)
{
    maxx=x;
}

void Interpolator::SetMinX(double x)
{
    minx=x;
}

void Interpolator::SetOutOfRangeErrors(bool er)
{
    out_of_range_errors=er;
}
