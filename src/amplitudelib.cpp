/*
 * BK evolved dipole amplitude
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2019
 */

#include "amplitudelib.hpp"
#include "datafile.hpp"
#include "tools.hpp"
#include <string>
#include <cmath>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_min.h>
#include <string>
#include <sstream>


#include <algorithm>

using std::isinf;
using std::isnan;

inline double SQR(double x) { return x*x; }
#define LINEINFO __FILE__ << ":" << __LINE__

/*
 * Load data from a given file
 * Format is specified in file bk/README
 */
AmplitudeLib::AmplitudeLib(std::string datafile)
{
    // Read BK solution
    datafilename = datafile;
    DataFile data(datafile);
    data.GetData(n, yvals);

    // Initialize
    out_of_range_errors=true;
    interpolation_method = SPLINE_LINEAR;
    
    minr = data.MinR();
    rmultiplier = data.RMultiplier();
    rpoints = data.RPoints();
    x0 = data.X0();

    tmprarray = new double[data.RPoints()];
    tmpnarray = new double[data.RPoints()];
    for (int i=0; i<rpoints; i++)
    {
        double tmpr = std::log(minr*std::pow(rmultiplier, i));
        lnrvals.push_back(tmpr);
        rvals.push_back(std::exp(lnrvals[i]));
        tmprarray[i] = std::exp(lnrvals[i]);
    }


    interpolator_xbj=-1;  // if >=0, interpolator is initialized, must free
    // memory (delete tmprarray and tmpnarray at the end)


	std::stringstream ss;
    ss << "Data read from file " << datafile << ", minr: " << minr
        << " maxr: " << MaxR() << " rpoints: " << rpoints << " maxy "
        << yvals[yvals.size()-1] << " x0 " << X0()
        << " Q_{s,0}^2 = " << SaturationScale(x0, 0.393469) << " GeV^2 [ N(r^2=2/Q_s^2, x=x0) = 0.3934]"
        << " (AmplitudeLib v. " << AMPLITUDELIB_VERSION << ")" ;
    info_string = ss.str();
    
    
}

/*
 * Constructor, not loading datafile but data is given directly as an array
 *
 */
AmplitudeLib::AmplitudeLib(std::vector< std::vector< double > > data, std::vector<double> yvals_, std::vector<double> rvals_)
{
    interpolation_method = SPLINE_LINEAR;
    out_of_range_errors = true;
    minr = rvals_[0];
    rmultiplier = rvals_[1] / rvals_[0];
    rpoints = rvals_.size();
    
    rvals = rvals_;
    n = data;
    yvals = yvals_;
    
    x0 = 0.01;
    
    tmprarray = new double[rpoints];
    tmpnarray = new double[rpoints];
    for (int i=0; i<rpoints; i++)
    {
        double tmpr = std::log(minr*std::pow(rmultiplier, i));
        double tmpr2 = rvals[i];
        if ( abs(std::exp(tmpr)/tmpr2-1.0) > 0.01)
        {
            cerr << "WARNING: Dipole sizes does not form a logarithmic grid! tmpr " << tmpr << " r2 " << rvals[i] << " AmplitudeLib::AmplitudeLib()" << endl;
        }
        lnrvals.push_back(tmpr);
        tmprarray[i] = std::exp(lnrvals[i]);
    }
    
   
    interpolator_xbj=-1;
    
    
    std::stringstream ss;
    ss << "#AmplitudeLib initialized , minr: " << minr
    << " maxr: " << MaxR() << " rpoints: " << rpoints << " maxy "
    << yvals[yvals.size()-1] << " x0 " << X0()
    << " Q_{s,0}^2 = " << SaturationScale(x0, 0.393469) << " GeV^2 [ N(r^2=2/Q_s^2, x=x0) = 0.3934]"
    << " (AmplitudeLib v. " << AMPLITUDELIB_VERSION << ")" ;
    info_string = ss.str();
    
    cout << ss.str() << endl;
    
    
}


/*
 * Release reserved memory
 */

AmplitudeLib::~AmplitudeLib()
{
    if (interpolator_xbj>=0)
    {
        delete interpolator;
    }
    delete[] tmpnarray;
    delete[] tmprarray;
}

/*
 * Calculate amplitude interpolating data
 * Interpolate rapidity linearly and r using spline
 */
double AmplitudeLib::N(double r, double xbj)
{
    if (isnan(r) or isinf(r))
    {
        cerr << "r=" << r << " at AmplitudeLib::N(): " << LINEINFO << endl;
        exit(1);
    }
    
    if (r < MinR() or r > MaxR() )
    {
        if (out_of_range_errors)
            cerr << "r must be between limits [" << MinR() << ", " << MaxR() << "]"
                << " asked r=" << r << " " << LINEINFO
                << endl;
        if (r<MinR()) r=MinR()*1.000001; else if (r>MaxR()) r=MaxR()*0.999999;
    }

    double y = std::log(X0()/xbj);
    
    // A slighlty ugly hack, to avoid possible problems around max y
    if ( std::abs(1.0-y / yvals[yvals.size()-1]) < 0.0000001)
		y=yvals[yvals.size()-1];

    
    if (y<0 or y>yvals[yvals.size()-1] )
    {
        //if (out_of_range_errors)
            cerr << "y must be between limits [" << 0 << ", "
                << yvals[yvals.size()-1] << "], asked y=" << y << ", x=" << xbj << " (x0=" << X0() << ") datafile " << datafilename << " "
                << LINEINFO << endl;
            exit(1);
        //if (y < 0) y=0; else if (y>yvals[yvals.size()-1]) y=yvals[yvals.size()-1];
    }

    
    // Use already initialized interpolator
    if (InterpolatorInitialized(xbj))
    {
        double result=0;

        // Can't interpolate (too large dipole), return 1.0
        if (r >= maxr_interpolate and maxr_interpolate>0) { return 1.0; }

        result = interpolator->Evaluate(r);
        
        if (result>1) return 1;                // limit N(r)<=1
        if (result<0) return 0;              // limit N >= 0
        return result;

    }
    
    int rind = FindIndex(r, rvals);
    int yind = FindIndex(y, yvals);
    
    // Check interpolation method, use bspline if needed
    if (interpolation_method == LINEAR_LINEAR)
    {
        int rind2 = rind+1;
        int yind2 = yind+1;
        if (rind2 >= rvals.size()) return 1.0;
        if (yind2 >= yvals.size()) return n[yind][rind]; // this should never happen
        double r1_ = rvals[rind];
        double r2_ = rvals[rind2];
        double y1_ = yvals[yind];
        double y2_ = yvals[yind2];
        double bilin =  (1.0/( (r2_ - r1_)*(y2_ - y1_) ))*( (n[yind][rind])*(r2_ - r)*(y2_ - y) + (n[yind][rind2])*(r - r1_)*(y2_ - y) + (n[yind2][rind])*(r2_ - r)*(y - y1_) + (n[yind2][rind2])*(r - r1_)*(y - y1_) );
        
        return bilin;
        
    }

    /// Initialize new interpolator and use it
    int interpolation_points = INTERPOLATION_POINTS;

    int interpolation_start, interpolation_end;
    if (rind - interpolation_points/2 < 0)
    {
        interpolation_start=0;
        interpolation_end=interpolation_points;
    }
    else if (rind + interpolation_points/2 > RPoints()-1 )
    {
        interpolation_end = RPoints()-1;
        interpolation_start = RPoints()-interpolation_points-3;
    }
    else
    {
        interpolation_start = rind - interpolation_points/2;
        interpolation_end = rind + interpolation_points/2;
    }

    int interpo_points = interpolation_end - interpolation_start+1;
    double *tmparray = new double[interpo_points];
    double *tmpxarray = new double[interpo_points];
    for (int i=interpolation_start; i<= interpolation_end; i++)
    {
        tmpxarray[i-interpolation_start]=rvals[i];

        tmparray[i-interpolation_start] = n[yind][i];

        // Interpolate in y if possible
        if (yind < static_cast<int>(yvals.size()-1) )
        {
            tmparray[i-interpolation_start]=n[yind][i] 
                + (y - yvals[yind]) * (n[yind+1][i] - n[yind][i])
                / (yvals[yind+1]-yvals[yind]);
        }
    }
    
    Interpolator interp(tmpxarray, tmparray, interpo_points);
    interp.Initialize();
    double result=0;

    result = interp.Evaluate(r);

    
    delete[] tmparray;
    delete[] tmpxarray;

    return result;

}

bool AmplitudeLib::InterpolatorInitialized(double xbj)
{
    if (interpolator_xbj > 0 and std::abs(xbj - interpolator_xbj)/std::min(xbj, interpolator_xbj) < 0.001)
        return true;
    else
        return false;
}

double AmplitudeLib::S(double r, double xbj)
{
    double s = 1.0 - N(r,xbj);
    if (s<=0) return 0;
    if (s>=1.0) return 1.0;
    return s;
}

// Dipole amplitude as a functoin of rapidity, wrapper y -> xbj
double AmplitudeLib::Ny(double r, double y)
{
    return N(r, x0*std::exp(-y));
}

double AmplitudeLib::Sy(double r, double y)
{
    return S(r, x0*std::exp(-y));
}




/*
 * Calculate saturation scale, definition is
 * N(r_s, y) = Ns = (for example) 0.5
 */
struct SatscaleSolverHelper
{
    double xbj;
    double Ns;
    double gammac;
    AmplitudeLib* N;
};
double SatscaleSolverHelperf(double r, void* p)
{
    SatscaleSolverHelper* par = (SatscaleSolverHelper*)p;
    return par->N->N(r, par->xbj) - par->Ns;
}


double AmplitudeLib::SaturationScale(double xbj, double Ns)
{
    
    SatscaleSolverHelper helper;
    helper.xbj=xbj; helper.Ns=Ns; helper.N=this; helper.gammac=0.5;
    const int MAX_ITER = 1000;
    const double ROOTFINDACCURACY = 0.00001;
    gsl_function f;
    f.params = &helper;
        
    f.function = &SatscaleSolverHelperf;

    const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
        
    gsl_root_fsolver_set(s, &f, MinR()*1.0001, MaxR()*0.999);
    int iter=0; int status; double min,max;
    do
    {
        iter++;
        gsl_root_fsolver_iterate(s);
        min = gsl_root_fsolver_x_lower(s);
        max = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(min, max, 0, ROOTFINDACCURACY);    
    } while (status == GSL_CONTINUE and iter < MAX_ITER);

    if (iter>=MAX_ITER)
        cerr << "Solving failed at x=" << xbj << " " << LINEINFO << endl;


    double sat = gsl_root_fsolver_root(s);

    gsl_root_fsolver_free(s);

    return 2.0/(sat*sat);
    
}

/*
 * Amplitude in adjoint representation
 * Coordinate space: N_A(r) = 2N(r)-N(r)^2
 */
double AmplitudeLib::N_A(double r, double y)
{

    double result=0;

    double n = N(r,y);
    result = 2.0*n - n*n;
    if (result>1) result=1;
    if (result<0) result=0;
    return result;
    
    

    return result;
}




/*
 * Initializes interpolation method with all read data points at given y
 */
void AmplitudeLib::InitializeInterpolation(double xbj)
{
	if (xbj > X0())
	{
		cerr << "Asked to initialize interpolator with too large x=" << xbj
        << " (x0=" << X0() <<")! "
        << "Dont know what to do, panicking... " << LINEINFO << endl;
		exit(1);
	}

    double y = std::log(xbj/X0());
	if (y>MaxY())
	{
		cerr << "Asked to initialize interpolator with too large y=" << y <<", maxy=" << MaxY() << endl;
		exit(1);
	}
	
    if (InterpolatorInitialized(xbj)) return;    // Already done it
    
    if (interpolator_xbj>=0)
    {
        interpolator_xbj=-1;
        delete interpolator;
    }
    for (int i=0; i<rpoints; i++)
    {
        double tmpr = tmprarray[i];
        if (i==0) tmpr*=1.0001; if (i==rpoints-1) tmpr*=0.9999;
        tmpnarray[i] = N(tmpr, xbj);
    }
    interpolator = new Interpolator(tmprarray, tmpnarray, rpoints);

    interpolator->Initialize();
    interpolator_xbj = xbj;
    
    maxr_interpolate=-1;
    int iter=0;
    // Find r s.t. N(r)==1
    double step=2; double prevr=0.1;
    const int MAXITER=40;
    for (double r=0.01; r<=MaxR(); r+=step)
    {
        if (N(r,xbj)>=0.999999)		// check that this accuracy is the same as in rbk/src/solver.cpp
        {
            if (step<1e-2)
            {
                maxr_interpolate = r;
                break;
            }
            step /= 1.5;
            r=prevr;
        }
        prevr=r;
        iter++;
        if (iter > MAXITER)
        {
            //cerr << "Didn't find maxr_interpolate at y=" << y <<", ignoring it " << LINEINFO << endl;
            maxr_interpolate=-1;
            break;  // Didn't find, dont force any upper limit
        }
    }
}



/*
 * Initializes interpolator and returns it
 */
Interpolator* AmplitudeLib::MakeInterpolator(double xbj)
{
    std::vector<double> tmpnvals;
    for (int i=0; i<rpoints; i++)
    {
        tmpnvals.push_back(N(rvals[i], xbj));
    }
    Interpolator* inter = new Interpolator(rvals, tmpnvals);
    inter->Initialize();

    return inter;
        
}


int AmplitudeLib::RPoints()
{
    return rpoints;
}
double AmplitudeLib::MinR()
{
    return minr;
}

double AmplitudeLib::MaxR()
{
    return minr*std::pow(rmultiplier, rpoints-1);
}

int AmplitudeLib::YPoints()
{
    return yvals.size();
}

double AmplitudeLib::MaxY()
{
    return yvals[yvals.size()-1];
}

bool AmplitudeLib::SetOutOfRangeErrors(bool er)
{
    bool outofrange = out_of_range_errors;
    out_of_range_errors=er;
    return outofrange;    
}

double AmplitudeLib::X0()
{
    return x0;
}

void AmplitudeLib::SetX0(double x0_)
{
	x0=x0_;
}

std::string AmplitudeLib::GetString()
{
	return info_string;
}

std::string AmplitudeLib::Version()
{
	std::stringstream s;
	s << AMPLITUDELIB_VERSION << " (build " <<  __DATE__ << " " << __TIME__ << ")";
	return s.str();
}


/*
 * Return ith rapidity value
 */
double AmplitudeLib::YValue(int yind)
{
    if (yind >= YPoints() or yind<0)
    {
        cerr << "Asked rapidity of index " << yind <<", but the maximum index is " << YPoints()-1 << " " << LINEINFO << endl;
        exit(1);
    }
    return yvals[yind];
}
