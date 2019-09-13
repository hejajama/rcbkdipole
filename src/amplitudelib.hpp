/*
 * BK evolved dipole amplitude
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2019
 */

#ifndef _AMPLITUDELIB_H
#define _AMPLITUDELIB_H

#include "interpolation.hpp"
#include <vector>
#include <string>


enum AMPLITUDE_INTERPOLATION_METHOD
{
    SPLINE_LINEAR,      // spline in r, linear in y
    LINEAR_LINEAR       // linear in r and y
                // thread safe, no need to initialize interpolation
};



/**
 * Amplitude class
 *
 * Stores the loaded dipole amplitude (loaded by DataFile class), and
 * evaluates the amplitude at given r,xbj.
 * Also computes Fourier transfers between transverse coordinate and
 * momentum spaces.
 */
class AmplitudeLib
{
    public:
        /**
         * Loader.
         * 
         * Constructor, loads the solution for the BK equation from
         * the given file and stores the result in memory.
         */
        AmplitudeLib(std::string datafile);
    
        /* Load amplitude from the give data tables
         *
         * data[i][j] is dipole amplitude at rapidity yvals[i] for dipole size rvals[i]
         */
        AmplitudeLib(std::vector< std::vector< double > > data, std::vector<double> yvals_, std::vector<double> rvals_);
    
        ~AmplitudeLib();

        /**
         * Dipole amplitude
         *
         * Evaluates the dipole amplitude at given dipole size r and at
         * given Bjorken x (xbj). In momentum space the first argument is
         * transverse momentum.
         *
         * @param r dipole size in GeV^-1
         * @param xbj Bjorken x at which dipole is evaluated
         */
        double N(double r, double xbj);
    
        /**
         * Dipole amplitude at given rapidity y
         *
         * Evaluates the dipole amplitude at given dipole size r and at
         * given evolution rapidity y (y=0 corresponding to the initial condition).
         * In momentum space the first argument is
         * transverse momentum.
         *
         * @param r dipole size in GeV^-1
         * @param y Evolution rapidity at which the dipole is evaluated
         */
        double Ny(double r, double y);

        /**
         * Dipole amplitude in adjoint representation
         *
         * Evaluates the dipole amplitude in adjoint representation,
         * 2N(r)-N(r)^2, at given dipole size r and at
         * given Bjorken x (xbj). In momentum space the first argument is
         * transverse momentum.
         */
        double N_A(double r, double xbj);

        /**
         * Scattering matrix
         *
         * Evaluate the scattering matrix 1-N (forced between 0 and 1)
         * at given dipole size and Bjorken-x.
         *
         * @param r dipole size in GeV^-1
         * @param xbj Bjorken x at which dipole is evaluated
         */
         double S(double r, double xbj);
    
        /**
         * Scattering matrix at given evolution rapidity
         *
         * Evaluate the scattering matrix 1-N (forced between 0 and 1)
         * at given dipole size and evolution rapidity,
         * (y=0 is the initial condition).
         *
         * @param r dipole size in GeV^-1
         * @param y Evolution rapidity at which the dipole is evaluated
         */
        double Sy(double r, double y);
        
    
        /**
         * Initialize interpolation for certain Bjorken-x
         */
        void InitializeInterpolation(double xbj);

        /**
         * Test if interpolator is initialized at given Bjorken-x
         */
        bool InterpolatorInitialized(double xbj);

        /**
         * Create interpolator at given Bjorken-x and return it
         */
        Interpolator* MakeInterpolator(double xbj);

    
        /**
         * Saturation scale
         *
         * Solve saturation scale defined as N(r^2=2/Q_s^2) = N_s
         * 
         * @return Saturation scale in GeV^2
         */
        double SaturationScale(double xbj, double Ns);

        
        /**
         * Number or rapidity values in the BK solution
         */
        int YPoints();

        /**
         * Return ith rapidity value
         *
         * Can be used to avoid interpolation uncertainties 
         */
        double YValue(int yind);
         

        /**
         * Number of dipole sizes in the BK solution
         */
        int RPoints();

        /**
         * Minimum dipole size
         */
        double MinR();

        /**
         * Maximum dipole size
         */
        double MaxR();
        double MaxY();
        
        /**
         * Initial Bjorken-x in the BK solution
         *
         * Can we owerwritten using SetX0() method
         * @see SetX0
         */
        double X0();

        /**
         * Set initial Bjorken-x
         *
         * Override the Bjoken-x at initial condition for the BK solution
         * (x0).
         */
        void SetX0(double x0_);

        /**
         * Specify whether an error message to stderr is printed for too
         * small/alrge dipoles or small/large xbj values
         *
         * @return previous setting
         */

        bool SetOutOfRangeErrors(bool er);

        /**
         * Returns an info string describing the BK solution and AmplitudeLib version
         */
        std::string GetString();

    

        /**
         * Return version string
         */
        std::string Version();

    
        /**
         * Set interpolation method
         */
        void SetInterpolationMethod(AMPLITUDE_INTERPOLATION_METHOD m) { interpolation_method = m;}

    private:
        // [yind][r/kind]
        std::vector< std::vector<double> > n;
        std::vector<double> yvals;
        std::vector<double> lnrvals;
        std::vector<double> rvals;
        Interpolator *interpolator;

        double interpolator_xbj;  //! xbj at which the interpolator is initialized
        double* tmprarray;
        double* tmpnarray;

        double minr;
        double rmultiplier;

        /**
         * Max dipole size for interpolation
         *
         * When interpolating in coordinate space,
         * force N(r>maxr_interpolate)=1 (avoid osciallatory 
         * artifacts from interpolation code)
         */
        double maxr_interpolate;
    
        /**
         * Method to use to interpolate dipole amplitude
         *
         * If set to LINEAR_LINEAR, one does not use splines,
         * and the code is completely thread safe even if
         * amplitude is evaluated simultaenously at different
         * rapidities
         */
        AMPLITUDE_INTERPOLATION_METHOD interpolation_method;
        
        int rpoints;    //! Number of dipole sizes

        double x0;      //! Initial Bjorken-x
        
        

        bool out_of_range_errors;  //! If true, don't print "out of range" errors
        
        std::string info_string;
 
        std::string datafilename;   //! Name of the file where the amplitude is read, for error messages

        
};




const int INTERPOLATION_POINTS = 12;




const std::string AMPLITUDELIB_VERSION = "2019-09-public";

#endif
