/*
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2019
 */

#include "tools.hpp"
#include "datafile.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
using std::ifstream;
using std::getline;
using std::stringstream;
using std::string;
using std::cerr;
using std::endl;

#define LINEINFO __FILE__ << ":" << __LINE__

DataFile::DataFile(string fname)
{
    x0=0.01;    // Default value if not given in the file
    filename=fname;
    ifstream file(fname.c_str());
    if (!file.is_open())
    {
        cerr <<  "ERROR! Coudn't read BK solution file " << fname << " " << LINEINFO << endl ;
        exit(1);
    }
    int confid=0;
    while(!file.eof() and confid < 4)
    {
        string line;
        getline(file, line);
        if (line.substr(0, 3)=="###")
        {                    
            switch (confid)
            {
                case 0:
                    minr = StrToReal(line.substr(3,line.length()-3));
                    break;
                case 1:
                    r_multiplier = StrToReal(line.substr(3,line.length()-3));
                    break;
                case 2:
                    rpoints = StrToInt(line.substr(3,line.length()-3));
                    break;
                case 3:
                    x0 = StrToReal(line.substr(3,line.length()-3));
                    if (x0 < 1e-4 or x0>1)
                    {
                        cerr <<"Invalid x0! " << x0 <<", correct x0 must be set by SetX0() before usage! " << LINEINFO << endl;
                        x0=0;
                    }
                    break;
                default:
                    cerr << "File " << fname << " is formatted incorrectly!" << endl;
                    break;
            }
            confid++; 
        }
    }

    // Ok, configurations are read, then read all yvals
    double y=-1;
    std::vector<double> tmpvec;
    while (!file.eof())
    {
        string line;
        getline(file, line);

        // New rapidity?
        if (line.substr(0,3)=="###")
        {
            if (tmpvec.size()>0)
                data.push_back(tmpvec);

            if (tmpvec.size()>0 and tmpvec.size() != rpoints)
            {
                cerr << "File " << fname << ": read " << tmpvec.size() << " points, but "
                << "there should have been " << rpoints << " points, y=" 
                << y << ". " << LINEINFO << endl;
            }

            y = StrToReal(line.substr(3,line.length()-3));
            yvals.push_back(y);
            tmpvec.clear();
            continue;   // Next line is probably new amplitude value
        }
        else if (line.substr(0,1)=="#")
            continue;   // Comment

        // Ok, so this a new amplitude value
        tmpvec.push_back(StrToReal(line));
    }

    // Add last entry
    data.push_back(tmpvec);
    
    file.close();

    if (data[0].size() != rpoints)
    {
        cerr << "File " << fname << ": read " << data.size() << " rpoints, but "
        << "there should have been " << rpoints << " points???" << LINEINFO << endl;
    }
}

void DataFile::GetData(std::vector< std::vector<double> > &n,
                        std::vector<double> &rapidities)
{
    n.clear();
    rapidities.clear();
    // Return vector where indexes are vec[y][r] containing amplitude

    for (uint yind=0; yind < data.size(); yind++)
    {
        std::vector<double> tmpvec;
        for (uint rind=0; rind<rpoints; rind++)
        {
            tmpvec.push_back(data[yind][rind]);
        }
        n.push_back(tmpvec);

        rapidities.push_back(yvals[yind]);
        
    }
    
}

double DataFile::MinR()
{
    return minr;
}

double DataFile::RMultiplier()
{
    return r_multiplier;
}

int DataFile::RPoints()
{
        return rpoints;
}

double DataFile::MaxY()
{
    return yvals[yvals.size()-1];
}

double DataFile::X0()
{
    return x0;
}
