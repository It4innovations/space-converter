// zero_check.cpp
// This program searches for "0.0" positions in a tipsy file
// This is to check for potential I/O errors.

#include <iostream>

#include "TipsyReader.h"

#include <assert.h>

using namespace std;
using namespace Tipsy;

template <typename TPos, typename TVel>
void check_all_parts(TipsyReader &r, double dEps) 
{
    header hHeader = r.getHeader();

    int iZcount = 0;  // count of zeros
    for(int i = 0; i < hHeader.nsph; i++) {
        gas_particle_t<TPos, TVel> p;
        r.getNextGasParticle_t(p);
	assert(p.mass >= 0.0);
	assert(isfinite(p.pos.x));
	assert(isfinite(p.pos.y));
	assert(isfinite(p.pos.z));
        if((fabs(p.pos.x) < dEps) && (fabs(p.pos.y) < dEps) && (fabs(p.pos.z) < dEps))
            iZcount++;
    }
    for(int i = 0; i < hHeader.ndark; i++) {
        dark_particle_t<TPos, TVel> dp;
        r.getNextDarkParticle_t(dp);
	assert(dp.mass >= 0.0);
	assert(isfinite(dp.pos.x));
	assert(isfinite(dp.pos.y));
	assert(isfinite(dp.pos.z));
        if((fabs(dp.pos.x) < dEps) && (fabs(dp.pos.y) < dEps) && (fabs(dp.pos.z) < dEps))
            iZcount++;
    }
    for(int i = 0; i < hHeader.nstar; i++) {
        star_particle_t<TPos, TVel> p;
        r.getNextStarParticle_t(p);
	assert(p.mass >= 0.0);
	assert(isfinite(p.pos.x));
	assert(isfinite(p.pos.y));
	assert(isfinite(p.pos.z));
        if((fabs(p.pos.x) < dEps) && (fabs(p.pos.y) < dEps) && (fabs(p.pos.z) < dEps))
            iZcount++;
    }
    cout << "Total zeros " << iZcount << endl;
}

int main(int argc, char** argv) {
    if(argc < 3) {
        cerr << "Usage: zero_check eps tipsyfile" << endl;
        return 1;
    }
    TipsyReader r(argv[2]);
    double dEps = atof(argv[1]);
    
    if(!r.status()) {
        cerr << "Couldn't load the tipsy file properly!" << endl;
        return 3;
    }
    if(r.isDoublePos() && r.isDoubleVel())
        check_all_parts<double,double>(r, dEps);
    else if (r.isDoublePos())
        check_all_parts<double,float>(r, dEps);
    // else if (r.isDoubleVel())
    //    check_all_parts<float,double>(r, dEps);
    else
        check_all_parts<float,float>(r, dEps);
    return 0;
}
