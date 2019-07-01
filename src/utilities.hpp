#ifndef D__UTILITIES__HPP___
#define D__UTILITIES__HPP___
#include <ostream>
#include <iostream>
#include <map>
#include <boost/foreach.hpp>
#include <boost/numeric/interval.hpp>
#include <boost/assert.hpp>
#include <vector>
#include <algorithm>
#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>


using namespace std;
using namespace boost;

#ifdef USE_BOOST_INTERVAL

using namespace boost::numeric;
using namespace interval_lib;

// http://stackoverflow.com/questions/15860761/c-boost-interval-and-cos                                                                                                                                   
typedef interval<double, policies<save_state<rounded_transc_std<double> >,
                                  checking_base<double> > > i_double;
#else 

#include "mpfiWrapper.hpp"
typedef MpfiWrapper i_double;

#endif



typedef adjacency_list <vecS, vecS, undirectedS> Graph;

// http://stackoverflow.com/questions/15860761/c-boost-interval-and-cos


extern bool debug;
ostream & operator<< (ostream & os, i_double const & I);

bool setsIntersect(set<int> const &a, set<int> const & b);

/*-- Graph related utilities --*/

int findMaximumDegree(Graph const & gRet);
void outputGraphToDot(Graph const & g, ostream & os);
int findMaximumDegree(Graph const & g);
int  getComponentsOfGraph(Graph const & g, std::vector<int> & componentIDs);
void applySubGaussianCMIInequalitiesOverRange(ostream & os, i_double range, i_double expect, double cSquared, int nSubdivs);
void applyBernsteinCMIInequalityOverRange(ostream & os, i_double range, i_double expect, double var, double M, int nSubDivs);
#define MAX(a,b) ((a) >= (b)?(a):(b))

#endif
