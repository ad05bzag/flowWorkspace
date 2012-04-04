/*
 * GatingSet.hpp
 *
 *  Created on: Mar 15, 2012
 *      Author: wjiang2
 */



#ifndef GATINGSET_HPP_
#define GATINGSET_HPP_
#include "GatingHierarchy.hpp"
#include "flowJoWorkspace.hpp"
#include <string>
#include <map>
#include <boost/foreach.hpp>

using namespace std;

/*GatingSet is multiple GatingHierarchies that has the flow data associated and gated*/
class GatingSet{
	workspace * ws;
	map<string,GatingHierarchy *> ghs;
	/*
	 * we assume the GatingSet is mainly for read-access
	 * so we have a copy of samplenames here as a vector
	 * for the convienient access samplenames
	 * If there is need to avoid data-inconsistency,we can loop through the map
	 * and return the sampleList
	 */
	vector<string> sampleList;

public:
	~GatingSet();
	GatingSet(string);
	GatingHierarchy * getGatingHierarchy(string );
	GatingHierarchy * getGatingHierarchy(unsigned int);
	void parseWorkspace(unsigned short);
	vector<string> getSamples(){this->sampleList;};
};
#endif /* GATINGSET_HPP_ */