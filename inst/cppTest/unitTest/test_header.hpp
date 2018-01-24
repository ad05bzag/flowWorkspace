/*
 * test_header.hpp
 *
 *  Created on: May 15, 2012
 *      Author: wjiang2
 */

#ifndef TEST_HEADER_HPP_
#define TEST_HEADER_HPP_



#include "flowWorkspace/openWorkspace.hpp"
#include "flowWorkspace/flowJoWorkspace.hpp"
#include "cytolib/GatingSet.hpp"
#include "cytolib/GatingHierarchy.hpp"
#include "cytolib/transformation.hpp"
#include "cytolib/spline.hpp"
using namespace std;
using namespace cytolib;
using namespace cytoml;

struct testCase{
	string filename; //xml file name
	unsigned short wsType; //workspace type
	string colfile; // text file that records the compensated channel names deprecated
	string ncfile; // raw data stored in hdf format deprecated
	string fcs; // raw data stored in hdf format
	map<string,string> samples; // fcs file name vs sampleID
	SAMPLE_NAME_LOCATION sampNloc; // the location where the sample name to be parsed
	string archive; // archived gating set dat file
	vector<bool> isEqual; // the bool vector records the counts discrepancy (using cv) between flowJo and flowCore
	float tolerance; // the threshold for cv value
	bool isParseGate; //whether to parse gate from xml
	int xmlParserOption;//xml parser option passed down to libxml2
	bool isTemplate;// whether test the template copying feature
	bool isLoadArhive;// whether to load archived gs
	bool isSaveArchive;
	unsigned archiveFormat;
	bool archiveType;// boost or google
	map<string,float> gains;
	vector<VertexID> skipPops;
//	vector<double> times;//global variable to collect run time

} ;
typedef MemCytoFrame FRAMETYPE;

void gs_gating(GatingSet<FRAMETYPE> &gs,string curSample, string fcs, map<string,float> &gains);
void gh_counts(GatingHierarchy<FRAMETYPE>& gh,vector<bool> &isEqual, const float tolerance);
void clone_test(testCase myTest);
//void gs_parse(testCase,unsigned short,bool,bool);
void parser_test(testCase &);
void ncdf_test();
void compCalTbl();
void spline_test();
void cpConsTest();
#endif /* TEST_HEADER_HPP_ */
