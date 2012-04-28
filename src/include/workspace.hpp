/*
 * workspace.hpp
 *
 *  Created on: Mar 22, 2012
 *      Author: wjiang2
 */

#ifndef WORKSPACE_HPP_
#define WORKSPACE_HPP_
#include <vector>
#include <string>
#include <libxml/xpath.h>
#include "wsNode.hpp"
#include "transformation.hpp"
typedef vector<calibrationTable*> CALTBS;
using namespace std;
/*TODO: so far I will see the differenc between wind and max workspace in terms of xpath(like xpath of sample node)
 * if this is the case eventually we can try to use one template class (eliminate two derived classes )
 * with T structure that stores different versions of xpaths for win/mac,for example:
 *
 * struct winWorkspace{
 * xpath_sample=xxx
 * ....
 * }
 *
 * struct macWorkspace{
 * xpath_sample=xxx
 * ....
 * }
 *
 * this may potentially reduce the amount of code
 *
 */
class compensation{
public:
	string cid;
	string prefix;
	string suffix;
	string comment;// store "Acquisition-defined" when the spillOver matrix is not supplied and cid==-1
	vector<string> marker;
	vector<double> spillOver;
};


struct xpath{
	string group;
	string sampleRef;
	string sample;
	string sampleNode;
	string popNode;
};

class workspace{
public:
	 xpath nodePath;
//protected:

	 xmlDoc * doc;
	 unsigned short dMode;//debug mode passed from gatingset class
public:
	 ~workspace();
	 virtual string xPathSample(string sampleID)=0;
	 virtual Trans_map getTransformation(wsSampleNode,string,CALTBS *)=0;
	 virtual compensation getCompensation(wsSampleNode)=0;
	 virtual CALTBS getCalTbls()=0;
	 virtual vector <string> getSampleID(unsigned short)=0;
	 virtual string getSampleName(wsSampleNode &)=0;
	 virtual wsRootNode getRoot(wsSampleNode sampleNode)=0;
	 virtual wsPopNodeSet getSubPop(wsNode *)=0;
	 virtual gate * getGate(wsPopNode &)=0;//gate is dynamically allocated within this function,it is currently freed within gate pointer owner object nodeProperties
	 virtual nodeProperties * to_popNode(wsRootNode &)=0;
	 virtual nodeProperties * to_popNode(wsPopNode &,bool isGating)=0;
	 valarray<double> toArray(string sCalTable);
};


#endif /* WORKSPACE_HPP_ */

