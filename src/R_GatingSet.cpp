/*
 * R_GatingSet.cpp
 *
 *these are R APIs
 *
 *  Created on: Mar 30, 2012
 *      Author: wjiang2
 */

#include "flowWorkspace/openWorkspace.hpp"
#include <Rcpp.h>
using namespace Rcpp;
using namespace cytolib;
using namespace cytoml;
CYTOLIB_INIT()
CYTOML_INIT()

GatingSet * getGsPtr(SEXP _gsPtr){

	if(R_ExternalPtrAddr(_gsPtr)==0)
			throw(domain_error("Null GatingSet pointer!"));
	XPtr<GatingSet>gs(_gsPtr);

	return gs;
}
/*
 * can't use module for exposing overloaded methods
 */


//[[Rcpp::export(name=".cpp_getSamples")]]
StringVec get_sample_uids(XPtr<GatingSet> gsPtr) {

	return gsPtr->get_sample_uids();
}

/*
 * constructing GatingSet from existing gating hierarchy and new data
 */
//[[Rcpp::export(name=".cpp_NewGatingSet")]]
XPtr<GatingSet> NewGatingSet(XPtr<GatingSet> gsPtr
               ,string src_sample_uid
               ,StringVec new_sample_uids)
  {

		GatingHierarchy & gh=gsPtr->getGatingHierarchy(src_sample_uid);

		/*
		 * used gh as the template to clone multiple ghs in the new gs
		 */
		GatingSet * newGS=new GatingSet(gh,new_sample_uids);

		/*
		 * using default finalizer to delete gs,which is triggered by gc() when
		 * xptr is out of scope
		 */

		return XPtr<GatingSet>(newGS);

}

/*
 * constructing GatingSet with only root node for each sample
 */
//[[Rcpp::export(name=".cpp_NewGatingSet_rootOnly")]]
XPtr<GatingSet> NewGatingSet_rootOnly(StringVec new_sample_uids) {

		GatingSet * newGS=new GatingSet(new_sample_uids);

		return XPtr<GatingSet>(newGS);

}

/*
 * save/load GatingSet
 */
//[[Rcpp::export(name=".cpp_saveGatingSet")]]
void saveGatingSet(XPtr<GatingSet> gs, string fileName) {
			gs->serialize_pb(fileName);
}


//[[Rcpp::export(name=".cpp_loadGatingSet")]]
XPtr<GatingSet> loadGatingSet(string fileName) {
		GatingSet * gs=new GatingSet(fileName);
		return XPtr<GatingSet>(gs);

}


//[[Rcpp::export(name=".cpp_CloneGatingSet")]]
XPtr<GatingSet> CloneGatingSet(XPtr<GatingSet> gs,StringVec new_sample_uids) {

		GatingSet * gs_new=gs->clone(new_sample_uids);

		return XPtr<GatingSet>(gs_new);

}

//[[Rcpp::export(name=".cpp_combineGatingSet")]]
XPtr<GatingSet> combineGatingSet(Rcpp::List gsList,Rcpp::List sampleList) {

		GatingSet * newGS=new GatingSet();

		for(int i=0;i<gsList.size();i++)
		{
			GatingSet *	gs=getGsPtr((SEXP)gsList[i]);
			StringVec samples=as<StringVec>(sampleList[i]);
			newGS->add(*gs,samples);
		}


		return XPtr<GatingSet>(newGS);

}

/**
 * change sample name
 */
//[[Rcpp::export(name=".cpp_setSample")]]
void set_sample_uid(XPtr<GatingSet> gs,string oldName, string newName) {
	
		gs->setSample(oldName,newName);

}

