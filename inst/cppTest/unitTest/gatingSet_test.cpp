/*
 * gatingSet_test.cpp
 *
 *  Created on: May 15, 2012
 *      Author: wjiang2
 */

#include "test_header.hpp"
//#include <boost/dynamic_bitset.hpp>
//#include <boost/serialization/bitset.hpp>
//#include <bitset>
//#include <boost/serialization/array.hpp>
#include <boost/math/distributions/normal.hpp>

/*
	 * plot gating hierarchy tree
	 */


void plotGraph(GatingHierarchy& gh){

			gh.drawGraph("../output/test.dot");
			system("dot2gxl ../output/test.dot -o ../output/test.gxl");
}

void gh_accessor_test(GatingHierarchy& gh){
	/*
		 * getNodes by the T order
		 */

		cout<<endl<<"tsorted node list"<<endl;
		VertexID_vec vertices;
		vertices=gh.getVertices(TSORT);
		for(VertexID_vec::iterator it=vertices.begin();it!=vertices.end();it++)
		{
			nodeProperties &node=gh.getNodeProperty(*it);
			cout<<*it<<"."<<node.getName()<<endl;
		}

		/*
		 * getNodes BFS
		 */
		cout<<endl<<"bfs node list"<<endl;

		vertices=gh.getVertices(BFS);
		for(VertexID_vec::iterator it=vertices.begin();it!=vertices.end();it++)
		{
			nodeProperties &node=gh.getNodeProperty(*it);
			cout<<*it<<"."<<node.getName()<<endl;
		}


		cout<<endl<<"compensation info"<<endl;
		compensation comp=gh.get_compensation();
		cout<<"cid:"<<comp.cid<<endl;
		cout<<"comment:"<<comp.comment<<endl;
		/*
		 * getNodes by vertices ID order
		 * and get stats from each node
		 */

		cout<<endl<<"node list in regular order and stats,gate"<<endl;
		vertices=gh.getVertices(REGULAR);
		for(VertexID_vec::iterator it=vertices.begin();it!=vertices.end();it++)
		{
			VertexID u=*it;
			nodeProperties &node=gh.getNodeProperty(u);
			cout<<u<<"."<<node.getName()<<":";
			cout<<node.getStats(false)["count"]<<endl;
			if(u!=ROOTNODE)
			{
				gate * g=node.getGate();
				cout<<typeid(*g).name()<<endl;


			}
		}

		/*
		 * getPopNames with full path
		 */
		cout<<endl<<"node list with/without full path:"<<endl;
		vector<string> popNames=gh.getPopPaths(REGULAR,true,true);
		for(vector<string>::iterator it=popNames.begin();it!=popNames.end();it++)
			cout<<*it<<endl;
		popNames=gh.getPopPaths(REGULAR,false,true);
		for(vector<string>::iterator it=popNames.begin();it!=popNames.end();it++)
			cout<<*it<<endl;


		/*
		 * get children and parent node index
		 */

		cout<<endl<<"check parent node"<<endl;
		for(size_t i=0;i<vertices.size();i++)
		{
			if(i!=0)
			{
				VertexID parent=gh.getParent(i);
				cout<<i<<"<--"<<parent<<" ";
				cout<<endl;
			}
		}



		cout<<endl<<"check children node"<<endl;
		for(size_t i=0;i<vertices.size();i++)
		{
			VertexID_vec children=gh.getChildren(i);
			cout<<i<<"-->";
			for(VertexID_vec::iterator it=children.begin();it!=children.end();it++)
						cout<<*it<<",";
			cout<<endl;
		}
}

void gs_gating(GatingSet &gs,SampleInfo sample_info, bool is_fix_slash_in_channel_name, bool isH5, string h5_path){
	cout<<endl<<"do the gating after the parsing"<<endl;

	//read transformed data once for all nodes
	GatingHierarchy & gh=gs.getGatingHierarchy(sample_info.uid);

//	gh.loadData(curSample);//get flow data from cdf

	/*
	 * read flow data from cdf into memory first (mimic the R code)
	 */
	FCS_READ_PARAM config;
	config.header.is_fix_slash_in_channel_name = is_fix_slash_in_channel_name;

	CytoFrame * frm;
	if(isH5)
		frm = new H5CytoFrame(sample_info.fcs_path, config, false, h5_path);
	else
		frm = new MemCytoFrame(sample_info.fcs_path, config, false);

	frm->set_pheno_data("name", sample_info.name);

	gh.set_frame_ptr(frm);
	gh.load_fdata_cache();//
	if(sample_info.comp.marker.size()>0)
	{
		sample_info.comp.cid = "1";
		gh.set_compensation(sample_info.comp);
	}
	gh.compensate();
//	gh.adjustGate(gains);
	gh.transform_gate();
	gh.transform_data();
	gh.extendGate(0);
	gh.gating(0,false, true);
	//sync comp & trans data

	gh.release_fdata_cache(true);

}
void gh_counts(GatingHierarchy & gh,vector<bool> &isEqual, const float tolerance, const vector<VertexID> skipPops){
	cout<<endl<<"flowJo(flowcore) counts after gating"<<endl;
	VertexID_vec vertices=gh.getVertices(0);
	for(VertexID_vec::iterator it=vertices.begin();it!=vertices.end();it++)
	{
		VertexID u=*it;
		if(find(skipPops.begin(), skipPops.end(), u) == skipPops.end())//skip some pops that flowJo records the wrong counts
		{
			if(u!=ROOTNODE){
				nodeProperties &node=gh.getNodeProperty(u);
				int flowJoCount = node.getStats(false)["count"];
				if(flowJoCount != -1) //skip the unrecorded flowJo counts
				{

					int myCount = node.getStats(true)["count"];
					cout<<u<<"."<<node.getName()<<":";
					cout<< flowJoCount;
					cout<<"("<<myCount<<") "<< "cv = ";

					bool thisEqual;
					float thisCV ;
					float thisTol = tolerance;
					if(flowJoCount == myCount){
						thisEqual = true;
						thisCV = 0;
					}
					else
					{
						float min = flowJoCount>myCount?myCount:flowJoCount;
						float max = flowJoCount<myCount?myCount:flowJoCount;
						float mean = (min+max)/2;
						float sd = sqrt((pow((min-mean), 2) + pow((max-mean), 2))/2);
						boost::math::normal_distribution<> dist(mean,sd);
						float Q1 = quantile(dist,0.25);
						float Q3 = quantile(dist,0.75);
						float IQR = Q3 - Q1;
						thisCV = IQR/mean;
						thisTol = 1/max+thisTol;//add the weight of N of cells to make it more robust
						thisEqual = (thisCV < thisTol);
		//				cout << Q1 <<":" <<Q3 << " " << thisCV<<endl;
					}
					cout << thisCV << " tol = " << thisTol << endl;
					isEqual.push_back(thisEqual);
				}
			}
		}
	}
}

void gh_removeGate(GatingHierarchy& gh){
	gh.removeNode(5);

}
void clone_test(testCase myTest){
	string archive=myTest.archive;
	GatingSet *gs=new GatingSet(archive);
	gs->clone(gs->get_sample_uids());

}

void parser_test(testCase & myTest){
//	print_supported_workspace_version();
	bool isTemplate = myTest.isTemplate;
	bool isLoadArchive = myTest.isLoadArhive;
	unsigned format = myTest.archiveFormat;
	bool isParseGate = myTest.isParseGate;
	bool isSaveArchive = myTest.isSaveArchive;
	bool archiveType = myTest.archiveType;
	string archiveName = myTest.archive;
	map<string,float> gains = myTest.gains;
	bool is_fix_slash_in_channel_name = false;
	if(archiveType)
		archiveName = archiveName.append(".pb");
	else
		archiveName = archiveName.append(".dat");
//	unsigned short wsType = myTest.wsType;
		boost::scoped_ptr<GatingSet> gs;
		if(isLoadArchive)
		{

			gs.reset(new GatingSet(archiveName));

		}
		else
		{
			//parse a set of sampleIDs
			vector<string> sampleIDs;
			for(map<string,string>::iterator it=myTest.samples.begin();it!=myTest.samples.end();it++)
				sampleIDs.push_back(it->first);
			vector<string> sample_uids;
				for(map<string,string>::iterator it=myTest.samples.begin();it!=myTest.samples.end();it++)
					sample_uids.push_back(it->second);
			if(isTemplate)
				sampleIDs.erase(sampleIDs.begin());//remove the first sample,which is used for testing gating template feature

			if(!isLoadArchive)
			{
				workspace * ws = openWorkspace(myTest.filename, myTest.sampNloc,myTest.xmlParserOption);
				gs.reset(ws->ws2gs(myTest.sample_info,isParseGate));
				is_fix_slash_in_channel_name = ws->is_fix_slash_in_channel_name();
				delete ws;
			}


			cout<<endl<<"get sample names from gating set"<<endl;


		}


		vector<string> samples=gs->get_sample_uids();
		for(vector<string>::iterator it=samples.begin();it!=samples.end();it++)
			cout<<*it<<endl;



		SampleInfo curSample=myTest.sample_info.at(0);

		GatingHierarchy& gh=gs->getGatingHierarchy(curSample.uid);

//		gh_accessor_test(gh);

		if(!isLoadArchive)
			gs_gating(*gs,curSample, is_fix_slash_in_channel_name, myTest.isH5, myTest.h5_path);

		/*
		 * recompute the gate to check if the pop indices are restored properly
		 */
		if(isLoadArchive)
		{
			std::clock_t    start;

			 start = std::clock();
			 FCS_READ_PARAM config;


			gh.set_frame_ptr(new MemCytoFrame(myTest.fcs, config, true));
			gh.load_fdata_cache();//
			gh.transform_data();
			gh.gating(0, true);
			gh.release_fdata_cache(false);

			string filename = "timelog2";
			double runtime = (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
			ofstream timelog;
			timelog.open(filename.c_str(), ios_base::app);
//			myTest.times.push_back(runtime);
			timelog << runtime << ",";
			timelog.close();
			cout << runtime << endl;
		}

		gh_counts(gh, myTest.isEqual, myTest.tolerance, myTest.skipPops);
//		nodeProperties * np = gh.getNodeProperty(102);
//		vector<bool> thisInd = np->getIndices();
//		vector<bool> thisInd1(1024);
//		char thisInd[128];
//		boost::dynamic_bitset<>thisInd(346700);
//		bitset <346700> thisInd;

//		np->setIndices(10);
//		std::ofstream ofs(myTest.archive.c_str(),std::ios::out|std::ios::trunc|std::ios::binary);
//		boost::archive::text_oarchive oa(ofs);
//		boost::archive::xml_oarchive oa(ofs);
//		boost::archive::binary_oarchive oa(ofs);
//		oa << BOOST_SERIALIZATION_NVP(boost::serialization::make_array(&thisInd[0],thisInd.size()));

		if(isSaveArchive){
			if(archiveType == PB)
				gs->serialize_pb(archiveName);
			else
				throw("bs is no longer supported!");
//				gs->serialize_bs(archiveName,myTest.archiveFormat);
		}




}
