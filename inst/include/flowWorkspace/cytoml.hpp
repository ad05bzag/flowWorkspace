#include "wsNode.hpp"
using namespace cytolib;

enum class WS_TYPE {WS_WIN, WS_MAC, WS_MAC_3, WS_VX};

typedef unordered_map<string, WS_TYPE> WS_VERSION_TYPE_MAP;
extern WS_VERSION_TYPE_MAP ws_ver_type_map; //declare global var
#define CYTOML_INIT() \
	WS_VERSION_TYPE_MAP ws_ver_type_map = {\
											{"1.61", WS_TYPE::WS_WIN},\
											{"1.6", WS_TYPE::WS_WIN},\
											{"2.0", WS_TYPE::WS_MAC},\
											{"3.0", WS_TYPE::WS_MAC_3},\
											{"1.8", WS_TYPE::WS_VX},\
											{"20.0", WS_TYPE::WS_VX}\
											};\

namespace cytoml
{

	inline void print_supported_workspace_version()
	{
		vector<string> s_ws_type = {"WS_WIN", "WS_MAC", "WS_MAC_3", "WS_VX"};
		PRINT("Workspace versions vs types:");
		for(auto & it : ws_ver_type_map)
			PRINT(it.first + ": " + s_ws_type[static_cast<int>(it.second)] + "\n");
	}

	inline void add_workspace_version(const string & version, WS_TYPE ws_type)
	{
		ws_ver_type_map["version"] = ws_type;
	}

	inline WS_TYPE get_workspace_type(xmlDocPtr doc)
	{
		wsNode root(doc->children);
		xmlXPathObjectPtr res = root.xpath("/Workspace");
		wsNode curNode(res->nodesetval->nodeTab[0]);
		xmlXPathFreeObject(res);
		string version=curNode.getProperty("version");

		auto it = ws_ver_type_map.find(version);
		if(it==ws_ver_type_map.end())
			throw(domain_error("Unrecognized workspace version: " + version \
									+ "\n Please use 'print_supported_workspace_version()' to list the supported version." \
									+ "\n Or use 'add_workspace_version' to add the new version number into the list (if you are sure it is compatible to the existing workspace type)" \
									));

		return it->second;
	}

};
