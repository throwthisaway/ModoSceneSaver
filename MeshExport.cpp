// #pragma warning(disable: 4786)
#undef min
#undef max
#include <lxu_scene.hpp>
#include <lxidef.h>
#include <lxu_select.hpp>
#include <lxu_math.hpp>
#include <lxu_format.hpp>
#include <lx_layer.hpp>
#include <lx_item.hpp>
#include <lx_shade.hpp>

#include <map>
#include <string>
#include <iostream>
using namespace std;
enum class TextureTypes{kAlbedo, kNormal, kMetallic, kRoughness, kCount};
struct Texture {
    unsigned image, uv;
};
struct Material {
	unsigned count;
	string name;
    float rgb[3];
    float metallic, roughess;
	bitset<32> textureMask;
	unsigned uvCount;
	Texture textures[(int)TextureTypes::kCount];
};
struct LineFormat : CLxLineFormat {
	const char * lf_Separator  () override	{ return " | "; }
	void WritePNTS(unsigned count) {
		lf_Output("PNTS");lf_Delim();lf_Output((unsigned)(count * sizeof(float) * 3));lf_Delim();lf_Output(count); lf_Break();
	}
	void WritePoint(float * vec) {
		lf_Output(vec[0]);lf_Delim();lf_Output(vec[1]);lf_Delim();lf_Delim();lf_Output(vec[2]);lf_Break();
	}
	void WritePOLY(unsigned count) {
		lf_Output("POLY");lf_Delim();lf_Output((unsigned)(count * sizeof(unsigned) * 3));lf_Delim();lf_Output(count); lf_Break();
	}
	void WritePolygon(unsigned pnt0, unsigned pnt1, unsigned pnt2) {
		lf_Output(pnt0);lf_Delim();lf_Output(pnt1);lf_Delim();lf_Delim();lf_Output(pnt2);
	}
	void WriteTexture(const Texture& texture) {
		lf_Delim();lf_Output(texture.image);lf_Delim();lf_Output(texture.uv);
	}
	void WriteMaterials(const vector<Material>& materials) {
		for (auto& material : materials) {
			lf_Output(material.name.c_str());
			lf_Delim();
			lf_Output(material.rgb[0]);
			lf_Delim();
			lf_Output(material.rgb[1]);
			lf_Delim();
			lf_Output(material.rgb[2]);
			lf_Delim();
			lf_Output(material.metallic);
			lf_Delim();
			lf_Output(material.roughess);
			lf_Delim();
			lf_Output(material.textureMask.to_string().c_str());
			WriteTexture(material.textures[(int)TextureTypes::kAlbedo]);
			WriteTexture(material.textures[(int)TextureTypes::kNormal]);
			WriteTexture(material.textures[(int)TextureTypes::kMetallic]);
			WriteTexture(material.textures[(int)TextureTypes::kRoughness]);
			lf_Break();
		}
	}
	void WriteImages(const std::vector<string>& images) {
		for(const auto& str : images) {
			lf_Output(str.c_str()); lf_Break();
		}
	}
};
class MeshExport : public CLxSceneSaver {
public:
	static LXtTagInfoDesc descInfo[];
	//CLxBinaryFormat fileFormat;
    LineFormat fileFormat;
	CLxFileFormat* ss_Format() override { return &fileFormat; }

	void ss_Verify() override;
	LxResult ss_Save() override;
	void ss_Point() override;
	void ss_Polygon() override;

//    void EnumRender(CLxUser_Item& item);
    void GatherMaterials();
    bool GatherTexture(Texture& texture, const char* fx, unsigned& uvCount);

	LXtMatrix		 xfrm;
	LXtVector		 xfrmPos;
	map<LXtPointID, unsigned> pnt_index;
	struct Point {
		float vec[3];
	};
	vector<Point> points;
	map<string, unsigned> imageMap;
	std::vector<string> images;
	unsigned	polyCount;
	string prevMaterial;
	map<string, unsigned> materialMap;
	vector<Material> materials;
	map<string, unsigned> uvMap;
	std::vector<string> uvs;
	bool			 get_matr;
};

//void MeshExport::EnumRender(CLxUser_Item& render) {
//    CLxUser_Item item;
//    bool res = item.GetSubItem(0, item);
//    int i = 0;
//    while (res) {
//        item.
//        auto ch = ChanString(LXsICHAN_MASK_PTAG);
//        fileFormat.lf_Output(ch);
//        const char* name, *uniqueName;
//        item2.Name(&name);
//        item2.UniqueName(&uniqueName);
//        fileFormat.lf_Output("\t");
//        fileFormat.lf_Output(uniqueName);
//        //item2.Type();
//        fileFormat.lf_Output("|");
//        if (name) fileFormat.lf_Output(name);
//        fileFormat.lf_Output("|");
//        if (item2.GetTag(LXi_PTAG_MATR)) fileFormat.lf_Output(item2.GetTag(LXi_PTAG_MATR));
//        fileFormat.lf_Break();
//        EnumRender2(item2);
//        res = item.GetSubItem(++i, item2);
//    }
//}

void MeshExport::ss_Verify() {
	Message("common", 2020);
	MessageArg(1, "Export mesh from scene");
}
bool MeshExport::GatherTexture(Texture& texture, const char* name, unsigned& uvCount) {
    const char* fx = LayerEffect();
	CLxUser_Item item;
    if (fx && !strcmp(fx, name)) {
		GetItem(item);
        if (TxtrLocator()) {
            auto uv = ChanString(LXsICHAN_TEXTURELOC_UVMAP);
			if (uv) {
				auto res = uvMap.insert(make_pair(uv, (unsigned)uvs.size()));
				if (res.second) {
					texture.uv = uvCount;
					uvs.push_back(res.first->first);
					++uvCount;
				}
			}
            assert(ChanInt(LXsICHAN_TEXTURELOC_PROJTYPE) == LXi_TEXTURE_PROJ_MODE_UVMAP);
			SetItem(item);	// restore from texture locator item
        }
		const char* fname;
		if (TxtrImage () && (fname = ChanString(LXsICHAN_VIDEOSTILL_FILENAME))) {
			string str = fname;
			auto pos = str.rfind('/');
			if (pos == string::npos) pos = str.rfind('\\');
			if (pos == string::npos) pos = 0; else ++pos;
			str = str.substr(pos, string::npos);
			auto res = imageMap.insert(make_pair(str, (unsigned)images.size()));
			texture.image = res.first->second;
			if (res.second) images.push_back(res.first->first);
			SetItem(item);	// restore from texture image item
			return true;
		}
    }
	return false;
}
void MeshExport::GatherMaterials() {
    // same goes for other tag types like: Part, Selection set
    //!!! auto material_tag_of_the_polygon = ChanString(LXsICHAN_MASK_PTAG);
    for (auto& kv : materialMap) {
        if (!ScanMask(kv.first.c_str())) continue;
        while (NextLayer()) {
			auto& material = materials[kv.second];
            if (!strcmp(ItemType(), "unrealShader")) {
                material.rgb[0] = (float)ChanFloat(LXsICHAN_UNREALMATERIAL_BASECOL".R");
                material.rgb[1] = (float)ChanFloat(LXsICHAN_UNREALMATERIAL_BASECOL".G");
                material.rgb[2] = (float)ChanFloat(LXsICHAN_UNREALMATERIAL_BASECOL".B");
                material.metallic = (float)ChanFloat(LXsICHAN_UNREALMATERIAL_METALLIC);
                material.roughess = (float)ChanFloat(LXsICHAN_UNREALMATERIAL_ROUGH);
            }
            if (ChanInt(LXsICHAN_TEXTURELAYER_ENABLE)) {
                if (ItemIsA(LXsITYPE_IMAGEMAP)) {
					if (GatherTexture(material.textures[(int)TextureTypes::kAlbedo], "baseUE", material.uvCount)) material.textureMask.set((size_t)TextureTypes::kAlbedo);
					if (GatherTexture(material.textures[(int)TextureTypes::kNormal], "normalUE", material.uvCount)) material.textureMask.set((size_t)TextureTypes::kNormal);
					if (GatherTexture(material.textures[(int)TextureTypes::kMetallic], "metallicUE", material.uvCount)) material.textureMask.set((size_t)TextureTypes::kMetallic);
					if (GatherTexture(material.textures[(int)TextureTypes::kRoughness], "roughUE", material.uvCount)) material.textureMask.set((size_t)TextureTypes::kRoughness);
				}
			}
        }
    }
}

LxResult MeshExport::ss_Save() {
	unsigned					npts;
	CLxUser_SelectionService	selSvc;
	double						currTime;

	currTime = selSvc.GetTime();
	/*
	 * Count points in all meshes.
	 */
	npts = 0;
	StartScan();
	while (NextMesh())
		npts += PointCount();

	/*
	 * Tabulate the colors of all the materials.
	 */
	get_matr = true;
	StartScan();
	while (NextMesh())
		WritePolys();

    GatherMaterials();

	///*
	// * Write the sync line and the number of points.
	// */
	//lf_Output("3DG1");
	//lf_Break();
	//lf_Output(npts);
	//lf_Break();

	/*
	 * Write point positions.
	 */
	fileFormat.WriteMaterials(materials);
	fileFormat.WriteImages(images);
	fileFormat.WritePNTS((unsigned)points.size());
	StartScan();
	while (NextMesh())
	{
//		// Get world transformation
//		if (!WorldXform(xfrm, xfrmPos))
//		{
//			lx::MatrixIdent(xfrm);
//			std::fill(xfrmPos, xfrmPos + 3, 0.0);
//		}

		SetMeshTime(currTime);
		WritePoints();
	}
	for (auto& p : points) {
		fileFormat.WritePoint(p.vec);
	}
	/*
	 * Write polygons.
	 */
	get_matr = false;
	fileFormat.WritePOLY(polyCount);
	polyCount = 0;
	prevMaterial.clear();
	StartScan();
	while (NextMesh())
	{
		SetMeshTime(currTime);
		WritePolys(0, true);
	}

//    bool hasUvs;
//    hasUvs = SetSelMap (LXi_VMAP_TEXTUREUV);
//    if (!hasUvs) {
//        hasUvs = SetMap (LXi_VMAP_TEXTUREUV);
//    }
//    if (hasUvs) {
//        for (auto& kv : materials) {
//            auto res = SetMap (LXi_VMAP_TEXTUREUV, kv.second.uv.c_str());
//            int i = 0;
//        }
//    }
	/*
	 * Clear any persistent state.
	 */
	materials.clear();
	pnt_index.clear();
	points.clear();
	materialMap.clear();
	images.clear();
	imageMap.clear();
	uvMap.clear();
	uvs.clear();
	return LXe_OK;
}
void MeshExport::ss_Polygon() {
	if (get_matr) {
		unsigned n = PolyNumVerts();
		if (n == 3) {
			++polyCount;
			auto res = pnt_index.insert(make_pair(PolyVertex(0), points.size()));
			if (res.second) points.push_back({});
			res = pnt_index.insert(make_pair(PolyVertex(1), points.size()));
			if (res.second) points.push_back({});
			res = pnt_index.insert(make_pair(PolyVertex(2), points.size()));
			if (res.second) points.push_back({});
			const char *mask = PolyTag(LXi_PTAG_MATR);
			if (mask) {
				auto res = materialMap.emplace(mask, (unsigned)materials.size());
				if (res.second) materials.push_back({0, mask});
				else ++materials[res.first->second].count;
			}
		}
		return;
	}

	auto material = materialMap.find(PolyTag(LXi_PTAG_MATR));
	assert(material != materialMap.end());
	if (prevMaterial != material->first) {
		fileFormat.lf_Output(material->second);
		fileFormat.lf_Delim();
		fileFormat.lf_Output(material->first.c_str());
		fileFormat.lf_Delim();
		fileFormat.lf_Output(materials[material->second].count);
		fileFormat.lf_Break();
		prevMaterial = material->first;
	}
	unsigned n = PolyNumVerts();
	if (n != 3) return;
	++polyCount;
	fileFormat.WritePolygon(pnt_index[PolyVertex(0)], pnt_index[PolyVertex(1)], pnt_index[PolyVertex(2)]);
	for (auto& material : materials) {
		std::vector<bool> uvWritten(material.uvCount);
		for (unsigned tex = 0; tex < (unsigned)TextureTypes::kCount; ++tex) {
			if (material.textureMask.test(tex)) {
				for (unsigned uv = 0; uv < material.uvCount; ++uv) {
					if(!uvWritten[uv] && SetMap(LXi_VMAP_TEXTUREUV, uvs[uv].c_str())) {
						if (n != 3) continue;
						uvWritten[uv] = true;
						for (unsigned i = 0; i < n; i++) {
							LXtPointID vrt = PolyVertex(i);
							float f[2];
							if (PolyMapValue(f, vrt)) {
								fileFormat.lf_Output(f[0]);
								fileFormat.lf_Delim();
								fileFormat.lf_Output(f[1]);
								fileFormat.lf_Delim();
							}
						}
					}
				}
			}
		}
	}
	fileFormat.lf_Break();
}
void MeshExport::ss_Point() {
	auto res = pnt_index.find(PntID());
	if (res == pnt_index.end()) return;
	float			 opos[3];
//	float			 vec[3];

	PntPosition(opos);

//	lx::MatrixMultiply(vec, xfrm, opos);
//	vec[0] += xfrmPos[0];
//	vec[1] += xfrmPos[1];
//	vec[2] += xfrmPos[2];
	points[res->second].vec[0] = opos[0];points[res->second].vec[1] = opos[1];points[res->second].vec[2] = opos[2];
}

LXtTagInfoDesc	 MeshExport::descInfo[] = {
		{ LXsSAV_OUTCLASS,	LXa_SCENE	},
		{ LXsSAV_DOSTYPE,	"mesh"		},
		{ LXsSRV_USERNAME,	"Mesh Export"},
		{ 0 }
};
void initialize() {
	LXx_ADD_SERVER(Saver, MeshExport, "vs_Mesh");
}

