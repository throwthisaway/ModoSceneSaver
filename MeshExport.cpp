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
#include <sstream>

#include <meshoptimizer.h>
#include <ModoMeshLoader.h>

using namespace std;
using namespace ModoMeshLoader;
namespace {
template<typename T>
struct FileFormat {
	T ff;
	void WriteVERT(unsigned size, unsigned count) {
		ff.lf_Output(Tag(VERT));ff.lf_Output(size);ff.lf_Output(count);ff.lf_Break();
	}
	void WriteVertices(float * vec, unsigned count, unsigned elementCount) {
		for (unsigned i = 0; i < count; ++i) {
			ff.lf_Output(vec[i]);
			if (!((i + 1) % elementCount)) ff.lf_Break();
		}
		ff.lf_Break();
	}
	void WritePOLY(unsigned count) {
		ff.lf_Output(Tag(POLY));ff.lf_Output((unsigned)(count * sizeof(index_t)));ff.lf_Output(count); ff.lf_Break();
	}
	void WritePolygon(unsigned pnt0, unsigned pnt1, unsigned pnt2) {
		ff.lf_Output(pnt0);ff.lf_Output(pnt1);ff.lf_Output(pnt2);
	}
	void WriteTexture(const Texture& texture) {
		ff.lf_Output(texture.id);ff.lf_Output(texture.uv);
	}
	void WriteMaterials(const vector<Material>& materials) {
		ff.lf_Output(Tag(MATR)); ff.lf_Output((unsigned)(materials.size() * sizeof(materials[0])));ff.lf_Output((unsigned)materials.size());ff.lf_Break();
		for (auto& material : materials) {
			//ff.lf_Output(material.name.c_str());
			ff.lf_Output(material.indexByteOffset);
			ff.lf_Output(material.vertexByteOffset);
			ff.lf_Output(material.stride);
			ff.lf_Output(material.rgb[0]);
			ff.lf_Output(material.rgb[1]);
			ff.lf_Output(material.rgb[2]);
			ff.lf_Output(material.metallic);
			ff.lf_Output(material.roughness);
			ff.lf_Output(material.textureMask);
			ff.lf_Output(material.uvCount);
			WriteTexture(material.textures[(int)TextureTypes::kAlbedo]);
			WriteTexture(material.textures[(int)TextureTypes::kNormal]);
			WriteTexture(material.textures[(int)TextureTypes::kMetallic]);
			WriteTexture(material.textures[(int)TextureTypes::kRoughness]);
			ff.lf_Output(material.count);
			ff.lf_Break();
		}
	}
	void WriteImages(const std::vector<string>& images) {
		unsigned size = 0;
		for(const auto& str : images) {
			size += (unsigned)str.length() + 1; ff.lf_Break();
		}
		ff.lf_Output(Tag(IMAG));ff.lf_Output(size);ff.lf_Output((unsigned)images.size());ff.lf_Break();
		for(const auto& str : images) {
			ff.lf_Output(str.c_str()); ff.lf_Break();
		}
	}
	void WritePolygons(uint16_t * indices, uint16_t count) {
		for (uint16_t j = 0; j < count; j += 2) {
			ff.lf_Output((indices[j + 1] << 16) | indices[j]);
			if (!((j + 1) % kVertPerPoly)) ff.lf_Break();
		}
	}
};

struct LineFormat : CLxLineFormat {
	const char * lf_Separator  () override	{ return " | "; }
};
struct BinaryFormat : CLxBinaryFormat {
	void lf_Break() {}
};
}
class MeshExport : public CLxSceneSaver {
public:
	static LXtTagInfoDesc descInfo[];
	FileFormat<BinaryFormat> fileFormat;
	//FileFormat<LineFormat> fileFormat;
	CLxFileFormat* ss_Format() override { return &fileFormat.ff; }

	void ss_Verify() override;
	LxResult ss_Save() override;
	void ss_Point() override;
	void ss_Polygon() override;

//    void EnumRender(CLxUser_Item& item);
    void GatherMaterials();
    bool GatherTexture(Texture& texture, const char* fx, unsigned& uvCount);
	void Copy(const char* src, const char * dst);
	size_t FindLastSeparator(const string& str);
	LXtMatrix		 xfrm;
	LXtVector		 xfrmPos;
	map<LXtPointID, LXtFVector> points;
	map<string, unsigned> imageMap;
	std::vector<string> images;
	map<string, unsigned> materialMap;
	vector<Material> materials;
	vector<std::vector<float>> vertexData;	// vertices per material
	map<string, unsigned> uvMap;
	std::vector<string> uvs;
	struct Out {
		std::vector<index_t> indices;
		std::vector<float> vertices;
	};
	std::vector<Out> out;
	bool			 get_matr;
	int polyCount = 0;
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
void MeshExport::Copy(const char* src, const char * dst) {
	ostringstream cp;
	cp << "/bin/cp -rf \"" << src << "\" \"" << dst << "\"";
	int res = system(cp.str().c_str());
	assert(!res);
}
size_t MeshExport::FindLastSeparator(const string& str) {
	auto pos = str.rfind('/');
	if (pos == string::npos) pos = str.rfind('\\');
	if (pos == string::npos) pos = 0; else ++pos;
	return pos;
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
			{
				string out(fileFormat.ff.file_name);
				out = out.substr(0, FindLastSeparator(out));
				Copy(fname, out.c_str());
			}
			string str = fname;
			auto pos = FindLastSeparator(str);
			str = str.substr(pos, string::npos);
			auto res = imageMap.insert(make_pair(str, (unsigned)images.size()));
			texture.id = res.first->second;
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
		bitset<32> mask;
		auto& material = materials[kv.second];
        while (NextLayer()) {
            if (!strcmp(ItemType(), "unrealShader")) {
                material.rgb[0] = (float)ChanFloat(LXsICHAN_UNREALMATERIAL_BASECOL".R");
                material.rgb[1] = (float)ChanFloat(LXsICHAN_UNREALMATERIAL_BASECOL".G");
                material.rgb[2] = (float)ChanFloat(LXsICHAN_UNREALMATERIAL_BASECOL".B");
                material.metallic = (float)ChanFloat(LXsICHAN_UNREALMATERIAL_METALLIC);
                material.roughness = (float)ChanFloat(LXsICHAN_UNREALMATERIAL_ROUGH);
            }
            if (ChanInt(LXsICHAN_TEXTURELAYER_ENABLE)) {
                if (ItemIsA(LXsITYPE_IMAGEMAP)) {
					if (GatherTexture(material.textures[(int)TextureTypes::kAlbedo], "baseUE", material.uvCount)) mask.set((size_t)TextureTypes::kAlbedo);
					if (GatherTexture(material.textures[(int)TextureTypes::kNormal], "normalUE", material.uvCount)) mask.set((size_t)TextureTypes::kNormal);
					if (GatherTexture(material.textures[(int)TextureTypes::kMetallic], "metallicUE", material.uvCount)) mask.set((size_t)TextureTypes::kMetallic);
					if (GatherTexture(material.textures[(int)TextureTypes::kRoughness], "roughUE", material.uvCount)) mask.set((size_t)TextureTypes::kRoughness);
				}
			}
        }
		material.textureMask = (unsigned)mask.to_ulong();
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
	polyCount = 0;
	get_matr = true;
	StartScan();
	while (NextMesh())
		WritePolys();

	GatherMaterials();

	StartScan();
	while (NextMesh())
	{
		SetMeshTime(currTime);
		WritePoints();
	}

	get_matr = false;
	polyCount = 0;
	StartScan();
	while (NextMesh())
	{
		SetMeshTime(currTime);
		WritePolys(0, true);
	}

	// meshoptimizer
	int indexByteOffset = 0, vertexByteOffset = 0;
	int iCount = 0, vCount = 0;
	for (int i = 0; i < vertexData.size(); ++i) {
		std::vector<index_t> indices;
		std::vector<float> vertices;
		auto& vd = vertexData[i];
		const auto& material = materials[i];
		int elementSize = 3 * 2 + 2 * material.uvCount;
		int byteStride = sizeof(float)  * elementSize;
		const size_t elementCount = vd.size() / elementSize;
		std::vector<uint32_t> remap(elementCount);
		size_t vertexCount = meshopt_generateVertexRemap(remap.data(), nullptr, elementCount,vd.data(), elementCount, byteStride);
		indices.resize(elementCount);
		meshopt_remapIndexBuffer<index_t>(indices.data(), nullptr, elementCount, remap.data());
		vertices.resize(vertexCount * elementSize);
		meshopt_remapVertexBuffer(vertices.data(), vd.data(), elementCount, byteStride, remap.data());
		materials[i].indexByteOffset = indexByteOffset; materials[i].vertexByteOffset = vertexByteOffset;
		materials[i].stride = byteStride;
		materials[i].count = (int)indices.size();
		indexByteOffset += indices.size() * sizeof(indices[0]);
		vertexByteOffset += vertices.size() * sizeof(vertices[0]);
		iCount += indices.size();
		vCount += vertices.size();
		out.push_back({std::move(indices), std::move(vertices)});
	}

	// actual serialization
	fileFormat.WriteMaterials(materials);
	fileFormat.WriteImages(images);
	fileFormat.WriteVERT(vertexByteOffset, vCount);
	for (int i = 0; i < out.size(); ++i) {
		int elementSize = 3 * 2 + 2 * materials[i].uvCount;
		fileFormat.WriteVertices(out[i].vertices.data(), (unsigned)out[i].vertices.size(), elementSize);
	}
	fileFormat.WritePOLY(iCount);
	for (int i = 0; i < out.size(); ++i)
		fileFormat.WritePolygons(out[i].indices.data(), (uint16_t)out[i].indices.size());

	vertexData.clear();
	materials.clear();
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
		if (PolyNumVerts() != kVertPerPoly) return;
		const char *mask = PolyTag(LXi_PTAG_MATR);
		if (mask) {
			auto res = materialMap.emplace(mask, (unsigned)materials.size());
			if (res.second) materials.push_back({});
			//else ++materials[res.first->second].count;
		}
		return;
	}
	++polyCount;
	vertexData.resize(materials.size());
	auto& vertices = vertexData[materialMap[PolyTag(LXi_PTAG_MATR)]];
	auto& material = materials[materialMap[PolyTag(LXi_PTAG_MATR)]];
	if (PolyNumVerts() != kVertPerPoly) return;
	LXtVector n;
	for (int i = kVertPerPoly - 1; i >=0 ; --i) {
		auto& pt = points[PolyVertex(i)];
		vertices.push_back(pt[0]);vertices.push_back(pt[1]);vertices.push_back(pt[2]);
		if (PolyNormal(n, PolyVertex(i))) {
			vertices.push_back((float)n[0]);vertices.push_back((float)n[1]);vertices.push_back((float)n[2]);
		} else {
			vertices.push_back(0.f); vertices.push_back(0.f); vertices.push_back(0.f);
		}
		for (int uv = 0; uv < material.uvCount; ++uv) {
			if(!SetMap(LXi_VMAP_TEXTUREUV, uvs[uv].c_str())) {
				vertices.push_back(0.f); vertices.push_back(0.f);
				continue;
			}
			float f[2];
			if (PolyMapValue(f, PolyVertex(i))) {
				vertices.push_back(f[0]); vertices.push_back(f[1]);
			} else {
				vertices.push_back(0.f); vertices.push_back(0.f);
			}
		}
	}
}
void MeshExport::ss_Point() {
	PntPosition(points[PntID()]);
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

