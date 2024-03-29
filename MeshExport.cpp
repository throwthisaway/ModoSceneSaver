// #pragma warning(disable: 4786)
#include "Platform.h"
#ifdef PLATFORM_WIN
#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
// Windows Header Files
#include <windows.h>
#undef max
#undef min
#endif
#include <lxu_scene.hpp>
#include <lxidef.h>
#include <lxu_select.hpp>
#include <lxu_math.hpp>
#include <lxu_format.hpp>
#include <lx_layer.hpp>
#include <lx_item.hpp>
#include <lx_shade.hpp>

#include <algorithm>
#include <map>
#include <unordered_map>
#include <string>
#include <iostream>
#include <sstream>

#include <meshoptimizer.h>
#include <ModoMeshLoader.h>
#include <bitset>

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
	void WriteSubmeshes(const vector<Submesh>& submeshes) {
		const unsigned size = 88;
		ff.lf_Output(Tag(MATR)); ff.lf_Output((unsigned)(submeshes.size() * size));ff.lf_Output((unsigned)submeshes.size());ff.lf_Break();
		for (auto& submesh : submeshes) {
			//ff.lf_Output(material.name.c_str());
			ff.lf_Output(submesh.indexByteOffset); // 4
			ff.lf_Output(submesh.vertexByteOffset); // 8
			ff.lf_Output(submesh.stride); // 12
			ff.lf_Output(submesh.material.diffuse.r); // 16
			ff.lf_Output(submesh.material.diffuse.g); // 20
			ff.lf_Output(submesh.material.diffuse.b); // 24
			ff.lf_Output(submesh.material.metallic_roughness.r); // 28
			ff.lf_Output(submesh.material.metallic_roughness.g); // 32
			ff.lf_Output(submesh.vertexType); // 36
			ff.lf_Output(submesh.textureMask); // 40
			ff.lf_Output(submesh.uvCount); // 44
			WriteTexture(submesh.textures[(int)TextureTypes::kAlbedo]); // 52
			WriteTexture(submesh.textures[(int)TextureTypes::kNormal]); // 60
			WriteTexture(submesh.textures[(int)TextureTypes::kMetallic]); // 68
			WriteTexture(submesh.textures[(int)TextureTypes::kRoughness]); // 76
			WriteTexture(submesh.textures[(int)TextureTypes::kBump]); // 84
			ff.lf_Output(submesh.count); // 88
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
	void WritePolygons(index_t * indices, uint32_t count) {
		for (uint32_t j = 0; j < count; j += 2) {
			ff.lf_Output((indices[j + 1] << 16) | indices[j]);
			if (!((j + 1) % kVertPerPoly)) ff.lf_Break();
		}
	}
	void WriteHeader(const Header& header) {
		ff.lf_Output(Tag(HEAD)); ff.lf_Output((unsigned)sizeof(Header)); ff.lf_Output((unsigned)1); ff.lf_Break();
		ff.lf_Output(header.version);ff.lf_Break();
		ff.lf_Output(header.r);ff.lf_Break();
		ff.lf_Output(header.aabb.min.x);ff.lf_Output(header.aabb.min.y);ff.lf_Output(header.aabb.min.z);ff.lf_Break();
		ff.lf_Output(header.aabb.max.x);ff.lf_Output(header.aabb.max.y);ff.lf_Output(header.aabb.max.z);ff.lf_Break();
	}
};

struct LineFormat : CLxLineFormat {
	const char * lf_Separator  () override	{ return " | "; }
};
struct BinaryFormat : CLxBinaryFormat {
	void lf_Break() {}
};

inline bool HasTangent(const Submesh& submesh) {
	return submesh.textureMask & (1<<(int)TextureTypes::kNormal);
}
inline int CalcElementCount(const Submesh& submesh) {
	int tangent = HasTangent(submesh) ? 3 : 0;
	return 3*2 /*pos + normal*/ + tangent + 2 * submesh.uvCount;
}
void Copy(const char* src, const char * dst) {
	if (!strcmp(src, dst)) return;	// TODO:: proper path comparison
#ifdef PLATFORM_WIN
	bool res = ::CopyFileA(src, dst, FALSE);
	auto err = GetLastError();
	int ii = 0;
#else
	ostringstream cp;
	cp << "/bin/cp -rf \"" << src << "\" \"" << dst << "\"";
	int res = system(cp.str().c_str());
	assert(!res);
#endif
}
size_t FindLastSeparator(const string& str) {
	auto pos = str.rfind('/');
	if (pos == string::npos) pos = str.rfind('\\');
	if (pos == string::npos) pos = 0; else ++pos;
	return pos;
}
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
 
	std::string GetOutDir() {
		std::string out(fileFormat.ff.file_name);
		out = out.substr(0, FindLastSeparator(out));
		return out;
	}
	void SelectMesh(int index) {
		StartScan();
		for (int i = 0; i <= index; ++i) NextMesh();
	}
	glm::vec3 CalcPointTangent(LXtPointID ptID, int submeshIndex);
	glm::vec3 CalcPolyTangent();
	bool Equals(const glm::dvec3& n, const glm::vec2& uv, LXtPointID pt);

	LXtMatrix		 xfrm;
	LXtVector		 xfrmPos;
	struct Context {
		void Serialize(const std::string& outdir) {
			std::string fname(outdir);
			fname.append("\\");
			fname.append(name);
			fname.append(".mesh");
			FileFormat<LineFormat> fileFormat;
			fileFormat.ff.ff_Open(fname.c_str());

			fileFormat.WriteHeader(header);
			fileFormat.WriteSubmeshes(submeshes);
			fileFormat.WriteImages(images);
			fileFormat.WriteVERT((unsigned)vertexByteOffset, (unsigned)vCount);
			for (int i = 0; i < out.size(); ++i) {
				int elementSize = CalcElementCount(submeshes[i]);
				fileFormat.WriteVertices(out[i].vertices.data(), (unsigned)out[i].vertices.size(), elementSize);
			}
			fileFormat.WritePOLY((unsigned)iCount);
			for (int i = 0; i < out.size(); ++i)
				fileFormat.WritePolygons(out[i].indices.data(), (uint32_t)out[i].indices.size());
		}
		void Optimize();
		void GatherSubmeshes(CLxSceneSaver& saver, const std::string& outdir);
		bool GatherTexture(CLxSceneSaver& saver, const std::string& outdir, Texture& texture, const char* fx, unsigned& uvCount);
		map<LXtPointID, LXtFVector> points;

		map<string, unsigned> imageMap;
		std::vector<string> images;

		map<string, unsigned> submeshMap;
		vector<Submesh> submeshes;

		std::vector<std::unordered_map<LXtPointID, std::vector<LXtPolygonID>>> pointsToPolyMapPerSubmesh;
		vector<std::vector<float>> vertexData;	// vertices per material

		map<string, unsigned> uvMap;
		std::vector<string> uvs;

		Header header;
		struct Out {
			std::vector<index_t> indices;
			std::vector<float> vertices;
		};
		size_t iCount = 0, vCount = 0;
		size_t vertexByteOffset = 0;
		std::vector<Out> out;
		std::string name;
	};
	Context context;
	int meshCount = 0;
	bool get_matr;
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
bool MeshExport::Context::GatherTexture(CLxSceneSaver& saver, const std::string& outdir, Texture& texture, const char* name, unsigned& uvCount) {
    const char* fx = saver.LayerEffect();
	CLxUser_Item item;
    if (fx && !strcmp(fx, name)) {
		saver.GetItem(item);
        if (saver.TxtrLocator()) {
            auto uv = saver.ChanString(LXsICHAN_TEXTURELOC_UVMAP);
			if (uv) {
				auto res = uvMap.insert(make_pair(uv, (unsigned)uvs.size()));
				if (res.second) {
					texture.uv = uvCount;
					uvs.push_back(res.first->first);
					++uvCount;
				}
			}
            assert(saver.ChanInt(LXsICHAN_TEXTURELOC_PROJTYPE) == LXi_TEXTURE_PROJ_MODE_UVMAP);
			saver.SetItem(item);	// restore from texture locator item
        }
		const char* fname;
		if (saver.TxtrImage () && (fname = saver.ChanString(LXsICHAN_VIDEOSTILL_FILENAME))) {
			string str = fname;
			auto pos = FindLastSeparator(str);
			str = str.substr(pos, string::npos);
			Copy(fname, (outdir + str).c_str());
			auto res = imageMap.insert(make_pair(str, (unsigned)images.size()));
			texture.id = res.first->second;
			if (res.second) images.push_back(res.first->first);
			saver.SetItem(item);	// restore from texture image item
			return true;
		}
    }
	return false;
}
void MeshExport::Context::GatherSubmeshes(CLxSceneSaver& saver, const std::string& outdir) {
    // same goes for other tag types like: Part, Selection set
    //!!! auto material_tag_of_the_polygon = ChanString(LXsICHAN_MASK_PTAG);
    for (auto& kv : submeshMap) {
        if (!saver.ScanMask(kv.first.c_str())) continue;
		bitset<32> mask;
		auto& submesh = submeshes[kv.second];
        while (saver.NextLayer()) {
            if (!strcmp(saver.ItemType(), "unrealShader")) {
                submesh.material.diffuse[0] = (float)saver.ChanFloat(LXsICHAN_UNREALMATERIAL_BASECOL".R");
                submesh.material.diffuse[1] = (float)saver.ChanFloat(LXsICHAN_UNREALMATERIAL_BASECOL".G");
                submesh.material.diffuse[2] = (float)saver.ChanFloat(LXsICHAN_UNREALMATERIAL_BASECOL".B");
                submesh.material.metallic_roughness.x = (float)saver.ChanFloat(LXsICHAN_UNREALMATERIAL_METALLIC);
                submesh.material.metallic_roughness.y = (float)saver.ChanFloat(LXsICHAN_UNREALMATERIAL_ROUGH);
            }
            if (saver.ChanInt(LXsICHAN_TEXTURELAYER_ENABLE)) {
                if (saver.ItemIsA(LXsITYPE_IMAGEMAP)) {
					if (GatherTexture(saver, outdir, submesh.textures[(int)TextureTypes::kAlbedo], "baseUE", submesh.uvCount)) mask.set((size_t)TextureTypes::kAlbedo);
					if (GatherTexture(saver, outdir, submesh.textures[(int)TextureTypes::kNormal], "normalUE", submesh.uvCount)) mask.set((size_t)TextureTypes::kNormal);
					if (GatherTexture(saver, outdir, submesh.textures[(int)TextureTypes::kMetallic], "metallicUE", submesh.uvCount)) mask.set((size_t)TextureTypes::kMetallic);
					if (GatherTexture(saver, outdir, submesh.textures[(int)TextureTypes::kRoughness], "roughUE", submesh.uvCount)) mask.set((size_t)TextureTypes::kRoughness);
					if (GatherTexture(saver, outdir, submesh.textures[(int)TextureTypes::kBump], "bumpUE", submesh.uvCount)) mask.set((size_t)TextureTypes::kBump);
				}
			}
        }
		submesh.textureMask = (unsigned)mask.to_ulong();
    }
}
void MeshExport::Context::Optimize() {
	size_t indexByteOffset = 0;
	for (int i = 0; i < vertexData.size(); ++i) {
		std::vector<index_t> indices;
		std::vector<float> vertices;
		auto& vd = vertexData[i];
		auto& submesh = submeshes[i];
		int elementSize = CalcElementCount(submesh);
		int byteStride = sizeof(float) * elementSize;
		const size_t elementCount = vd.size() / elementSize;
		std::vector<uint32_t> remap(elementCount);
		size_t vertexCount = meshopt_generateVertexRemap(remap.data(), nullptr, elementCount, vd.data(), elementCount, byteStride);
		indices.resize(elementCount);
		meshopt_remapIndexBuffer<index_t>(indices.data(), nullptr, elementCount, remap.data());
		vertices.resize(vertexCount * elementSize);
		meshopt_remapVertexBuffer(vertices.data(), vd.data(), elementCount, byteStride, remap.data());
		submesh.indexByteOffset = (uint32_t)indexByteOffset; submesh.vertexByteOffset = (uint32_t)vertexByteOffset;
		submesh.stride = byteStride;
		submesh.count = (int)indices.size();
		indexByteOffset += indices.size() * sizeof(indices[0]);
		vertexByteOffset += vertices.size() * sizeof(vertices[0]);
		iCount += indices.size();
		vCount += vertices.size();
		out.push_back({ std::move(indices), std::move(vertices) });
	}
}
LxResult MeshExport::ss_Save() {
	CLxUser_SelectionService	selSvc;
	double						currTime;

	currTime = selSvc.GetTime();
	meshCount = 0;
	StartScan();
	while (NextMesh()) ++meshCount;

	for (int i = 0; i < meshCount; ++i) {
		SelectMesh(i);
		unsigned npts = PointCount();
		context = Context();
		context.name = ItemName();
		/*
		 * build material map, count points to polygons
		 */
		get_matr = true;
		WritePolys();

		SetMeshTime(currTime);
		WritePoints();

		context.GatherSubmeshes(*this, GetOutDir());

		SelectMesh(i);
		get_matr = false;
		SetMeshTime(currTime);
		WritePolys(0, true);

		// meshoptimizer
		context.Optimize();

		// actual serialization
		context.Serialize(GetOutDir());

		//vertexData.clear();
		//submeshes.clear();
		//points.clear();
		//submeshMap.clear();
		//images.clear();
		//imageMap.clear();
		//uvMap.clear();
		//uvs.clear();
	}
	return LXe_OK;
}

bool MeshExport::Equals(const glm::dvec3& n, const glm::vec2& uv, LXtPointID ptID) {
	glm::dvec3 otherN;
	PolyNormal((double*)&otherN, ptID);
	glm::vec2 otherUV;
	PolyMapValue((float*)&otherUV, ptID);
	return n == otherN && uv == otherUV;
}

glm::vec3 MeshExport::CalcPointTangent(LXtPointID ptID, int submeshIndex) {
	const auto& polys = context.pointsToPolyMapPerSubmesh[submeshIndex][ptID];
	glm::vec3 result{0};
	glm::dvec3 n;
	PolyNormal((double*)&n, ptID);
	glm::vec2 uv;
	PolyMapValue((float*)&uv, ptID);

	int count = 0;
	for (const auto& poly : polys) {
		PolySet(poly);
		if (Equals(n, uv, PolyVertex(0)) || Equals(n, uv, PolyVertex(1)) || Equals(n, uv, PolyVertex(2))) {
			result += CalcPolyTangent();
			++count;
		}
	}
	return result / (float)count;
}
glm::vec3 MeshExport::CalcPolyTangent() {
	glm::vec3 result;
	glm::vec2 uv0, uv1, uv2;
	auto res = PolyMapValue((float*)&uv0, PolyVertex(0));
	assert(res);
	res = PolyMapValue((float*)&uv1, PolyVertex(1));
	assert(res);
	res = PolyMapValue((float*)&uv2, PolyVertex(2));
	assert(res);
	auto pt = context.points[PolyVertex(0)];
	glm::vec3 p0 = glm::vec3(pt[0], pt[1], pt[2]);
	pt = context.points[PolyVertex(1)];
	glm::vec3 p1 = glm::vec3(pt[0], pt[1], pt[2]);
	pt = context.points[PolyVertex(2)];
	glm::vec3 p2 = glm::vec3(pt[0], pt[1], pt[2]);
	auto e0 = p1 - p0, e1 = p2 - p0;
	auto duv0 = uv1 - uv0, duv1 = uv2 - uv0;
	float f = 1.f / (duv0.x * duv1.y - duv0.y * duv1.x);
	result = glm::normalize(f * (duv1.y * e0 - duv0.y * e1));
	return result;

}
void MeshExport::ss_Polygon() {
	if (get_matr) {
		if (PolyNumVerts() != kVertPerPoly) return;
		const char *mask = PolyTag(LXi_PTAG_MATR);
		if (mask) {
			auto res = context.submeshMap.emplace(mask, (unsigned)context.submeshes.size());
			if (res.second) context.submeshes.push_back({});
			//else ++materials[res.first->second].count;

			// record points to polygons per submesh
			if (res.second) context.pointsToPolyMapPerSubmesh.push_back({});
			
			auto& pointsToPolyMap = context.pointsToPolyMapPerSubmesh[res.first->second];
			for (int i = 0; i < kVertPerPoly; ++i) {
				auto res = pointsToPolyMap.emplace(PolyVertex(i), std::vector<LXtPolygonID>{});
				res.first->second.push_back(PolyID());
			}
		}
		return;
	}
	assert(!context.submeshes.empty() && "no material with triangles");
	context.vertexData.resize(context.submeshes.size());
	auto& vertices = context.vertexData[context.submeshMap[PolyTag(LXi_PTAG_MATR)]];
	int submeshIndex = context.submeshMap[PolyTag(LXi_PTAG_MATR)];
	auto& submesh = context.submeshes[submeshIndex];
	if (PolyNumVerts() != kVertPerPoly) return;
	submesh.vertexType = 1 << (int)VertexFields::kNormal;
	if (HasTangent(submesh)) submesh.vertexType |= 1 << (int)VertexFields::kTangent;
	for (unsigned uv = 0; uv < submesh.uvCount; ++uv) submesh.vertexType |= 1 << ((int)VertexFields::kUV0 + uv);
	LXtVector n;
	for (int i = 0; i < kVertPerPoly ; ++i) {
		auto& pt = context.points[PolyVertex(i)];
		vertices.push_back(pt[0]);vertices.push_back(pt[1]);vertices.push_back(pt[2]);
		if (PolyNormal(n, PolyVertex(i))) {
			vertices.push_back((float)n[0]);vertices.push_back((float)n[1]);vertices.push_back((float)n[2]);
		} else {
			vertices.push_back(255.f); vertices.push_back(255.f); vertices.push_back(255.f);
		}
		// if has normal map, calculate and store tangents
		if (HasTangent(submesh)) {
			unsigned normalUVMapIndex = submesh.textures[(int)TextureTypes::kNormal].uv;
			if(!SetMap(LXi_VMAP_TEXTUREUV, context.uvs[normalUVMapIndex].c_str())) { vertices.push_back(255.f);vertices.push_back(255.f);vertices.push_back(255.f);}
			
			auto polyID = PolyID(); //save
			auto tangent = CalcPointTangent(PolyVertex(i), submeshIndex);
			vertices.push_back(tangent.x); vertices.push_back(tangent.y); vertices.push_back(tangent.z);
			PolySet(polyID); // restore
		}
		// distinct uv map indices
		std::set<int> uvIndices;
		for (int i = 0; i < (int)TextureTypes::kCount; ++i)
			if (submesh.textureMask & (1 << i)) uvIndices.insert(submesh.textures[i].uv);

		std::vector<int> sortedUvIndices;
		for (auto uv : uvIndices) sortedUvIndices.push_back(uv);
		// the order uv coordinates references uv maps by their indices
		std::sort(std::begin(sortedUvIndices), std::end(sortedUvIndices));
		for (auto uv : sortedUvIndices) {
			if(!SetMap(LXi_VMAP_TEXTUREUV, context.uvs[uv].c_str())) {
				vertices.push_back(255.f); vertices.push_back(255.f);
				continue;
			}
			float f[2];
			if (PolyMapValue(f, PolyVertex(i))) {
				vertices.push_back(f[0]); vertices.push_back(f[1]);
			} else {
				vertices.push_back(255.f); vertices.push_back(255.f);
			}
		}
	}
}
void MeshExport::ss_Point() {
	LXtFVector& pt = context.points[PntID()];
	auto& header = context.header;
	PntPosition(pt);
	header.r = std::max(header.r, std::max(std::abs(pt[0]), std::max(std::abs(pt[1]), std::abs(pt[2]))));
	if (pt[0] < header.aabb.min.x) header.aabb.min.x = pt[0];
	if (pt[0] > header.aabb.max.x) header.aabb.max.x = pt[0];
	if (pt[1] < header.aabb.min.y) header.aabb.min.y = pt[1];
	if (pt[1] > header.aabb.max.y) header.aabb.max.y = pt[1];
	if (pt[2] < header.aabb.min.z) header.aabb.min.z = pt[2];
	if (pt[2] > header.aabb.max.z) header.aabb.max.z = pt[2];
	
}

LXtTagInfoDesc	 MeshExport::descInfo[] = {
		{ LXsSAV_OUTCLASS,	LXa_SCENE	},
		{ LXsSAV_DOSTYPE,	"scene"		},
		{ LXsSRV_USERNAME,	"Mesh Scene Export"},
		{ 0 }
};
void initialize() {
	LXx_ADD_SERVER(Saver, MeshExport, "vs_MeshScene");
}

