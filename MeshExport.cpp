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

struct Texture {
    std::string fname;
    std::string uv;
};
struct Material {
    float rgb[3];
    float metallic, roughess;
    struct {
        Texture color, metallic, roughness, normal;
    }textures;
};
class MeshExport : public CLxSceneSaver {
public:
	static LXtTagInfoDesc descInfo[];
	//CLxBinaryFormat fileFormat;
    CLxLineFormat fileFormat;
	CLxFileFormat* ss_Format() override { return &fileFormat; }

	void ss_Verify() override;
	LxResult ss_Save() override;
	void ss_Point() override;
	void ss_Polygon() override;

//    void EnumRender(CLxUser_Item& item);
    void GatherMaterials();
    void WriteMaterials();
    void GatherTexture(Texture& texture, const char* fx);
    void WriteTexture(const Texture& texture);

	LXtMatrix		 xfrm;
	LXtVector		 xfrmPos;
	map<LXtPointID, unsigned> pnt_index;
	unsigned		 pnt_count;
	map<string, Material>	 materials;
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
void MeshExport::GatherTexture(Texture& texture, const char* name) {
    const char* fx = LayerEffect();
	CLxUser_Item item;
    if (fx && !strcmp(fx, name)) {
		GetItem(item);
        if (TxtrLocator()) {
            auto uv = ChanString(LXsICHAN_TEXTURELOC_UVMAP);
            if (uv) texture.uv = uv;
            assert(ChanInt(LXsICHAN_TEXTURELOC_PROJTYPE) == LXi_TEXTURE_PROJ_MODE_UVMAP);
			SetItem(item);	// restore from texture locator item
//            CLxUser_Item map;
//            GetItem(map);
//            const char* fname;
//            if (map.test () && SetItem(map) &&
//                TxtrImage () && (fname = ChanString (LXsICHAN_VIDEOSTILL_FILENAME)))
//                texture.fname = fname;
            // TODO:: texture file name

        }
		const char* fname;
		if (TxtrImage () && (fname = ChanString (LXsICHAN_VIDEOSTILL_FILENAME))) {
			texture.fname = fname;
			SetItem(item);	// restore from texture image item
		}
    }
}
void MeshExport::GatherMaterials() {
    // same goes for other tag types like: Part, Selection set
    //!!! auto material_tag_of_the_polygon = ChanString(LXsICHAN_MASK_PTAG);
    for (auto& kv : materials) {
        if (!ScanMask(kv.first.c_str())) continue;
        while (NextLayer()) {
            if (!strcmp(ItemType(), "unrealShader")) {
                kv.second.rgb[0] = (float)ChanFloat(LXsICHAN_UNREALMATERIAL_BASECOL".R");
                kv.second.rgb[1] = (float)ChanFloat(LXsICHAN_UNREALMATERIAL_BASECOL".G");
                kv.second.rgb[2] = (float)ChanFloat(LXsICHAN_UNREALMATERIAL_BASECOL".B");
                kv.second.metallic = (float)ChanFloat(LXsICHAN_UNREALMATERIAL_METALLIC);
                kv.second.roughess = (float)ChanFloat(LXsICHAN_UNREALMATERIAL_ROUGH);
            }
            if (ChanInt(LXsICHAN_TEXTURELAYER_ENABLE)) {
                if (ItemIsA(LXsITYPE_IMAGEMAP)) {
                    GatherTexture(kv.second.textures.color, "baseUE");
                    GatherTexture(kv.second.textures.metallic, "metallicUE");
                    GatherTexture(kv.second.textures.normal, "normalUE");
                    GatherTexture(kv.second.textures.roughness, "roughUE");
                }
                       // auto fname = ChanString(LXsICHAN_VIDEOSTILL_FILENAME);
                        //                        if (name) kv.second.uv = name;
//                        GetItem(cmap);
//                        cmap.GetUniqueName(kv.second.albedoUnique);
//                        const char * name;
//                        cmap.Name(&name);
            }
//            CLxUser_Item item;
//            if (GetItem(item)) {
//                CLxUser_Item parent;
//                if (item.GetParent(parent)) {
//                    parent.
//                }
//            }
        }
    }
}
void MeshExport::WriteTexture(const Texture& texture) {
    if (!texture.uv.empty()) {
        fileFormat.lf_Output("\t\t");
        fileFormat.lf_Output(texture.uv.c_str());
    }
    if (!texture.fname.empty()) {
        fileFormat.lf_Output(" | ");
        fileFormat.lf_Output(texture.fname.c_str());
    }
}
void MeshExport::WriteMaterials() {
    for (auto& kv : materials) {
        fileFormat.lf_Output(kv.first.c_str());
        fileFormat.lf_Output("\t");
        fileFormat.lf_Output(kv.second.rgb[0]);
        fileFormat.lf_Output(",");
        fileFormat.lf_Output(kv.second.rgb[1]);
        fileFormat.lf_Output(",");
        fileFormat.lf_Output(kv.second.rgb[2]);
        fileFormat.lf_Output("\t");
        fileFormat.lf_Output(kv.second.metallic);
        fileFormat.lf_Output("\t");
        fileFormat.lf_Output(kv.second.roughess);
        WriteTexture(kv.second.textures.color);
        WriteTexture(kv.second.textures.metallic);
        WriteTexture(kv.second.textures.roughness);
        WriteTexture(kv.second.textures.normal);
        fileFormat.lf_Break();
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
    WriteMaterials();
	pnt_count = 0;
	StartScan();
	while (NextMesh())
	{
		// Get world transformation
		if (!WorldXform(xfrm, xfrmPos))
		{
			lx::MatrixIdent(xfrm);
			std::fill(xfrmPos, xfrmPos + 3, 0.0);
		}

		SetMeshTime(currTime);
		WritePoints();
	}

	/*
	 * Write polygons.
	 */
	get_matr = false;
	StartScan();
	while (NextMesh())
	{
		SetMeshTime(currTime);
		WritePolys();
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
	return LXe_OK;
}
void MeshExport::ss_Polygon() {
	if (get_matr) {
		const char *mask = PolyTag(LXi_PTAG_MATR);
		if (mask)
            materials[mask] = {};
		return;
	}

	unsigned		 i, n;
	char			 buf[32];

	n = PolyNumVerts();
	fileFormat.lf_Output(n);

	for (i = 0; i < n; i++) {
		fileFormat.lf_Output(pnt_index[PolyVertex(i)]);
        fileFormat.lf_Output(",");
	}

	//sprintf(buf, "0x%06X", materials[PolyTag(LXi_PTAG_MATR)]);
	fileFormat.lf_Output(PolyTag(LXi_PTAG_MATR));
	fileFormat.lf_Output(",");
	for (auto& kv : materials) {
		if (!kv.second.textures.color.uv.empty()) {
			if(SetMap(LXi_VMAP_TEXTUREUV, kv.second.textures.color.uv.c_str())) {
				fileFormat.lf_Output(kv.second.textures.color.uv.c_str());
				fileFormat.lf_Output(":");
				for (i = 0; i < n; i++) {
					LXtPointID vrt = PolyVertex(i);
					float f[2];
					if (PolyMapValue(f, vrt)) {
						fileFormat.lf_Output(f[0]);
						fileFormat.lf_Output(",");
						fileFormat.lf_Output(f[1]);
						fileFormat.lf_Output(";");
					}
				}

			}
		}
	}
	fileFormat.lf_Break();
}
void MeshExport::ss_Point() {
	float			 opos[3];
	float			 vec[3];

	PntPosition(opos);

	lx::MatrixMultiply(vec, xfrm, opos);
	vec[0] += xfrmPos[0];
	vec[1] += xfrmPos[1];
	vec[2] += xfrmPos[2];

    fileFormat.lf_Output(vec[0]);
	fileFormat.lf_Output(",");
    fileFormat.lf_Output(vec[1]);
	fileFormat.lf_Output(",");
    fileFormat.lf_Output(vec[2]);
    fileFormat.lf_Break();

	pnt_index[PntID()] = pnt_count++;
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

