// #pragma warning(disable: 4786)
#undef min
#undef max
#include <lxu_scene.hpp>
#include <lxidef.h>
#include <lxu_select.hpp>
#include <lxu_math.hpp>
#include <lxu_format.hpp>

#include <map>
#include <string>
using namespace std;


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



	LXtMatrix		 xfrm;
	LXtVector		 xfrmPos;
	map<LXtPointID, unsigned> pnt_index;
	unsigned		 pnt_count;
	map<string, unsigned>	 matr_color;
	bool			 get_matr;
};

void MeshExport::ss_Verify() {
	Message("common", 2020);
	MessageArg(1, "Export mesh from scene");
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

	//GatherColors();

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

	/*
	 * Clear any persistent state.
	 */
	matr_color.clear();
	pnt_index.clear();
	return LXe_OK;
}
void MeshExport::ss_Polygon() {
	if (get_matr) {
		const char *mask = PolyTag(LXi_PTAG_MATR);
		if (mask)
			matr_color[mask] = 0;
		return;
	}

	unsigned		 i, n;
	char			 buf[32];

	n = PolyNumVerts();
	fileFormat.lf_Output(n);

	for (i = 0; i < n; i++)
		fileFormat.lf_Output(pnt_index[PolyVertex(i)]);

	sprintf(buf, "0x%06X", matr_color[PolyTag(LXi_PTAG_MATR)]);
	fileFormat.lf_Output(buf);
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
    fileFormat.lf_Output(vec[1]);
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

