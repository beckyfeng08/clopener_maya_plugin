#include "clopenercmd.h"

#include <maya/MArgList.h>
#include <maya/MGlobal.h>
#include <maya/MStringArray.h>
#include <list>
#include <sstream>
clopenercmd::clopenercmd() : MPxCommand()
{
}

clopenercmd::~clopenercmd()
{
}

// // Plugin doIt function
MStatus clopenercmd::doIt(const MArgList& argList)
{
    MStatus status;

	int iterations = 1;
	double edgelen = 0.5;
	double radius = 0.5;
	bool isOpen = false; // true = open operation, false = closing operation
	bool isMesh = true; // true = selected mesh, false = selected faces
	MStringArray faces;
	MString mesh = "";


	// parse arguments
	for (unsigned int i = 0; i < argList.length(); i++)
	{
		MString arg = argList.asString(i, &status);

		if (arg == "-i") {
			iterations = argList.asInt(i + 1, &status);
		}
		else if (arg == "-e") {
			edgelen = argList.asDouble(i + 1, &status);
		}
		else if (arg == "-r") {
			radius = argList.asDouble(i + 1, &status);
		}
		else if (arg == "-o") {
			isOpen = (argList.asInt(i + 1, &status) == 1);
		}
		else if (arg == "-d") {
			isMesh = (argList.asInt(i + 1, &status) == 1);
		}
		else if (arg == "-m") {
			mesh = argList.asString(i + 1, &status);
		}
		else if (arg == "-f") {
			MString facesStr = argList.asString(i + 1, &status);

			MStringArray tokens;
			facesStr.split(' ', tokens);

			faces = tokens;
		}
	}

	// clopen algorithm, call clopener

	MSelectionList sel;
	sel.add(mesh);

	MDagPath dagPath;
	sel.getDagPath(0, dagPath);
	//dagPath.extendToShape();
   
	// render mesh with new geometry
	Eigen::MatrixXd V = getMeshVertices(dagPath);
	Eigen::MatrixXi F = getMeshFaces(dagPath);
	MGlobal::displayInfo("Vertices: " + MString() + V.rows());
	MGlobal::displayInfo("Triangles: " + MString() + F.rows());
    
	return MS::kSuccess;
}

Eigen::MatrixXd
clopenercmd::getMeshVertices(const MDagPath& meshDagPath) {
	// create a Eigen::MatrixXd to hold the vertex positions
	MFnMesh fnMesh(meshDagPath);

	MPointArray points;
	fnMesh.getPoints(points, MSpace::kWorld);

	int n = points.length();
	Eigen::MatrixXd V(n, 3);

	for (int i = 0; i < n; i++) {
		V(i, 0) = points[i].x;
		V(i, 1) = points[i].y;
		V(i, 2) = points[i].z;
	}

	return V;
}

Eigen::MatrixXi
clopenercmd::getMeshFaces(const MDagPath& meshDagPath) {
	// create a Eigen::MatrixXi to hold the face indices

	MFnMesh fnMesh(meshDagPath);

	MIntArray triCounts, triIndices;
	fnMesh.getTriangles(triCounts, triIndices);

	// Each triangle has 3 indices
	int numTris = triIndices.length() / 3;

	Eigen::MatrixXi F(numTris, 3);

	for (int i = 0; i < numTris; i++) {
		F(i, 0) = triIndices[3 * i];
		F(i, 1) = triIndices[3 * i + 1];
		F(i, 2) = triIndices[3 * i + 2];
	}

	return F;
}