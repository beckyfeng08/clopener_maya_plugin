#include "clopenercmd.h"

#include <maya/MArgList.h>
#include <maya/MGlobal.h>
#include <maya/MStringArray.h>
#include <maya/MFnSet.h>
#include <list>
#include <sstream>
#include <iostream>

#include "closing_flow.h"

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

	// set stuff to have stuff
	int iterations = 5;
	double edgelen = 0.03;
	double radius = 0.12;
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
			isOpen = false;
			int val = argList.asInt(i + 1, &status);

			if (val == 1) {
				isOpen = true;

			}
		}
		else if (arg == "-d") {
			isMesh = false;

			if (argList.asInt(i + 1, &status) == 1) {
				isMesh = true;
			}
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

	// ACTUAL CLOSING FLOW CALL

	ClosingFlowParams params; // create inputs to closing flow
	params.maxiter = iterations;
	params.h = edgelen;
	params.bd = 1.0 /radius;
	params.opening = isOpen;

	Eigen::MatrixXd Vout;
	Eigen::MatrixXi Fout;
	// 

	bool old = false;
	bool closed = false;
	if (old) {
		closed = closing_flow(V, F, params, Vout, Fout); // running closing operations on original vertices and faces

		if (!closed) {
			std::cerr << "closing_flow failed\n";
			return MS::kFailure;
		}

		V = std::move(Vout);
		F = std::move(Fout);

		MGlobal::displayInfo("Vertices: " + MString() + V.rows()); // remeshed
		MGlobal::displayInfo("Triangles: " + MString() + F.rows());
		status = createNewMesh(V, F, dagPath);
	} else {
		ClosingFlow cf = ClosingFlow(V, F, params);
		for (int i = 0; i < iterations; i++) {
			if (!cf.step()) {
				std::cerr << "Flow converged at iteration " << i << "\n";
				closed = true;
				break;
			}
			// show mesh, create new maya mesh object essentially
			Vout = cf.current_V();
			Fout = cf.current_F();
			status = createNewMesh(Vout, Fout, dagPath);
		}

		if (!closed) {
			std::cerr << "closing_flow failed\n";
			return MS::kFailure;
		}

		V = std::move(Vout);
		F = std::move(Fout);

		MGlobal::displayInfo("Vertices: " + MString() + V.rows()); // remeshed
		MGlobal::displayInfo("Triangles: " + MString() + F.rows());
	}
	
	
	MGlobal::displayInfo("end");

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

MStatus
clopenercmd::createNewMesh(Eigen::MatrixXd V, Eigen::MatrixXi F, const MDagPath& dagPath) { // for some reason this doesn't output a mesh onto the viewport?
	MStatus status;

	
	// MESH CREATION GIVEN EIGEN VALUES
	MPointArray points;
	points.setLength(V.rows());

	for (int i = 0; i < V.rows(); i++) {
		points[i] = MPoint(V(i, 0), V(i, 1), V(i, 2));
	}
	MIntArray polygonCounts;
	MIntArray polygonConnects;

	int numFaces = F.rows();

	polygonCounts.setLength(numFaces);
	polygonConnects.setLength(numFaces * 3);

	for (int i = 0; i < numFaces; i++) {
		polygonCounts[i] = 3; // triangles

		polygonConnects[3 * i + 0] = F(i, 0);
		polygonConnects[3 * i + 1] = F(i, 1);
		polygonConnects[3 * i + 2] = F(i, 2);
	}

	// create the new mesh object
	MObject newMeshObj = MFnMesh().create(
		V.rows(),
		numFaces,
		points,
		polygonCounts,
		polygonConnects
	);

	MFnMesh newMeshFn(newMeshObj); // rebind

	// assign new mesh to initialShadingGroup for it to be visible
	MObjectArray sets, comps;
	MFnMesh oldMeshFn(dagPath);
	oldMeshFn.getConnectedSetsAndMembers(0, sets, comps, true);

	if (sets.length() > 0) {
		MFnSet fnSet(sets[0]);
		fnSet.addMember(newMeshObj);
	}
	
	return MS::kSuccess;
}