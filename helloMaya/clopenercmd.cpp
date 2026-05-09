#include "clopenercmd.h"
#include <maya/MGlobal.h>
#include <maya/MStringArray.h>
#include <maya/MFnSet.h>
#include <maya/MComputation.h>
#include <maya/MProgressWindow.h>
#include <maya/MFnAnimCurve.h>
#include <maya/MPlug.h>
#include <maya/MPlugArray.h>
#include <maya/MTime.h>
#include <maya/MFnDependencyNode.h>

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
	bool useRelative = false;
	double radius = 0.12;
	bool isOpen = false; // true = open operation, false = closing operation
	Eigen::VectorXi verts;
	MString mesh = "";

	// parsing args

	for (unsigned int i = 0; i < argList.length(); i++)
	{
		MString arg = argList.asString(i, &status);

		if (arg == "-i") {
			iterations = argList.asInt(i + 1, &status);
		}
		else if (arg == "-e") {
			edgelen = argList.asDouble(i + 1, &status);
		}
		else if (arg == "-er") {
			int val = argList.asInt(i + 1, &status);
			if (val == 1) { // relative is chosen
				useRelative = true; 
			}
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
		else if (arg == "-m") {
			mesh = argList.asString(i + 1, &status);
		}
		else if (arg == "-v") {
			MString vertStr = argList.asString(i + 1, &status);

			MStringArray tokens;
			vertStr.split(' ', tokens);

			std::vector<int> temp;

			for (unsigned int k = 0; k < tokens.length(); ++k)
			{
				if (tokens[k].length() > 0)
					temp.push_back(tokens[k].asInt());
			}

			verts = Eigen::Map<Eigen::VectorXi>(temp.data(), temp.size());
		}
	}


	MComputation computation; // this is here to check for user interruptions (pressing esc key)
	computation.beginComputation();

	// --- Setup progress bar
	MProgressWindow::reserve();
	MProgressWindow::setTitle("Running Clopener...");
	MProgressWindow::setProgressRange(0, iterations);
	MProgressWindow::setInterruptable(true);
	MProgressWindow::setProgressStatus("Press Escape to stop...");

	if (!MProgressWindow::startProgress())
	{
		MGlobal::displayError("Could not start progress window.");
		computation.endComputation();
		return MS::kFailure;
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

	// setting params
	ClosingFlowParams params; // create inputs to closing flow
	params.maxiter = iterations;
	params.use_relative = useRelative;
	params.frac_of_avg_edge = edgelen;
	params.h = edgelen;
	params.bd = 1.0 / radius;
	params.opening = isOpen;
	params.selection = verts;
	
	Eigen::MatrixXd Vout;
	Eigen::MatrixXi Fout;

	// ACTUAL CLOSING FLOW CALL
	bool closed = false;

	ClosingFlow cf = ClosingFlow(V, F, params);
	for (int i = 0; i < iterations; i++) {

		// --- Check for user cancel (ESC)
		if (computation.isInterruptRequested() || MProgressWindow::isCancelled())
		{
			MGlobal::displayWarning("Operation cancelled by user.");
			status = MS::kFailure;
			break;
		}

		if (!cf.step()) {
			std::cerr << "Flow converged at iteration " << i << "\n";
			closed = true;
			break;
		}

		// show mesh, create new maya mesh object essentially
		Vout = cf.current_V();
		Fout = cf.current_F();
		status = createNewMesh(Vout, Fout, dagPath, i);

		// refresh the GUI, and the loading bar
		MGlobal::executeCommand("refresh -f");
		MProgressWindow::setProgress(i);
		if (computation.isInterruptRequested() || MProgressWindow::isCancelled())
		{
			MGlobal::displayWarning("Operation cancelled by user.");
			status = MS::kFailure;
			break;
		}
	}
	

	V = std::move(Vout);
	F = std::move(Fout);

	MGlobal::displayInfo("Vertices: " + MString() + V.rows()); // remeshed
	MGlobal::displayInfo("Triangles: " + MString() + F.rows());
	
	// --- Cleanup (VERY IMPORTANT)
	MProgressWindow::endProgress();
	computation.endComputation();

	if (!closed) {
		std::cerr << "closing_flow is not complete\n";
		return MS::kFailure;
	}

	//TODO: key visibility (one mesh, per frame)?

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
clopenercmd::createNewMesh(Eigen::MatrixXd V, Eigen::MatrixXi F, const MDagPath& dagPath, int frame) { // for some reason this doesn't output a mesh onto the viewport?
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
	// key animation

	MFnDagNode dagNode(newMeshObj);

	MString shapeName = dagNode.name();
	MString fullShapePath = dagNode.fullPathName();
	if (frame > 0) {
		keyVisibility(fullShapePath, frame - 1, false);
	}
	keyVisibility(fullShapePath, frame, true);
	keyVisibility(fullShapePath, frame + 1, false);

	MGlobal::displayInfo("Shape just created: " + dagNode.name());
	return status;
}

void clopenercmd::keyVisibility(const MString& meshName,
	int frame,
	bool visible)
{
	MString cmd =
		"setKeyframe -attribute visibility -t " +
		MString() + frame +
		" -v " +
		MString() + (visible ? 1 : 0) +
		" \"" +
		meshName +
		"\"";

	MGlobal::executeCommand(cmd);

	cmd.format(
		"keyTangent -itt step -ott step \"^1s\"",
		meshName
	);

	MGlobal::executeCommand(cmd);
}
