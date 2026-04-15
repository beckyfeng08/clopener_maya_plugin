#pragma once
//#include "clopener.h"

#include <maya/MPxCommand.h>
#include <maya/MStatus.h> 
#include <maya/MArgList.h> 
#include <string>

#include <maya/MSelectionList.h>
#include <maya/MFnMesh.h>
#include <maya/MPointArray.h>
#include <maya/MIntArray.h>
#include <maya/MDagPath.h>

#include <Eigen/Dense>

class clopenercmd : public MPxCommand
{
public:
    clopenercmd();
    virtual ~clopenercmd();

    static void* creator() { return new clopenercmd(); }

    // Explicit return type now recognized
    MStatus doIt(const MArgList& args) override;

	Eigen::MatrixXd getMeshVertices(const MDagPath& meshDagPath);
	Eigen::MatrixXi getMeshFaces(const MDagPath& meshDagPath);
    MStatus createNewMesh(Eigen::MatrixXd V, Eigen::MatrixXi F, const MDagPath& meshDagPath);
};