#pragma once
#include <maya/MPxNode.h>
#include <maya/MPlug.h>
#include <maya/MDataBlock.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnStringData.h>
#include <maya/MFnMeshData.h>
#include <maya/MFnMesh.h>
//#include "clopener.h"
#include <maya/MPoint.h>


class clopener_node : public MPxNode
{
	clopener_node() {};
	virtual ~clopener_node() {};

	virtual MStatus compute(const MPlug& plug, MDataBlock& data);
	static void* creator();
	static MStatus initialize();

	static MObject inputMesh;
	static MObject numIters;
	static MObject edgelen;
	static MObject close;
	static MObject outputMesh;
	static MTypeId id;

private:
	//Clopener l;

};

