#include "clopener_node.h"

MObject clopener_node::inputMesh;
MObject clopener_node::numIters;
MObject clopener_node::edgelen;
MObject clopener_node::close;
MObject clopener_node::outputMesh;
MTypeId clopener_node::id(0x90000);

void* clopener_node::creator()
{
	return new clopener_node;
}

MStatus clopener_node::initialize()
{
	MStatus status;

	return status;
}

MStatus clopener_node::compute(const MPlug& plug, MDataBlock& data) {
	MStatus status;


	return status;
}