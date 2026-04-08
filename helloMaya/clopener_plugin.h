#pragma once
#include <maya/MArgList.h>
#include <maya/MObject.h>
#include <maya/MGlobal.h>
#include <maya/MPxCommand.h>
#include <maya/MSyntax.h>

// custom Maya command
class clopener : public MPxCommand
{
public:
	clopener() {};
	static void* creator();
	static MSyntax newSyntax();

};
