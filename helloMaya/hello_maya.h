#ifndef HELLOMAYA_H
#define HELLOMAYA_H
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
	virtual MStatus doIt(const MArgList& args);
	static void* creator();
	static MSyntax newSyntax();

};
#endif