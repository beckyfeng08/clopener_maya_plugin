#include "clopenercmd.h"

#include <maya/MArgList.h>
#include <maya/MGlobal.h>
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

    MString cmd = "print(\"piped through code\");";
    MGlobal::executeCommand(cmd);


   
    return MS::kSuccess;
}