#include "clopener_plugin.h"
#include <maya/MFnPlugin.h>
#include <maya/MGlobal.h>

#include "clopener_cmd.h"
#include "clopener_node.h"

// define EXPORT for exporting dll functions
#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif// Maya Plugin creator function

static const char* kMenuName = "clopenerMenu";

void* clopener::creator()
{
	return new clopener;
}
MSyntax clopener::newSyntax()
{
	MSyntax syntax;
	return syntax;
}

// Initialize Maya Plugin upon loading
EXPORT MStatus initializePlugin(MObject obj)
{
    MStatus status;
    MFnPlugin plugin(obj, "CIS660", "1.0", "Any");

    status = plugin.registerCommand(
        "clopener",
        clopener::creator,
        clopener::newSyntax
    );

    if (!status)
        status.perror("registerCommand failed");

    MString cmd;
    cmd += "global string $gMainWindow;";
    cmd += "if (`menu -exists ";
    cmd += kMenuName;
    cmd += "`) deleteUI ";
    cmd += kMenuName;
    cmd += ";";
    cmd += "menu -label \"Clopener\" -parent $gMainWindow ";
    cmd += kMenuName;
    cmd += ";";
    cmd += "menuItem -label \"Open Window\" -command \"clopener\" -parent ";
    cmd += kMenuName;
    cmd += ";";

    MGlobal::executeCommand(cmd);

    return status;
}
// Cleanup Plugin upon unloading
EXPORT MStatus uninitializePlugin(MObject obj)
{
    MStatus status;
    MFnPlugin plugin(obj);

    status = plugin.deregisterCommand("clopener");
    if (!status)
        status.perror("deregisterCommand failed");

    MString cmd;
    cmd += "if (`menu -exists ";
    cmd += kMenuName;
    cmd += "`) deleteUI ";
    cmd += kMenuName;
    cmd += ";";

    MGlobal::executeCommand(cmd);

    return status;
}