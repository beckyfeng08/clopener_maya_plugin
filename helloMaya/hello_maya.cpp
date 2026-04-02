#include "hello_maya.h"
#include <maya/MFnPlugin.h>
#include <maya/MGlobal.h>
#include <igl/massmatrix.h>
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
// Plugin doIt function
MStatus clopener::doIt(const MArgList& argList)
{
	MStatus status;


	MString command =
        "if (`window -exists clopenerWindow`) deleteUI clopenerWindow;"
        "window -title \"Clopener\" -widthHeight 360 240 clopenerWindow;"

        // Outer padding layout
        "columnLayout -adjustableColumn true -rowSpacing 12 "
        "-columnAlign \"center\" "
        "-columnAttach \"both\" 20;"

        // Mesh input
        "text -label \"Input Mesh\" -align \"center\";"
        "textFieldButtonGrp "
        "-label \"Mesh\" "
        "-buttonLabel \"Select\" "
        "-columnAlign3 \"center\" \"center\" \"center\" "
        "clopenerMeshField;"

        // Operation
        "text -label \"Operation\" -align \"center\";"
        "radioButtonGrp "
        "-numberOfRadioButtons 2 "
        "-labelArray2 \"Opening\" \"Closing\" "
        "-select 1 "
        "-columnAlign2 \"center\" \"center\" "
        "clopenerOp;"

        // Radius
        "text -label \"Radius\" -align \"center\";"
        "floatSliderGrp "
        "-label \"Radius\" "
        "-field true "
        "-minValue 0 "
        "-maxValue 10 "
        "-value 0 "
        "-columnAlign3 \"center\" \"center\" \"center\" "
        "clopenerRadius;"

        // Spacer
        "separator -style \"none\" -height 10;"

        // Bottom buttons spanning width
        "rowLayout "
        "-numberOfColumns 2 "
        "-adjustableColumn 2 "
        "-columnAttach2 \"both\" \"both\" "
        "-columnOffset2 0 0 "
        "-columnWidth2 180 180;"

        "button -label \"Cancel\" -command \"deleteUI clopenerWindow\";"
        "button -label \"Apply\" -command \"clopener\";"

        "setParent ..;"
        "showWindow clopenerWindow;";
	MGlobal::executeCommand(command);

	return MS::kSuccess;
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