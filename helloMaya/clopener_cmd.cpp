#include "clopener_cmd.h"

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

    MString command =
        "global proc clopenerToggleMode() {"
        "int $mode = `radioButtonGrp -q -select clopenerSelectMode`;"
        "if($mode == 1){"
        "textFieldButtonGrp -e -enable true clopenerMeshField;"
        "textFieldButtonGrp -e -enable false clopenerFaceField;"
        "} else {"
        "textFieldButtonGrp -e -enable false clopenerMeshField;"
        "textFieldButtonGrp -e -enable true clopenerFaceField;"
        "}"
        "};"

        "if (`window -exists clopenerWindow`) deleteUI clopenerWindow;"
        "window -title \"Clopener\" -widthHeight 400 500 clopenerWindow;"

        // Outer padding layout
        "columnLayout -adjustableColumn true -rowSpacing 12 "
        "-columnAlign \"center\" "
        "-columnAttach \"both\" 20;"

        // Selection mode
        "text -label \"Selection Mode\" -align \"center\";"
        "radioButtonGrp "
        "-numberOfRadioButtons 2 "
        "-labelArray2 \"Mesh\" \"Faces\" "
        "-select 1 "
        "-changeCommand \"clopenerToggleMode\" "
        "-columnAlign2 \"center\" \"center\" "
        "clopenerSelectMode;"

        // Mesh input
        "text -label \"Input Mesh\" -align \"center\";"
        "textFieldButtonGrp "
        "-label \"Mesh\" "
        "-buttonLabel \"Use Selected\" "
        "-buttonCommand \"string \\$sel[] = `ls -sl -type transform`; "
        "if(size(\\$sel)>0) textFieldButtonGrp -e -text \\$sel[0] clopenerMeshField;\" "
        "-columnAlign3 \"center\" \"center\" \"center\" "
        "clopenerMeshField;"

        // Face selection
        "text -label \"Selected Faces\" -align \"center\";"
        "textFieldButtonGrp "
        "-label \"Faces\" "
        "-buttonLabel \"Use Selected\" "
        "-buttonCommand \"string \\$faces[] = `filterExpand -sm 34`; "
        "if(size(\\$faces)==0){ warning \\\"No faces selected\\\"; } "
        "string \\$mesh = \\\"\\\"; "
        "int \\$valid = 1; "
        "for(\\$f in \\$faces){ "
        "string \\$parts[]; tokenize \\$f \\\".\\\" \\$parts; "
        "if(\\$mesh==\\\"\\\") \\$mesh = \\$parts[0]; "
        "if(\\$parts[0] != \\$mesh) \\$valid = 0; "
        "} "
        "if(!\\$valid){ warning \\\"Faces must belong to the same mesh\\\"; } "
        "else{ "
        "string \\$faceStr = stringArrayToString(\\$faces, \\\", \\\"); "
        "textFieldButtonGrp -e -text \\$faceStr clopenerFaceField; "
        "}\" "
        "-columnAlign3 \"center\" \"center\" \"center\" "
        "clopenerFaceField;"

        // Operation
        "text -label \"Operation\" -align \"center\";"
        "radioButtonGrp "
        "-numberOfRadioButtons 2 "
        "-labelArray2 \"Opening\" \"Closing\" "
        "-select 1 "
        "-columnAlign2 \"center\" \"center\" "
        "clopenerOp;"

        // Radius
        "text -label \"Radius (sensitivity)\" -align \"center\";"
        "floatSliderGrp "
        "-label \"Radius\" "
        "-field true "
        "-minValue 0 "
        "-maxValue 10 "
        "-value 0 "
        "-columnAlign3 \"center\" \"center\" \"center\" "
        "clopenerRadius;"

        // Target Edge Length
        "text -label \"Target Edge Length\" -align \"center\";"
        "floatSliderGrp "
        "-label \"Edge Length\" "
        "-field true "
        "-minValue 0.0 "
        "-maxValue 10.0 "
        "-value 1.0 "
        "-columnAlign3 \"center\" \"center\" \"center\" "
        "clopenerEdgeLength;"

        // Iterations
        "text -label \"Iterations\" -align \"center\";"
        "intSliderGrp "
        "-label \"Iterations\" "
        "-field true "
        "-minValue 1 "
        "-maxValue 300 "
        "-value 5 "
        "-columnAlign3 \"center\" \"center\" \"center\" "
        "clopenerIterations;"

        // Spacer
        "separator -style \"none\" -height 10;"

        // Bottom buttons
        "rowLayout "
        "-numberOfColumns 2 "
        "-adjustableColumn 2 "
        "-columnAttach2 \"both\" \"both\" "
        "-columnOffset2 0 0 "
        "-columnWidth2 180 180;"

        "button -label \"Cancel\" -command \"deleteUI clopenerWindow\";"
        "button -label \"Apply\" -command \"clopener\";"

        "setParent ..;"
        "clopenerToggleMode();"
        "showWindow clopenerWindow;";

    MGlobal::executeCommand(command);

    return MS::kSuccess;
}