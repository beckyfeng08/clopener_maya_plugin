//#include "clopener_plugin.h"
#include <maya/MFnPlugin.h>
#include <maya/MGlobal.h>
#include <maya/MDGModifier.h>
#include <maya/MPlugArray.h>
#include <maya/MVector.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MStringArray.h>
#include <list>

#include "clopenercmd.h"

// define EXPORT for exporting dll functions
#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif// Maya Plugin creator function

static const char* kMenuName = "clopenerMenu";

// Initialize Maya Plugin upon loading
EXPORT MStatus initializePlugin(MObject obj)
{
    MStatus status;
    MFnPlugin plugin(obj, "Clopener", "1.0", "Any");

    // register command
    status = plugin.registerCommand(
        "clopenercmd",
        clopenercmd::creator);

    if (!status) {
        status.perror("registerCommand failed");
        return status;
    }


    MGlobal::displayInfo("clopenercmd loaded");

    MString cmd;
    cmd = R"(
global proc clopenerToggleVertexField()
{
    int $useVerts = `checkBox -q -value clopenerUseVerts`;
    textFieldButtonGrp -e -enable $useVerts clopenerVertexField;
}

global proc clopenerToggleEdgeMode()
{
    int $mode = `radioButtonGrp -q -select clopenerEdgeMode`;

    // 1 = Absolute, 2 = Relative
    floatSliderGrp -e -enable ($mode == 1) clopenerEdgeLengthAbs;
    floatSliderGrp -e -enable ($mode == 2) clopenerEdgeLengthRel;
}

global proc clopenerFillVertices()
{
    string $verts[] = `filterExpand -sm 31`;

    if(size($verts) == 0){
        warning "No vertices selected";
        return;
    }

    string $mesh = "";
    int $valid = 1;

    string $indices[];

    for($v in $verts){
        string $parts[];
        tokenize $v "." $parts;

        if($mesh == "")
            $mesh = $parts[0];

        if($parts[0] != $mesh)
            $valid = 0;

        // Extract index inside [ ]
        string $idxParts[];
        tokenize $parts[1] "[]" $idxParts;

        // idxParts[1] should be the number
        $indices[size($indices)] = $idxParts[1];
    }

    if(!$valid){
        warning "Vertices must belong to the same mesh";
        return;
    }

    string $displayStr = stringArrayToString($indices, " ");
    textFieldButtonGrp -e -text $displayStr clopenerVertexField;
}

global proc openWindow()
{
    if (`window -exists clopenerWindow`)
        deleteUI clopenerWindow;

    window -title "Clopener" -widthHeight 400 500 clopenerWindow;

    columnLayout -adjustableColumn true -rowSpacing 12
        -columnAlign "center"
        -columnAttach "both" 20;

    // Mesh input (REQUIRED)
    text -label "Input Mesh" -align "center";
    textFieldButtonGrp
        -label "Mesh"
        -buttonLabel "Use Selected"
        -buttonCommand "string $sel[]=`ls -sl -type transform`; if(size($sel)>0) textFieldButtonGrp -e -text $sel[0] clopenerMeshField;"
        -columnAlign3 "center" "center" "center"
        clopenerMeshField;

    // Vertex toggle
    checkBox
        -label "Select Vertices"
        -value 0
        -changeCommand "clopenerToggleVertexField() "
        clopenerUseVerts;

    // Vertex input (OPTIONAL)
    textFieldButtonGrp
        - label "Vertices"
        - buttonLabel "Use Selected"
        - buttonCommand "clopenerFillVertices() "
        - enable 0
        - columnAlign3 "center" "center" "center"
        clopenerVertexField;

    // Operation
    text - label "Operation" - align "center";
    radioButtonGrp
        - numberOfRadioButtons 2
        - labelArray2 "Opening" "Closing"
        - select 1
        - columnAlign2 "center" "center"
        clopenerOp;

    // Radius
    text - label "Radius (Intensity of Effect) " - align "center";
    floatSliderGrp
        - label "Radius"
        - field true
        - minValue 0.001
        - maxValue 10
        - value 0.2
        - precision 3
        - columnAlign3 "center" "center" "center"
        clopenerRadius;

    // Target Edge Length
text -label "Target Edge Length (Coarseness) " - align "center";

    // Mode selection
    radioButtonGrp
        - numberOfRadioButtons 2
        - labelArray2 "Absolute" "Relative (Fraction of Average Edge Length) "
        - select 1
        - changeCommand "clopenerToggleEdgeMode() "
        - columnAlign2 "center" "center"
        clopenerEdgeMode;

    // Absolute edge length
    floatSliderGrp
        - label "Absolute Length"
        - field true
        - minValue 0.001
        - maxValue 2.0
        - value 0.5
        - precision 3
        - columnAlign3 "center" "center" "center"
        clopenerEdgeLengthAbs;

    // Relative edge length (fraction)
    floatSliderGrp
        - label "Relative Fraction"
        - field true
        - minValue 0.01
        - maxValue 1.0
        - value 0.5
        - precision 2
        - enable 0
        - columnAlign3 "center" "center" "center"
        clopenerEdgeLengthRel;

    // Iterations
    text - label "Iterations" - align "center";
    intSliderGrp
        - label "Iterations"
        - field true
        - minValue 1
        - maxValue 100
        - value 5
        - columnAlign3 "center" "center" "center"
        clopenerIterations;

    separator - style "none" - height 10;

    // Buttons
    rowLayout
        - numberOfColumns 2
        - adjustableColumn 2
        - columnAttach2 "both" "both"
        - columnOffset2 0 0
        - columnWidth2 180 180;

    button - label "Cancel" - command "deleteUI clopenerWindow";
    button - label "Apply" - command "processData() ";

    setParent ..;

    showWindow clopenerWindow;
}

global proc processData()
{
    int $iterations = `intSliderGrp - q - value clopenerIterations`;
        int $edgeMode = `radioButtonGrp - q - select clopenerEdgeMode`;

        float $edgelen;
    int $useRelative = 0;

    if ($edgeMode == 1) {
        $edgelen = `floatSliderGrp - q - value clopenerEdgeLengthAbs`;
            $useRelative = 0;
    }
    else {
        $edgelen = `floatSliderGrp - q - value clopenerEdgeLengthRel`;
            $useRelative = 1;
    }
    float $radius = `floatSliderGrp - q - value clopenerRadius`;
        int $op = `radioButtonGrp - q - select clopenerOp`;
        int $useVerts = `checkBox - q - value clopenerUseVerts`;

        string $mesh = `textFieldButtonGrp - q - text clopenerMeshField`;
        string $verts = "";

    if ($mesh == "") {
        error "Mesh is required.";
        return;
    }

    if ($useVerts) {
        $verts = `textFieldButtonGrp - q - text clopenerVertexField`;
            if ($verts == "") {
                warning "Vertex mode enabled but no vertices provided.";
            }
    }

    print("Mesh: " + $mesh + "\n");
    print("Use Vertices: " + $useVerts + "\n");
    print("Vertices: " + $verts + "\n");

    clopenercmd
        - i $iterations
        - e $edgelen
        - er $useRelative
        - r $radius
        - o $op
        - m $mesh
        - v $verts;
}

global proc addClopenerMenu()
{
    global string $gMainWindow;

    if (!`menu - exists clopenerMenu`)
        menu - label "Clopener" - parent $gMainWindow clopenerMenu;

        if (!`menuItem - exists clopenerMenuItem`)
            menuItem
            - label "Clopener Window"
            - command "openWindow() "
            - parent clopenerMenu
            clopenerMenuItem;
}

addClopenerMenu();
    )";


    MGlobal::executeCommand(cmd);
    MGlobal::displayInfo("Menu item successfully generated");

    return status;
}
// Cleanup Plugin upon unloading
EXPORT MStatus uninitializePlugin(MObject obj)
{
    MStatus status;
    MFnPlugin plugin(obj);

    status = plugin.deregisterCommand("clopenercmd");
    if (!status) {
        status.perror("deregisterCommand failed");
        return status;
    }



    MString cmd = "if (`menu -exists  clopenerMenu`) deleteUI clopenerMenu;";

    MGlobal::executeCommand(cmd);
    MGlobal::displayInfo("clopenercmd unloaded");

    return status;
}