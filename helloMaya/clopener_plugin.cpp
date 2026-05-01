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
        global proc clopenerToggleMode()
        {
            int $mode = `radioButtonGrp -q -select clopenerSelectMode`;

            if($mode == 1){
                textFieldButtonGrp -e -enable true clopenerMeshField;
                textFieldButtonGrp -e -enable false clopenerFaceField;
            }
            else{
                textFieldButtonGrp -e -enable false clopenerMeshField;
                textFieldButtonGrp -e -enable true clopenerFaceField;
            }
        }

        global proc clopenerFillFaces()
        {
            string $faces[] = `filterExpand -sm 34`;

            if(size($faces) == 0){
                warning "No faces selected";
                return;
            }

            string $mesh = "";
            int $valid = 1;

            for($f in $faces){
                string $parts[];
                tokenize $f "." $parts;

                if($mesh == "")
                    $mesh = $parts[0];

                if($parts[0] != $mesh)
                    $valid = 0;
            }

            if(!$valid){
                warning "Faces must belong to the same mesh";
                return;
            }

            string $faceStr = stringArrayToString($faces, " ");;
            textFieldButtonGrp -e -text $faceStr clopenerFaceField;
        }

        global proc openWindow()
        {
            if (`window -exists clopenerWindow`)
                deleteUI clopenerWindow;

            window -title "Clopener" -widthHeight 400 500 clopenerWindow;

            columnLayout -adjustableColumn true -rowSpacing 12
                -columnAlign "center"
                -columnAttach "both" 20;

            // Selection mode
            text -label "Selection Mode" -align "center";
            radioButtonGrp
                -numberOfRadioButtons 2
                -labelArray2 "Mesh" "Faces"
                -select 1
                -changeCommand "clopenerToggleMode"
                -columnAlign2 "center" "center"
                clopenerSelectMode;

            // Mesh input
            text -label "Input Mesh" -align "center";
            textFieldButtonGrp
                -label "Mesh"
                -buttonLabel "Use Selected"
                -buttonCommand "string $sel[]=`ls -sl -type transform`; if(size($sel)>0) textFieldButtonGrp -e -text $sel[0] clopenerMeshField;"
                -columnAlign3 "center" "center" "center"
                clopenerMeshField;

            // Face selection
            text -label "Selected Faces" -align "center";
            textFieldButtonGrp 
                -label "Faces"
                -buttonLabel "Use Selected"
                -buttonCommand "clopenerFillFaces() "
                    - columnAlign3 "center" "center" "center"
            clopenerFaceField;

            // Operation
            text -label "Operation" -align "center";
            radioButtonGrp
                -numberOfRadioButtons 2
                -labelArray2 "Opening" "Closing"
                -select 1
                -columnAlign2 "center" "center"
                clopenerOp;

            // Radius
            text -label "Radius (sensitivity) " - align "center";
            floatSliderGrp
                - label "Radius"
                - field true
                - minValue 0.001
                - maxValue 1
                - value 0.2
                 -precision 3
                - columnAlign3 "center" "center" "center"
                clopenerRadius;

            // Target Edge Length
            text - label "Target Edge Length" - align "center";
            floatSliderGrp
                - label "Edge Length"
                - field true
                - minValue 0.0
                - maxValue 1.0
                - value 0.5
                -precision 2
                - columnAlign3 "center" "center" "center"
                clopenerEdgeLength;

            // Iterations
            text - label "Iterations" - align "center";
            intSliderGrp
                - label "Iterations"
                - field true
                - minValue 1
                - maxValue 300
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
            button 
                - label "Apply" 
                - command "processData() ";

            setParent ..;

            clopenerToggleMode();
            showWindow clopenerWindow;
        }

        global proc processData()
        {
            int $iterations = `intSliderGrp -q -value clopenerIterations`;
            float $edgelen = `floatSliderGrp -q -value clopenerEdgeLength`;
            float $radius = `floatSliderGrp -q -value clopenerRadius`;
            int $op = `radioButtonGrp -q -select clopenerOp`;
            int $mode = `radioButtonGrp -q -select clopenerSelectMode`;
            string $mesh = `textFieldButtonGrp -q -text clopenerMeshField`;
            string $faces = `textFieldButtonGrp -q -text clopenerFaceField`;

            print("Iterations: " + $iterations + "\n");
            print("Edge Length: " + $edgelen + "\n");
            print("Radius: " + $radius + "\n");
            print("Operation: " + $op + "\n");
            print("Mode: " + $mode + "\n");
            print("Mesh: " + $mesh + "\n");
            print("Faces: " + $faces + "\n");

            clopenercmd -i $iterations -e $edgelen -r $radius -o $op -d $mode -m $mesh -f $faces;
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