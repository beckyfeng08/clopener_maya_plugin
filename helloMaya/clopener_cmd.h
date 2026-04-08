#pragma once
//#include "clopener.h"

#include <maya/MPxCommand.h>
#include <maya/MStatus.h>    // <--- needed for MStatus
#include <maya/MArgList.h>   // <--- needed for MArgList
#include <string>

class clopenercmd : public MPxCommand
{
public:
    clopenercmd();
    virtual ~clopenercmd();

    static void* creator() { return new clopenercmd(); }

    // Explicit return type now recognized
    MStatus doIt(const MArgList& args) override;
};