## Opening and Closing Surfaces Maya Plug-In

Based on paper [Opening and Closing Surfaces SIGGRAPH Asia 2020 SILVIA SELLÁN, JACOB KESTEN, ANG YAN SHENG, ALEC JACOBSON](https://www.dgp.toronto.edu/projects/opening-and-closing-surfaces/)

Compile and build on Visual Studio 2022

To clone:

`git clone --recurse-submodules git@github.com:beckyfeng08/clopener_maya_plugin.git`

Checkout to libigl 2.4.0
`cd clopener_maya_plugin; cd libigl; git checkout v2.4.0`

## How to build the plugin file 

Open Visual Studio 2022, and open the project via its .sln file. Under the Solution Explorer, right click on "helloMaya", navigate to "Properties".
Under Configuration Properties > General > Configuration Type, set to "Dynamic Library (.dll)". Click "Apply".
<img width="1171" height="801" alt="Screenshot 2026-04-02 200557" src="https://github.com/user-attachments/assets/d3099f78-efd0-408c-93e3-732bac1a63e5" />


Under C/C++ > General > Additional Include Directories, 
 make sure you add the eigen-3.4.1 directory, libigl\include directory, and your Autodesk Maya include directories in. Click "Apply":
 
 <img width="671" height="861" alt="Screenshot 2026-04-02 200316" src="https://github.com/user-attachments/assets/5952f9a4-86f0-42c5-b0e6-a2ae0b1877a8" />

Make sure that under C/C++ > Preprocessor > Preprocessor Definitions is set to "NT_PLUGIN". Click "Apply".
<img width="1194" height="658" alt="Screenshot 2026-04-02 200624" src="https://github.com/user-attachments/assets/691868f0-9209-4a24-b4f9-25e1afe032be" />

Under Linker > General, make sure you have the Output File secton set to "$(OutDir)$(ProjectName).mll" and under Additional Library Directories, set to "C:\Program Files\Autodesk\Maya2024\lib". Click "Apply"

<img width="1192" height="785" alt="Screenshot 2026-04-02 201117" src="https://github.com/user-attachments/assets/114c5d34-bd11-45b0-b7dc-1f1e3a3bb0e2" />

Under Linker > Input > Additional Dependencies, be sure that you list the following:

<img width="2066" height="801" alt="Screenshot 2026-04-02 201247" src="https://github.com/user-attachments/assets/55e31df2-3903-4865-902c-0aa9de94582e" />

To build the project and create a .dll file, right click on helloMaya under the Solution Explorer and click "Rebuild"
