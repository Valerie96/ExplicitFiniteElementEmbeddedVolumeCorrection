# ExplicitFiniteElementEmbeddedVolumeCorrection
Explicit dynamic finite element code that provides an option to remove the volume redunancy inherit to the embedded element method
Also known as Flagshyp Modified (FM) in some work. This code is an updaded version of FM. 
The main function to run the code is ExplicitVolumeCorrectionEmbeddedFiniteElement.m
Code documentation and changes made from the orginial FLagshyp code base are descibed in the Word document ExplicitVolumeCorrectionEmbeddedFiniteElement-Documentation.docx
Functions for creating input files for this code are found in the InputFileCreation folder. This folder includes functions that read Abaqus input files and add embedded truss elements to extruded geometries (EmbeddedElements_090Truss_multipleparts.m) and to a specific geometry of curved plates (EmbeddCurvedPlate.m). More instructions for those functions are found in the function commments.

Additional information on the creation of this code can be found in my thesis (when it is published, I'll let you know). 
