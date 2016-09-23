# vectorConnectedComponentFilter
Vector-based connected component analysis for image processing

This project provides code to perform vector-based connected component labelling of binary images based on a user-defined heuristic. It is based on the itkVectorConnectedComponentImageFilter. The original connected component labelling algorithm labels binary volumes, based on voxel connectivity. The vector-based connected component analysis provides the user with more flexibilty on the decisions for connectivity. In addition to voxel connectivity, the user can define an extra set of functions. For example, the difference of voxel intensities or the dot product of eigenvectors can be used.

Note 1) only dotproduct of eigenvectors is hard-coded, more flexibility for user-defined functions will be included in later     releases

Note 2) C++ mex files will be included in later releases to deal with processing speed in large 3D medical volumes

    Syntax: [result,eqTable] = vectorConnectedComponentFilter(img,mask)
    
    Inputs
    img  - image to apply decision function
    mask - binary volume of img
    
    Outputs
    result  - labelled volume of img/mask
    eqTable - equivalency table


