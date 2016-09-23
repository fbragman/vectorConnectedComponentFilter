# vectorConnectedComponentFilter
Vector-based connected component analysis for image processing

This project provides code to perform vector-based connected component labelling of binary images based on a user-defined heuristic. It is based on the itkVectorConnectedComponentImageFilter. The original connected component labelling algorithm labels binary volumes, based on voxel connectivity. The vector-based connected component analysis provides the user with more flexibilty on the decisions for connectivity. In addition to voxel connectivity, the user can define an extra set of functions. For example, the difference of voxel intensities or the dot product of eigenvectors can be used.

Syntax: [result,eqTable] = vectorConnectedComponentFilter(img,mask)


