Rădulescu Matei

Using matching square algorithm to draw conteour curves parallel implementation description.

In order to parallelize the algorithm the program uses a struct called arg_images to pass the required arguments
to the function used ran by each of the threads implemnted using pthreads C library.

Memory is allocated once containing the image read from the memory, the resized image, the countour map and the grid
needed for the altgorithm, as well as a barrier which make all the threads waits after building each of these as each
subsequent requires something from the previous one. The arguments are passed in a bidimensional vector which all have
pointers pointing to these memories, which are used as a shared space so the threads modify the same image. The 
structure also contains the threads id which is used to devide the steps of the loops used in the algorithm between them.

In order to check if a resize is needed, a flag is sent also in the argument structure which is set to 1 only when
a downwards rescale is needed, if it is not this step can be skipped entirely.

In the function ran by the threads a barrier is used to wait for all of the threads to reach the following points:
1) after the contour map is build 2)after the image is rescaled if it is necessary  3) after the grid is built
This is done to ensure that these operations are fully completed before moving to the next step of the algorithm.

Thank you for your time reading about the implementation of this program.