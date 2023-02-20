# FastBCC
This repository contains a simple parallel algorithm for computing closeness centrality in uncertainty graph called FastBCC. 
FaseBCC implements techniques such as mapping vertices to 0-n, ordering neighboring vertices by subscript, bitmap to see whether vertices are accessed, and calculating the distance from one vertex to other vertices in parallel. 
FaseBCC is used as a baseline algorithm to compare with MGMS-BCC algorithm. MGMS-BCC implements techniques proposed in "Closeness Centrality on Uncertain Graphs" (TWEB 23),   
 

## Complie Commend
make


## Run Commend
[executaion file] [datasets file] [is_direccted] [number of source vertices] [file of weight] ([number of thread]) ([sampling number])

The file for parameter [datasets file] is a list of vertex pairs(u v) representing the edges in the graph.
The file for parameter [file of weight] is a list of floating-point numbers corresponding to the probability of each edge in the graph.  
Parameter [is_directed] can only enter 0 and 1, where 0 refers to the undirected graph and 1 refers to the directed graph.
Parameter [file of weight] is the file of edges' probability in graph.
If parameter [number of thread] and [sampling number] are not entered, the default number of 16 threads and 1024 samples are set.

For example:
-----------------------------------------------
./FastBCC ./test.txt 1 5 ./sample_test.txt
./FastBCC ./test.txt 1 5 ./sample_test.txt 1
./FastBCC ./test.txt 1 10 ./sample_test.txt 8 1024
-----------------------------------------------


## Output
[datasets file] [total time] [computation time] [sample time] [number of source vertices] [sample number] [thread number] [algorithm name]

If you want to see the result of closeness centrality, you can uncomment "#define PRINTRESULT" in the code.
