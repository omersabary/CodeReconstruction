# CodeReconstruction
Reconstruction Algorithms for DNA Storage Systems
In the trace reconstruction problem a length-n string x yields a collection of noisy copies, called traces, y1, …, yt where each yi is independently obtained from x by passing through a deletion channel, which deletes every symbol with some fixed probability. The main goal under this paradigm is to determine the required minimum number of i.i.d traces in order to reconstruct x with high probability. The trace reconstruction problem can be extended to the model where each trace is a result of x passing through a deletion-insertion-substitution channel, which introduces also insertions and substitutions. Motivated by the storage channel of DNA, this work is focused on another variation of the trace reconstruction problem, which is referred by the DNA reconstruction problem. A DNA reconstruction algorithm is a mapping Embedded Image which receives t traces y1, …, yt as an input and produces Embedded Image, an estimation of x. The goal in the DNA reconstruction problem is to minimize the edit distance Embedded Image between the original string and the algorithm’s estimation. For the deletion channel case, the problem is referred by the deletion DNA reconstruction problem and the goal is to minimize the Levenshtein distance Embedded Image.

In this work, we present several new algorithms for these reconstruction problems. Our algorithms look globally on the entire sequence of the traces and use dynamic programming algorithms, which are used for the shortest common supersequence and the longest common subsequence problems, in order to decode the original sequence. Our algorithms do not require any limitations on the input and the number of traces, and more than that, they perform well even for error probabilities as high as 0.27. The algorithms have been tested on simulated data as well as on data from previous DNA experiments and are shown to outperform all previous algorithms.
