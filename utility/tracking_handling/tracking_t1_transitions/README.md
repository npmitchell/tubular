# tracking-t1-pipeline
detect_t1_transitions.m reads tracking results from "labeled_groundTruth" and produces a list of pairs of cells that are neighbors in at least one time frame, then query these pairs to get rid of pairs that are always neighbors since such pairs do not involve the creation and elimination of cell-cell junctions. The result is stored in "not_always_neighbors.mat" in the form of a matrix of binaries that describes the state of each pair of cells (1 being neighbors and 0 otherwise). The first index is the pair id, and the second index is the time point.

screenPotentialT1.m screens the previously identified potential T1 events by showing to the user the movies of colored pairs of cells. We assume that when rows in "not_always_neighbors" matrix matches pattern "0011" or "0011," a junction is created or eliminated. To remove any artifacts caused by the code or edge cases, we examine these videos and save the results to "confirmed_t1.mat"

"identify_t1_overlap.m" removes the T1 events that are double counted from the previous stage since a single T1 event could be detected twice from merging and splitting and lets the user inspect the end results.
