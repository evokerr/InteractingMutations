# InteractingMutations
This program (written in R), takes files (csv) of mutations detected through whole genome sequencing over time on a pair of co-cultured microbes and determines whether a mutation in one species (the "initiator") may be setting up a second mutation in the other species (the "responder"). The code includes functions to filter/find candidate pairs of mutations, as well as give a simple read-out of the fraction of co-cultures in which the responder mutation did not precede the initiator mutation. Also, we have written a permutation test function to give a statistical test on whether the ordering of the mutations is unusually extreme. Preliminary comments on each function are available in the code.
