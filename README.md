# Master_Thesis

### Questions:
- What do you think about the organisation of the sections ? Should we do another section for permutation testing ? 
- Persistence Silhouette: interesting to see a lot of samples computed, how many samples ? Though of also comparing median/,ax/min on same gplot for boths groups
- Should we try on another data set ? for example, different cells within the same race (e.g. **reconstruced** dataset)
- Focus on apical dendrite: will add in the appendix results for basal dendrite as well

### General Structure of the report
- Introduction: General introduction on TDA, persistent Homology, why it's useful to biology (especially studies of neurons) and the context of this project
- Theorical Foundation, emphasize on functional summaries 
- Package Overview: introduce TMD and Giotto-TDA
- Application: Mouse vs human data set:
    - general overview of the method confidence band -> permutation test within/in between groups
    - results for different functional summaries: silhouette, landscape, entropy


### Next steps
- *Problem with analysing the axons* $\to$ Told me they used a subset of axons, did only apical dendrite (will send me something to filter)

### Oral presentation

**Tips**

- Give some general context
- could do a whole point on what is persistent diagram, how it works, how we get the silhouette, etc (explain what is the object we work with, also needed definition and/or thm, etc)
- have some graphs to explain the definition over
- have some example where my results could be used (either examples to motivate the project at the beginning, or applications we could use it for at the end) (could for example speak about classification of different type of cell)


**To look into**
- Persistent cohomology enriches point clouds (that are poorly topological)
- The more far a point is from the diagonal, the longer it lives (persistent diagram)

**For the report**
- Have a general introduction on TDA and Homology
- Could also try other functional summaries than persistent silhouette (and maybe justify why one is better than the other) (landscape and entropy)
- Could see which part of the neurons gives the best results for the permutation tests
- permutation test inside the groups itself

