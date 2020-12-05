### Distribution Mapping based Transfer Learning

This is the repo for distribution matching based transfer learning.  
    * source_set, target_set => estimate distribution mapping
    * target_set => Use DM to move to source space
    * source_set => Train predictive model => Mapped target_set => predict response in the source space
    * Predicted target response => Use DM to map back to target space

Currently under work...