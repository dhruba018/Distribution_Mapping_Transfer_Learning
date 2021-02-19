### Distribution Mapping based Transfer Learning

This is the repo for distribution mapping based transfer learning.  
* `source_set`, `target_set` => estimate distribution mapping => `DM`
* `target_set` => Use `DM` to move to source space => `target_set_mapped`
* `source_set` => Train predictive model => `Model`
* `target_set_mapped`=> predict response in the source space using `Model` => `target_response_pred_source`
* `target_response_pred_source` => Use `DM` to map back to target space => `target_response_pred`

* ---- Always under work ----*
For the CRAN released R package, take a look at the [DMTL](https://github.com/dhruba018/DMTL) repo. 
