# metamoRph

A framework (or "guardrails") for projecting new data onto an existing PCA space. This is a two step process where the user runs our `run_pca` function (which wraps `prcomp` and provides some sensible defaults and enhanced outputs. After `run_pca` the user can then use `metamoRph` to project new data onto the existing PCA. 

We (will) provide a small set of "reference" PCA embeddings of various ocular and body RNA-seq datasets for quick use.
