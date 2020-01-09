# Parallel-Computing
## Project of Parallel Computing for ENSAI MSSD

In this project, I implemented MLE of the logistic regression by Newton's Method.
Further, CV and stepwise selection are implemented.
Finally, the goal was to optimize it by vectorization and parallelization of the code.

Report.pdf is a report of the project, which I suggest as entry point.

Here is a file description:

* aux_fcts.R are differents stages of implementation of likelihood,gradient,hessian and Newton-Rhapson algorithm.
* basics.R builds on aux_fcts.R: Cross-validation, Model selection by stepwise algorithm.
* datagen.R generates data 
* run.R and run.html are different runs of the algorithm and their profiling with profviz.


