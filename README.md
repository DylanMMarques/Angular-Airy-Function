# Angular Airy function

The scripts presented in this paper were used to compute the results from the paper:

If you publish scientific results based on this scripts, please consider citing the paper.

The results are based on an optical model of Fabry-Pérot (FP) etalons illuminated with a focused beam. All scripts have been implemented in [Julia](https://julialang.org/)

The 4 scripts available do:
* angularAiryFunction.jl - implementation of the equations shown in the paper;
* modelCompairisonExperimental.jl - compute the data shown in fig. 2;
* intuitiveUnderstanding.jl - compute the data shown in fig. 3;
* modelCompairison.jl - show that the angular Airy function matches a multilayer model (data not shown in the paper);

To compute the results, we use [Jolab.jl](https://github.com/DylanMMarques/Jolab.jl) which (indirectly) implements the angular Airy function. The same results can be computed based on the scripts available in angularAiryFunction.jl
