# Power-Simulations

These syntax files can be ued to simulate data and conduct power analyses.

1. 'Simulation of mediation analyses.R' simulates a) the total effect of a mediator on a continuous outcome given a
binary treatment and b) the indirect effect of that mediator provided varying levels of mediation (quantified as
percentage of total effect explained by mediator).

2. 'Simulation of repeated measures power.R' simulates repeated measures data provided number of level-2 units,
number of level-1 units per level-2 unit, ICC, effect size, and rate of missingness. Note that this syntax can be
used to estimate power to detect time X treatment effects.

4. 'Simulation of repeated measures power for within-person design.R' simulates repeated measures data provided number
of level-2 units, number of level-1 units per level-2 unit, ICC, effect size, and rate of missingness. Note that this
estimate does not consider cross-level interactions, for instance, in the case of a time X treatment effect.
