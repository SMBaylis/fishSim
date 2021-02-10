# fishSim
R tool for simulation of fish population dynamics

Simple tool for fish population dynamics, written so that I have full
flexibility on what I simulate. Handles generating an initial
population, markovian movement between populations, mating,
sex-switching, aging, and flexible mortality patterns. Assigns
parentage at all breeding events to allow reconstruction of kin
relationships.  Also includes archiving and data-management tools.

Mating can be set to have an age at first breeding (and 'flat'
fecundity per-mating at greater ages), fully age-specific fecundity
per mating, or fully age- and sex-specific fecundity per mating.

Mortality can be set to bring a population down to a specific size by
randomly killing individuals, to randomly kill individuals with a set
probability (without setting a target population size, allowing
stochastic population growth), or to randomly kill individuals with a
stock-specific probability, age-specific probability, or
stock-and-age-specific probability.
