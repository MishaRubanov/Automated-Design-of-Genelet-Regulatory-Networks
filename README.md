# Automated Design of Genelet Regulatory Networks
## Automated design using multiobjective optimization (MOO) through PyMOO for genelet network concentrations and topology

### Code descriptions:
* **GeneralGeneletModel**: a modified version of the general genelet model from Schaffer et al: https://assets.researchsquare.com/files/rs-247740/v1_covered.pdf?c=1631857099 
* **Autoamplifiermodel**: a set of functions for running multi-objective optimization of the general genelet model to discover, in this case, amplifiers with high gain and minimal leak
* **Cascading_amplifier**: an example script that runs an optimizer to discover a set of concentrations that leads to a genelet network with a 'good' amplifier
