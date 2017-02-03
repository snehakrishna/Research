# Identifying a Ranking of Plant Preferences for a Pollinator

This project aims to model a pollinator's interaction with plant species in a meadow. This code base is the final code base used in the undergraduate thesis titled "Identifying a Ranking of Plant Preferences for a Pollinator".

This project was conducted by Sneha Krishna Kumaran working under Prof. Rebecca Hutchinson and Prof. Tom Dietterich at Oregon State University. 

We modeled a pollinator’s interaction with plant species in two ways: 
  1. by using a probabilistic multinomial approach to fit a preference score to each plant and 
  2. explaining our findings from (1) tby using the traits of the flowers themselves. 
Our findings showed that a model with preferences performs better than a model which does not have preferences. While this model shows potential, it does not fully explain the distribution of plant-pollinator interactions. 

To explain the interactions more fully, we incorporated the traits of the plants into the score of the plant. We found that the traits do have some effect on the score of the plant, but this model again does not fully explain the interactions.

Our model learned pollinator preferences by utilizing four years of data collected at the HJ Andrews Experimental Forest.

Code written for the pupose of the thesis can be found in the folder Code_datasets_generateddata
