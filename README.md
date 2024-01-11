# syntheic_radar_simulation

Satellite imaging is becoming a major asset in effectively seeing how statistics such as vegetation cover
are evolving over time. This is becoming increasingly important given the impacts of climate change.
Unfortunately, satellite images are often masked by cloud cover. This paper looks into classifying
areas masked by clouds, as vegetation or no vegetation. We do this through data collected from other
satellite images of the same area that are not affected by cloud cover. Prior to modelling we created
a new regional variable, scaled down the dataset via re-projection and re-sampling and imputed any
missing values. We subsequently used clustering, support vector machines, logistic regression and
neural networks to attempt to classify vegetated areas. We evaluated the models by recreating images
on test datasets as well as evaluating a 0-1 and brier loss function for synthesised cloud covered areas.
Despite processing limitations the report concludes that our neural network does the best at classifying
vegetated areas, however, only once a similar image has been used in the training dataset. The paper
recommends using logistic regression prior to this.
