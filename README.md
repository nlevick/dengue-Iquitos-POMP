# dengue-Iquitos-POMP
Model Code for "Mechanistic plug-and-play methods for understanding the impact of control and climate on seasonal dengue dynamics"
----------
### UNM Master's Thesis in Mathematics -- Nathan Levick ###
### Advisor -- Dr. Helen Wearing ###

----------


Model Descriptions
------
The repository possess five models with the following `model_name`:

1. **EIP** -- Temperature-dependent extrinsic incubation period (EIP).
2. **Control** -- Control measures.
3. **ControlG** -- Control measures + altered probability of mosquito exposure.
4. **ControlGH** -- Control measures + altered probability of mosquito & human exposure.
5. **Control2b0** -- Control measures + altered probability of mosquito & human exposure + two average value transmission parameters.


***
Order of Execution
------------
In each Models/`model_name` folder, there is are six files:

 - **1_`model_name`\_run.r** -- The first half of the MIF process
 - **2_`model_name`\_run.r** -- The second half of the MIF process
 - **3_`model_name`\_model.r** -- Plot and evaluate model states
 - **4_`model_name`\_profile.r** -- Plot model profiles
 - **`model_name`_rprocess.txt** -- The associated `rprocess` for the model.
 - **`model_name`.rds** -- The resulting R data information from 2\_`model_name`_run.r

The first two models should be modified to perform new MIF data rds files.
***
General Tips
------------
Change the 15 in the `ddply` to alter the number of top particles selected on each profiled value. 
```{r}
bakedRDS %<>%
    # remove failed particles
    subset(nfail.max==0) %>%
    mutate(S_0=exp(signif(log(S_0),5))) %>%
    # change number of best particles for each profile value
    ddply(~S_0,subset,rank(-loglik)<=15)
```
For instance, if there are 15 values for each free parameter and 15 values for the profiled parameter, there is 225 total particles. To see the top 150 particles set
```{r}
    ddply(~S_0,subset,rank(-loglik)<=10)
```  
***
