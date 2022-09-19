# PlateLayout.jl
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/jensenlab/PlateLayout/blob/main/LICENSE)
# Contents 
[Description](#description) \
[Instructions to Use PlateLayout](#instructions-to-use-platelayout) 


# Description 
Automatically generate plate layouts that minimize the number of manual pipetting operations for experimental designs. Currently supports 96 well plates.  

# Instructions to Use PlateLayout 
 Requires installation of  [Julia](https://julialang.org/downloads/). Once Julia is installed. Install PlateLayout by navigating to package mode:  

```julia 
    add https://github.com/jensenlab/PlateLayout
```

# FactorAssignGA 
FactorAssignGA takes an experimental design and uses a genetic algorithm to place experiments on a set of plates. Users must specify which factors are included in the design by supplying an array of the relevant column names. The user also specifies the "type" of factor that each factor is. Factors can be "auto","col","wafer","plate" level. A col factor indicates that all experiments in a column have the same level for that factor. The same holds for wafer factors (an 8x4 wafer for a breakaway 96-well plate) and plate factors. Auto factors are automatically dispensed onto the plate, therefore, their placement does not mater for manual pipetting. This function returns a plate or an array of plates indicating the well location for each experiment. 


```julia
    FactorAssignGA(design,factors,types;num_generations=50,popsize=30)
``` 

### Arguments 
- `design`: An experimental design dataframe. Usually generated in R using functions in the RSM package 
- `factors`: An array of the factor names in the design dataframe that should be included
- `types`: An array of factor types for each factor. Must be one of: "auto", "col", "wafer", "plate".  
### Keyword Arugments 
- `num_generations`: Number of generations in the genetic algorithm. Using more generations improves the chances of converging on a local optimum but is more costly.
- `popsize`: Number of candidate solutions in the population. Using more candidates improves solution diversity but is more costly.  

# UpdateDesign
Converts plates back into an array of design dataframes. One for each plate used.  Even if only one plate is used, it returns an array. The designs now include the well locations in the format used for the R package DispenseFormulatrix.  

```julia
    UpdateDesign(design,plate)
```


### Arguments 
- `design`: The original experimental design dataframe
- `plate`: the solution obtained using FactorAssignGA 






