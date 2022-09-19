module PlateLayout

using 
    CSV,
    DataFrames, 
    Statistics, 
    StatsBase, 
    Random,
    Multisets
export #FactorAssignGA
    FactorAssignGA,
    UpdateDesign
    
include("FactorAssignGA.jl")


end # module
