using CSV,DataFrames, Statistics, StatsBase, Random
include("Plate_OperationsRL.jl")



function get_levels(design::DataFrame,factor::String)
    return sort(unique(design[:,factor]))
end 

function wafer_priority(design,factor,levels)
    level_counts=Int.(zeros(length(levels)))
    for i in eachindex(level_counts)
        level_counts[i]=count(j->j==levels[i],design[:,factor])
    end 
    idxs=sortperm(level_counts;rev=true)
    types=["manual" for i in eachindex(level_counts)]
    for i=1:3
        types[idxs[i]]="wafer$i"
    end
    return types
end
    
function level_types(design,factors,types)
    levels=get_levels.((design,),factors)
    level_types=[]
    for i in  eachindex(factors)
        if types[i]=="auto"
            push!(level_types,["auto" for j in 1:length(levels[i])])
        elseif types[i]=="manual"
            push!(level_types,["manual" for j in 1:length(levels[i])])
        elseif types[i]=="wafer"
            push!(level_types,wafer_priority(design,factors[i],levels[i]))
        end 
    end
    return level_types
end

function find_valid_plate(current_plate,design, factors, levels,type_levels)
    R,C=size(current_plate)
    set_wells=[]
    for i=1:R
        for j=1:C
            if current_plate[i,j]!=0
                push!(set_wells,(i,j))
            end 
        end 
    end 
    plate=deepcopy(current_plate)
    wells=[(i,j)  for j=1:C for i =1:R]
    wells=setdiff(wells,set_wells)
    runs_accounted_for=setdiff(unique(current_plate),[0])
    total_runs=collect(1:nrow(design))
    runs_needed=setdiff(total_runs,runs_accounted_for)
    for i in runs_needed
        well=(0,0)
        run=design[i,factors]
        for j in eachindex(factors)
            level=run[factors[j]]
            #println("$level")
            idx=findfirst(x->x==level,levels[j])
            #println("$idx")
            l_type=type_levels[j][idx]
            #println("$l_type")
            if l_type=="wafer1"
                well=sample(wells[findall(x->x[2]<=4,wells)])
                #println("$well")
                #println("$wells")
                break
            elseif l_type=="wafer2"
                well=sample(wells[findall(x->x[2]>=5 && x[2]<=8,wells)])
                break
            elseif l_type=="wafer3"
                well=sample(wells[findall(x->x[2]>=9,wells)])
                break
            else
                well=sample(wells)
            end 
        end
        plate[well[1],well[2]]=i
        wells=setdiff(wells,[well])
    end 
    return plate
end  



function evaluate_plate(plate, design,factors,levels,type_levels;kwargs...)
    n_operations=0
    R,C=size(plate)
    plates=[]
    for i in eachindex(factors)
        for j in eachindex(type_levels[i])
            if type_levels[i][j]=="manual"
                factor=factors[i]
                level=levels[i][j]
                level_runs=[]
                for run in 1:nrow(design)
                    if design[run,factor]==level
                        push!(level_runs,run)
                    end 
                end 
                level_runs_plate=Int.(zeros(R,C))
                for run in level_runs
                    level_runs_plate[plate.==run].=1
                end
                solution=plate_operations_rollout(level_runs_plate;kwargs...)
                n_operations+=maximum(solution)
                push!(plates,level_runs_plate)
            end 
        end 
    end 
    return n_operations
end 

function available_wells(plate)
    R,C=size(plate)
    wells=[(i,j) for j in 1:C for i in 1:R ]
    used=[]
    for i=1:R
        for j=1:C
            if plate[i,j]!=0
                push!(used,(i,j))
            end 
        end 
    end 
    wells=setdiff(wells,used)
    return wells
end 

function action_space(plate,run,factors,levels,type_levels)
    wells=available_wells(plate)
    for i in eachindex(factors)
        level=run[factors[i]]
        idx=findfirst(x->x==level,levels[i])
        l_type=type_levels[i][idx]
        if l_type=="wafer1"
            wells=wells[findall(x->x[2]<=4,wells)]
            #println("$well")
            #println("$wells")
            break
        elseif l_type=="wafer2"
            wells=wells[findall(x->x[2]>=5 && x[2]<=8,wells)]
            break
        elseif l_type=="wafer3"
            wells=wells[findall(x->x[2]>=9,wells)]
            break
        end
    end 
    return wells
end 

function simulate_experiments_random(current_plate,design, factors, levels,type_levels;nsims_plate=100,kwargs...)
    steps=zeros(nsims_plate)
    for i=1:nsims_plate
        steps[i]=evaluate_plate(find_valid_plate(current_plate,design, factors, levels,type_levels),design,factors,levels,type_levels)
    end 
    return mean(steps)
end 

function plate_layout_rollout(design,factors,types;R=8,C=12,simulate=simulate_experiments_random,kwargs...)
    levels=get_levels.((design,),factors)
    type_levels=level_types(design,factors,types)
    nexp=nrow(design)
    plate=Int.(zeros(R,C))
    for exp in 1:nexp
        best_well=(0,0)
        best_operations=length(factors)*R*C
        actions=shuffle(action_space(plate,design[exp,factors],factors,levels,type_levels))
        for well in actions
            working_plate=deepcopy(plate)
            working_plate[well[1],well[2]]=exp
            operations=simulate(working_plate,design,factors,levels,type_levels;kwargs...)
            if operations < best_operations
                best_well=well
                best_operations=operations
            end 
        end 
        plate[best_well[1],best_well[2]]=exp
    end 
    return plate
end 

design = CSV.read("test_design.csv",DataFrame)

factors=["A","B","C","D"]

types=["auto","auto","wafer","manual"] #Auto, Wafer, or Manual

levels=get_levels.((design,),factors)
type_levels=level_types(design,factors,types)
current_plate=[
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 1 2 3 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
]
plate=find_valid_plate(current_plate,design,factors,levels,type_levels)
            

@timed ops=evaluate_plate(plate,design,factors,levels,type_levels)


wells=action_space(current_plate,design[1,factors],factors,levels,type_levels)        

plate=plate_layout_rollout(design,factors,types;nsims=10,nsims_plate=10)
