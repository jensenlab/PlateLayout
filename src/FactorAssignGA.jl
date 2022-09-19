using CSV,DataFrames, Statistics, StatsBase, Random,Multisets


function get_levels(design::DataFrame,factor::String)
    return sort(unique(design[:,factor]))
end 


function plates_upper_bound(design::DataFrame,factors::Array{String,1},types::Array{String,1})
    non_auto_types=types[types.!="auto"]
    non_auto_factors=factors[types.!="auto"]
    if length(non_auto_factors)==0
        n_exps=[nrow(design)]
        n_plates=cld.(nrow(design),96)
        return n_plates
    else 
        non_auto_design=[Array(design[i,non_auto_factors]) for i=1:nrow(design)]
        combos=unique(non_auto_design)
        n_exps=Int.(zeros(length(combos)))
        for i in eachindex(combos)
            n_exps[i]=count(x->x==combos[i],non_auto_design)
        end 
        n_cols=cld.(n_exps,8)
        col_design=vcat(fill.(combos,n_cols)...)
        col_design_matrix=Int.(zeros(length(col_design),length(col_design[1])))
        for i=eachindex(col_design)
            for j=eachindex(col_design[1])
                col_design_matrix[i,j]=col_design[i][j]
            end
        end
        col_design=DataFrame(col_design_matrix)
        rename!(col_design,non_auto_factors)
        if !("wafer" in non_auto_types)
        n_full_cols=[sum(n_cols)]
        wp_combos=combos
        wafer_plate_factors=non_auto_factors
        wafer_plate_types=non_auto_types
        else
            wafer_plate_types=non_auto_types[(non_auto_types.=="wafer") .| (non_auto_types.=="plate")]
            wafer_plate_factors=non_auto_factors[(non_auto_types.=="wafer") .| (non_auto_types.=="plate")]
            wafer_plate_design=[Array(col_design[i,wafer_plate_factors]) for i=1:nrow(col_design)]
            wp_combos=unique(wafer_plate_design)
            n_full_cols=Int.(zeros(length(wp_combos)))
            for i in eachindex(wp_combos)
                n_full_cols[i]=count(x->x==wp_combos[i],wafer_plate_design)
            end
        end 
        n_wafers=cld.(n_full_cols,4)
        wafer_design=vcat(fill.(wp_combos,n_wafers)...)
        wafer_design_matrix=Int.(zeros(length(wafer_design),length(wafer_design[1])))
        for i=eachindex(wafer_design)
            for j=eachindex(wafer_design[1])
                wafer_design_matrix[i,j]=wafer_design[i][j]
            end
        end 
        wafer_design=DataFrame(wafer_design_matrix)
        rename!(wafer_design,wafer_plate_factors)
        if !("plate" in wafer_plate_types)
            n_full_wafers=[sum(n_wafers)]
        else 
            plate_types=wafer_plate_types[(wafer_plate_types.=="plate")]
            plate_factors=wafer_plate_factors[(wafer_plate_types.=="plate")]
            plate_design=[Array(wafer_design[i,plate_factors]) for i=1:nrow(wafer_design)]
            p_combos=unique(plate_design)
            n_full_wafers=Int.(zeros(length(p_combos)))
            for i in eachindex(p_combos)
                n_full_wafers[i]=count(x->x==p_combos[i],plate_design)
            end 
        end
        n_plates=cld.(n_full_wafers,3)
        return sum(n_plates)
    end 
end

function typeorder(type::String)
    if type=="plate"
        return 1
    elseif type=="wafer"
        return 2
    elseif type=="col"
        return 3
    end
end


function initialize_solution(factors::Array{String,1},types::Array{String,1}, n_plates::Int64)
    sol_factors::Array{String,1}=[]
    sol_types::Array{String,1}=[]
    for i in eachindex(types)
        if types[i] == "plate"
            sol_factors=vcat(sol_factors,repeat([factors[i]],n_plates))
            sol_types=vcat(sol_types,repeat([types[i]],n_plates))
        elseif types[i]== "wafer"
            sol_factors=vcat(sol_factors,repeat([factors[i]],n_plates*3))
            sol_types=vcat(sol_types,repeat([types[i]],n_plates*3))
        elseif types[i]== "col"
            sol_factors=vcat(sol_factors,repeat([factors[i]],n_plates*12))
            sol_types=vcat(sol_types,repeat([types[i]],n_plates*12))
        end
    end
    #=
    idxs=shuffle(1:length(sol_factors))
    sol_factors=sol_factors[idxs]
    sol_types=sol_types[idxs]
    =#
    return sol_factors,sol_types
end

function action_space(design::DataFrame,factor::String)
    levels=get_levels(design,factor)
    return collect(0:length(levels))
end 

function map_level(level,design,factor)
    levels=get_levels(design,factor)
    mappingdict=Dict()
    for i in eachindex(levels)
        mappingdict[levels[i]]=i
    end
    mapped=mappingdict[level]
    return mapped
end

function make_evaluate(design::DataFrame,sol_factors,factors,types,n_plates)
    non_auto_types=types[types.!="auto"]
    non_auto_factors=factors[types.!="auto"] 
    factor_dict=Dict()
    for i in eachindex(non_auto_factors)
    factor_dict[non_auto_factors[i]]=non_auto_types[i]
    end 
    type_number_dict=Dict("col" => 8,"wafer" => 32, "plate" => 96)
    non_auto_design=design[:,non_auto_factors]
    for factor in non_auto_factors
        for exp in 1:nrow(design)
            non_auto_design[exp,factor]=map_level(non_auto_design[exp,factor],design,factor)
        end
    end
    non_auto_design=[Array(non_auto_design[i,:]) for i=1:nrow(design)]
    function evaluate(solution)
        allowed_exps=[]
        allowed_wells=[]
        for factor in non_auto_factors
            push!(allowed_wells,repeat(solution[sol_factors.==factor],inner=type_number_dict[factor_dict[factor]]))
        end 
        for i=1:96*n_plates
            push!(allowed_exps,[allowed_wells[j][i] for j in eachindex(allowed_wells)])
        end 
        x=Multiset(non_auto_design)
        y=Multiset(allowed_exps)
        not_plated=length(x-y)
        return nrow(design)-not_plated,nrow(design)
    end
    return evaluate
end


function make_ops_reward(sol_types)
    max_ops=0
    for i in eachindex(sol_types)
        if sol_types[i]=="col"
            max_ops+=1
        elseif sol_types[i]=="wafer"
            max_ops+=4
        elseif sol_types[i]=="plate"
            max_ops+=12
        end 
    end 
    function ops_reward(solution)
        ops=0
        for i in eachindex(sol_types)
            if sol_types[i]=="col"
                ops+=1*solution[i]!=0
            elseif sol_types[i]=="wafer"
                ops+=4*solution[i]!=0
            elseif sol_types[i]=="plate"
                ops+=12*solution[i]!=0
            end 
        end 
        return max_ops-ops
    end 
    return ops_reward
end 


function build_plate(solution, design, sol_factors,sol_types,n_plates)
    n_spots=96*n_plates
    level_plate=reshape([Int.(zeros(length(unique(sol_factors)))) for i =1:n_spots],8,12,n_plates)
    uniq_idx=findfirst.(isequal.(unique(sol_factors)),(sol_factors,))
    typedict=Dict()
    for i in uniq_idx
        typedict[sol_factors[i]]=sol_types[i]
    end
    non_auto_design=design[:,unique(sol_factors)]
    for factor in unique(sol_factors)
        for exp in 1:nrow(design)
            non_auto_design[exp,factor]=map_level(non_auto_design[exp,factor],design,factor)
        end
    end
    non_auto_design=[Array(non_auto_design[i,:]) for i=1:nrow(design)]
    for  i in eachindex(unique(sol_factors))
        factor=unique(sol_factors)[i]
        type=typedict[factor]
        sol_subset=solution[findall(x->x==factor,sol_factors)]
        println(sol_subset)
        if type =="col"
            level_input=repeat(sol_subset,inner=8)
        elseif type == "wafer"
            level_input=repeat(sol_subset,inner=32)
        elseif type == "plate"
            level_input=repeat(sol_subset, inner=96)
        end
        for j in eachindex(level_input)
            level_plate[j][i]=level_input[j]
        end
    end 
    
    plate=Int.(zeros(8,12,n_plates))
    spots=1:n_spots
    n_missing=0
    for i = 1:nrow(design)
        exp=non_auto_design[i]
        for well in spots
            if exp== level_plate[well]
                plate[well]=i
                spots=setdiff(spots,well)
                break
            end
        end 
        if !(exp in level_plate)
            n_missing +=1
        end 
    end
    println("$n_missing missing experiments")
    return plate
end 


function fitness(solution;evaluate=evaluate,ops_reward=ops_reward)
    n_plated,n_exps=evaluate(solution)
    if n_plated <n_exps
        return n_plated
    else
        return n_plated + ops_reward(solution)
    end 
end 


function make_mutate(actions)
    function mutate(solution;p_mutation=0.01,kwargs...)
        for i in eachindex(solution)
            if rand() < p_mutation
                solution[i]=sample(actions[i]) 
            end
        end 
        return solution
    end
    return mutate
end 

function crossover_independent(A,B;p_cross=0.5,kwargs...)
    C=deepcopy(A)
    for i in eachindex(C)
        if rand()>p_cross
            C[i]=A[i]
        else
            C[i]=B[i]
        end
    end
    return C
end

function softmax(population;k=1,kwargs...)
    fitnesses=fitness.(population;kwargs...)
    nums=exp.(k.*fitnesses)
    probs=nums/sum(nums)
    cumul_prob=zeros(length(population))
    for i in eachindex(population)
        cumul_prob=sum(probs[1:i])
    end
    randnum=rand()
    cutoff=minimum(findall(y-> y>=randnum,cumul_prob));
    return population[cutoff]
end 

function create_initial_population_random(actions,popsize=30;kwargs...)
    population=[Int.(zeros(length(actions))) for i =1:popsize]
    for i=1:popsize
        for j in eachindex(actions)
            population[i][j]=sample(actions[j])
        end
    end 
    return population
end 
"""
    FactorAssignGA(design,factors,types;num_generations=50,popsize=30)
Automatically place experiments on a plate to minimize manual pipetting

# Arguments
- `design`: the design dataframe
- `factors`: An array listing the relevant factors in the design
- `types`: An array listing the type of each factor. Must be "auto", "col","wafer", or "plate".
# Keyword Arguments
- `num_generations`: The number of generations the GA should run. 
- `popsize`: The GA population size.

"""
function FactorAssignGA(design::DataFrame,factors::Array{String,1},types::Array{String,1};num_generations::Int64=50,popsize::Int64=30,selection=softmax,crossover=crossover_independent,kwargs...)

    ##Initialize the solution format and initial population
    n_plates=plates_upper_bound(design,factors,types)
    sol_factors,sol_types=initialize_solution(factors,types,n_plates)
    actions=action_space.((design,),sol_factors)
    population=create_initial_population_random(actions;kwargs...)


    ##Define fitness, mutate, evaluate functions
    evaluate=make_evaluate(design,sol_factors,factors,types,n_plates)
    ops_reward=make_ops_reward(sol_types)
    mutate=make_mutate(actions)

    ## Initalize GA algorithm
    best=Int.(zeros(length(actions)))
    best_of_gen=[]
    for gen in 1:num_generations
        population=mutate.(population;kwargs...)
        ranks=sortperm(population,by=x->fitness(x;evaluate=evaluate,ops_reward=ops_reward),rev=true)
        population=population[ranks]
        population=population[1:popsize]
        best_cand=deepcopy(population[1])
        if fitness(best_cand;evaluate=evaluate,ops_reward=ops_reward)>fitness(best;evaluate=evaluate,ops_reward=ops_reward)
            best=deepcopy(best_cand)
        end 
        push!(best_of_gen,best_cand)
        new_population=[]
        for i =1:popsize
            ParentA=selection(population;evaluate=evaluate,ops_reward=ops_reward,kwargs...)
            ParentB=selection(population;evaluate=evaluate,ops_reward=ops_reward,kwargs...)
            Child=crossover(ParentA,ParentB;kwargs...)
            push!(new_population,Child)
        end 
        population=deepcopy(new_population)
    end 

    plate=build_plate(best,design,sol_factors,sol_types,n_plates)
    
    return plate
end 
"""
    UpdateDesign(design,plate)
Add the well locations for each experiment to the design dataframe

# Arguments
- `design`: the design dataframe
- `plate`: the array containing the well locations for each experiment.

"""
function UpdateDesign(design,plate)
    updated=[]
    row,col,n_plates=size(plate)
    alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for i in 1:n_plates
        rows=[]
        cols=[]
        expts=sort(setdiff(unique(plate[:,:,i]),[0]))
        for expt in expts
            idxs=findall(x->x==expt,plate[:,:,i])[1]
            push!(rows,alphabet[idxs[1]])
            push!(cols,idxs[2])
        end
        update=deepcopy(design)[expts,:]
        insert!(update,ncol(design)+1,rows,:Row)
        insert!(update,ncol(design)+2,cols,:Col)
        push!(updated,update)
    end 
    return updated
end

#design = CSV.read("test_design.csv",DataFrame)

#factors=["A","B","C","D"]
    
#types=["auto","col","auto","wafer"] #auto, col, wafer, plate



#plate=FactorAssignGA(design,factors,types)

#new_design=UpdateDesign(design,plate)