using Random, Statistics, Plots, StatsBase, DataFrames, CSV


plate=[
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 1 1 1 0 0 0 0 0 0 0 0;
    0 1 0 1 1 0 0 1 0 0 1 0;
    0 1 0 0 0 0 0 1 0 1 1 0;
    0 1 0 0 0 1 0 1 0 1 0 0;
    0 1 0 0 0 0 0 1 0 1 0 0;
    0 1 0 1 0 0 0 1 0 1 0 0;
    0 1 1 1 1 1 1 1 1 1 0 0 
]

function capture_chain(plate,row,col,direction)
    R,C=size(plate)
    chain=[];
    if plate[row,col]!=1
        return(chain)
    end
    if direction=="r"
        cand=plate[row,:]
        chain=[(row,col)]
        lower=col-1
        flagl=true
        upper=col+1
        flagu=true
        while lower>0 && flagl
            if cand[lower]==1
                push!(chain,(row,lower))
                lower-=1
                flagl=true
            elseif cand[lower]!=1
                flagl=false
            end 
        end 
        while upper<=C && flagu
            if cand[upper]==1
                push!(chain,(row,upper))
                upper+=1
                flagu=true
            elseif cand[upper]!=1
                flagu=false
            end
        end 
    elseif direction=="c"
        cand=plate[:,col]
        chain=[(row,col)]
        lower=row-1
        flagl=true
        upper=row+1
        flagu=true
        while lower>0 && flagl
            if cand[lower]==1
                push!(chain,(lower,col))
                lower-=1
                flagl=true
            elseif cand[lower]!=1
                flagl=false
            end 
        end 
        while upper<=R && flagu
            if cand[upper]==1
                push!(chain,(upper,col))
                upper+=1
                flagu=true
            elseif cand[upper]!=1
                flagu=false
            end
        end 
    end 
    if length(chain)>8
        chain=rank_distance((row,col),chain)[1:8]
    end
    return(chain)
end 

function rank_distance(cell::Tuple{Int64,Int64},chain::Array{Tuple{Int64,Int64},1})
    dists=p_norm.(chain,(cell,))
    ranked=chain[sortperm(dists)]
    return ranked
end 

function p_norm(x,y;p=1,kwargs...)
    x=collect(x)
    y=collect(y)
    return sum(abs.(x-y).^p)^(1/p)
end

function get_wells(plate)
    R,C=size(plate)
    wells=[]
    for i in 1:R
        for j in 1:C
            if plate[i,j]==1
                push!(wells,(i,j))
            end
        end
    end 
    return wells
end

function bit_flipper(orig_plate,wells)
    plate=deepcopy(orig_plate)
    for well in wells
        if plate[well[1],well[2]]==1
            plate[well[1],well[2]]=0
        elseif plate[well[1],well[2]]==0
            plate[well[1],well[2]]=1
        end 
    end 
    return plate
end 

function assign_solution_wells(solution,wells,step)
    for well in wells
        solution[well[1],well[2]]=step
    end 
    return solution
end 


function find_valid_plan(plate)
    wells=shuffle(get_wells(plate))
    rm_wells=[]
    step=0
    while length(wells)>0
        step+=1
        direction=sample(["r","c"])
        rm_chain=capture_chain(plate,wells[1][1],wells[1][2],direction)
        rm_chain=setdiff(rm_chain,rm_wells)
        rm_wells=vcat(rm_wells,rm_chain)
        wells=setdiff(wells,rm_chain)
    end 
    return step
end 

    
function simulate_random(plate;nsims=100)
    steps=zeros(nsims)
    for i=1:nsims
        steps[i]=find_valid_plan(plate)
    end 
    return mean(steps)
end

function plate_operations_rollout(plate;simulate=simulate_random,kwargs...)
    R,C=size(plate)
    working_plate=deepcopy(plate)
    solution=Int.(zeros(R,C))
    wells=shuffle(get_wells(plate))
    step=0
    best_steps=R*C
    while length(wells)>0
        step+=1
        
        best_wells=[]
        for well in wells
            for direction in ["r","c"]
                sim_wells=capture_chain(working_plate,well[1],well[2],direction)
                sim_plate=bit_flipper(working_plate,sim_wells)
                steps=simulate(sim_plate;kwargs...)
                if steps < best_steps
                    best_wells=sim_wells
                    best_steps=steps
                end 
            end 
        end
        working_plate=bit_flipper(working_plate,best_wells)
        wells=get_wells(working_plate)
        solution=assign_solution_wells(solution,best_wells,step)
    end 
    return solution 
end


function random_plate(chainlens)
    R=8
    C=12
    plate=Int.(zeros(R,C))
    for chain in chainlens
        if chain<9
            direction=sample(["r","c"])
            if direction =="r"
                col=sample(1:C)
                A=sample(1:R-chain+1)
                B=A+chain-1
                plate[A:B,col].=1
            elseif direction=="c"
                row=sample(1:R)
                A=sample(1:C-chain+1)
                B=A+chain-1
                plate[row,A:B].=1
            end
        else
            row=sample(1:R)
            A=sample(1:C-chain+1)
            B=A+chain-1
            plate[row,A:B].=1
        end 
    end
    return plate
end 





#=

function benchmark(outfile)
    nexps=1000
    sims=[5,10,20,50,100]
    experiment=[]
    operations=[]
    time=[]
    nsims=[]
    for exp in 1:nexps
        plate=random_plate([9,6,5,3,3,2,2,1,1,1])
        for sim in sims
            solution=@timed plate_operations_rollout(plate;nsims=sim)
            ops=maximum(solution.value)
            push!(experiment,exp)
            push!(operations,ops)
            push!(time,solution.time)
            push!(nsims,sim)
        end 
    end 
    df=DataFrame(experiment=experiment,operations=operations,time=time,nsims=nsims)
    CSV.write(outfile,df)
    tottime=sum(time)
    println("Total Time: $tottime")
end

benchmark("plate_operations_rollout_benchmark.csv")
=#
