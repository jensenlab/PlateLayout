using JuMP, Plots, Gurobi

function Plate_Operations(plate)
    (r,c)=size(plate)
    n=Int(sum(plate))
    model=Model(Gurobi.Optimizer)
    @variable(model,y[1:r,1:c,1:n],Bin);
    @variable(model,chains[1:n])

    @constraint(model, sum(y)==n);

    for k=1:n
        for i =1:r
            for j=1:c
                @constraint(model,chains[k]>=y[i,j,k])
            end 
        end 
        @constraint(model,chains[k]<=sum(y[:,:,k]))
        @constraint(model,chains[k]<=1)
    end 

    for i=1:r
        for j=1:c
            @constraint(model,sum(y[i,j,:])==plate[i,j])
        end 
    end 

    for i=1:r
        for j=1:c
            if plate[i,j]==1
                for k=1:n
                    if i==1 && j==1
                        @constraint(model,y[i+1,j,k]+y[i,j+1,k]<=1)
                    elseif i==1 && j>1 && j<c
                        @constraint(model, y[i+1,j,k]+y[i,j+1,k]<=1)
                        @constraint(model, y[i+1,j,k]+y[i,j-1,k]<=1)
                    elseif i==1 && j==c
                        @constraint(model, y[i+1,j,k]+y[i,j-1,k]<=1)
                    elseif i>1 && i<r && j==1
                        @constraint(model,y[i+1,j,k]+y[i,j+1,k]<=1)
                        @constraint(model,y[i-1,j,k]+y[i,j+1,k]<=1)
                    elseif i==r && j==1
                        @constraint(model,y[i-1,j,k]+y[i,j+1,k]<=1)
                    elseif i==r && j>1 && j<c 
                        @constraint(model,y[i-1,j,k]+y[i,j+1,k]<=1)
                        @constraint(model,y[i-1,j,k]+y[i,j-1,k]<=1)
                    elseif i==r && j==c
                        @constraint(model,y[i-1,j,k]+y[i,j-1,k]<=1)
                    elseif i>1 && i<r && j==c 
                        @constraint(model,y[i+1,j,k]+y[i,j-1,k]<=1)
                        @constraint(model,y[i-1,j,k]+y[i,j-1,k]<=1)
                    else
                        @constraint(model,y[i+1,j,k]+y[i,j+1,k]<=1)
                        @constraint(model,y[i-1,j,k]+y[i,j+1,k]<=1)
                        @constraint(model,y[i+1,j,k]+y[i,j-1,k]<=1)
                        @constraint(model,y[i-1,j,k]+y[i,j-1,k]<=1)
                    end 
                end 
            end 
        end
    end
    for k=1:n 
        @constraint(model, sum(y[:,:,k])<=8)
    end 
    @objective(model,Min,sum(chains));
    optimize!(model)
    return round.(JuMP.value.(y))
end


plate=[
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 1 1 1 0 0 0 0 0 0 0 0;
    0 1 0 1 1 0 0 1 0 0 1 0;
    0 1 0 0 0 0 0 1 0 1 1 0;
    0 1 0 0 0 1 0 1 0 1 0 0;
    0 1 0 0 0 0 0 1 0 1 0 0;
    0 1 0 1 0 0 0 1 0 1 0 0;
    0 1 1 1 1 1 1 1 1 0 0 0 
]

test=Plate_Operations(plate)                

                    






    

