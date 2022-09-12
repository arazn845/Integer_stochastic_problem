using JuMP, Gurobi

#####################################
# parameters
#####################################

ψ = [
    1 0 0 ;
    0 1 0 ;
    0 0 1 ;
    1 0 0 
]

χ = [
    0 0 0;
    0 0 0;
    0 0 0;
    0 0 0
]

H = 4
J = 3
Tᵈ = 1
Tˢ = 1

d = [
    0.3  1.2   
    0.2  0.2
    1.0  0.0
]

zᶜ = [7;
      7;
      9;
      11
] * 1000

#####################################
# sub
#####################################
sub = Model(Gurobi.Optimizer)
@variable(sub, 0 ≤ γ[1:J , 1:Tˢ])
@variable(sub, α[1:H , 1:J , 1:Tˢ], Bin)

@objective(sub, Min,  sum(zᶜ[j] * γ[j,t] for j in 1:J for t in 1:Tˢ  ) )

con1 = @constraint(sub, stage2_demand_met[j = 1:J , t =1:Tˢ], sum( α[h,j,t]*ψ[h,j]  for h in 1:H ) + sum(α[h,j,t]*χ[h,j] for h in 1:H) + γ[j,t] ≥ d[j,t] )
con2 = @constraint(sub, stage2_permanent_allocability[h in 1:H , j in 1:J , t in 1:Tˢ], ψ[h,j] + 10 * χ[h,j] ≥ α[h,j,t]  )
con3 = @constraint(sub, stage2_no_more_than_one_station[ h in 1:H , t in 1:Tˢ], 1 ≥ sum(α[h,j,t] for j in 1:J)  );

optimize!(sub);

s = value.(α)

#####################################
# sub_for_dual
#####################################

sub_for_dual = Model(Gurobi.Optimizer) 
@variable(sub_for_dual, 0 ≤ α[1:H,1:J,1:Tˢ])
@variable(sub_for_dual, 0 ≤ γ[1:J , 1:Tˢ])
@objective(sub_for_dual, Min,  sum(zᶜ[j] * γ[j,t] for j in 1:J for t in 1:Tˢ  ) )
for h in 1:H
        for j in 1:J
            for t in 1:Tˢ
                if s[h,j,t] == 0
                    con = @constraint(sub_for_dual, - α[h,j,t] ≥ 0)
                else
                    con = @constraint(sub_for_dual, α[h,j,t] ≥ 1)
                    con = @constraint(sub_for_dual, - α[h,j,t] ≥ -1)
                end
            end
        end
    end

con = @constraint(sub_for_dual, stage2_demand_met[j = 1:J , t =1:Tˢ], sum( α[h,j,t]*ψ[h,j]  for h in 1:H ) + sum(α[h,j,t]*χ[h,j] for h in 1:H) + γ[j,t] ≥ d[j,t] )
con = @constraint(sub_for_dual, stage2_permanent_allocability[h in 1:H , j in 1:J , t in 1:Tˢ], ψ[h,j] + 10 * χ[h,j] ≥ α[h,j,t] )
con = @constraint(sub_for_dual, stage2_no_more_than_one_station[ h in 1:H , t in 1:Tˢ], 1 ≥ sum(α[h,j,t] for j in 1:J)  )
optimize!(sub_for_dual)

JuMP.has_duals(sub_for_dual)

#####################################
# dual_vector
#####################################
all_consts = all_constraints(sub_for_dual , AffExpr , MOI.GreaterThan{Float64})
dual.(all_consts)
