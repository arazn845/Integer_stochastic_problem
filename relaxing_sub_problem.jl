sub_for_dual = Model(GLPK.Optimizer) 
@variable(sub_for_dual, α[1:H,1:J,1:Tˢ])
@variable(sub_for_dual, 0 ≤ γ[1:J , 1:Tˢ])
@objective(sub_for_dual, Min,  sum(zᶜ[j] * γ[j,t] for j in 1:J for t in 1:Tˢ  ) )
for h in 1:H
        for j in 1:J
            for t in 1:Tˢ
                if s[h,j,t] == 0
                    con = @constraint(sub_for_dual, α[h,j,t] == 0)
                else
                    con = @constraint(sub_for_dual, α[h,j,t] == 1)
                end
            end
        end
    end

con = @constraint(sub_for_dual, stage2_demand_met[j = 1:J , t =1:Tˢ], sum( α[h,j,t]*ψ[h,j]  for h in 1:H ) + sum(α[h,j,t]*χ[h,j] for h in 1:H) + γ[j,t] ≥ d[j,t] )
con = @constraint(sub_for_dual, stage2_permanent_allocability[h in 1:H , j in 1:J , t in 1:Tˢ], α[h,j,t] ≤ ψ[h,j] + 10 * χ[h,j] )
con = @constraint(sub_for_dual, stage2_no_more_than_one_station[ h in 1:H , t in 1:Tˢ], sum(α[h,j,t] for j in 1:J) ≤ 1 )
optimize!(sub_for_dual)

JuMP.has_duals(sub_for_dual)
