################################################
# sub solve   outout-> α γ o
################################################

function sub_solve(ψ , χ)
    sub = Model(GLPK.Optimizer)
    @variable(sub, α[1:H , 1:J , 1:Tˢ , 1:Ξ], Bin )
    @variable(sub, 0 ≤ γ[1:J , 1:Tˢ , 1:Ξ]);
    @objective(sub, Min,  (1/Ξ) * sum(zᶜ[j] * γ[j,t,ξ] for j in 1:J for t in 1:Tˢ for ξ in 1:Ξ) )
    con1_S2 = @constraint(sub, demand_met[ j in 1:J , t in 1:Tˢ , ξ in 1:Ξ],
     sum(α[h,j,t,ξ] * (ψ[h,j] + χ[h,j]) for h in 1:H ) + γ[j,t,ξ] ≥ dξ[j,t,ξ]  );
    con2_s2 = @constraint(sub, permanent_allocability[h in 1:H , j in 1:J , t in 1:Tˢ , ξ in 1:Ξ],
            α[h,j,t,ξ] ≤ ψ[h,j] + 10 * χ[h,j] );
    con3_s2 = @constraint(sub, no_more_than_one_station[h in 1:H , t in 1:Tˢ , ξ in 1:Ξ],
           sum( α[h,j,t,ξ] for j in 1:J) ≤ 1 );
    optimize!(sub)
    α = value.(α)
    γ = value.(γ)
    o = objective_value(sub)
    α = round.(α , digits = 2)
    γ = round.(γ , digits = 2)
    o = round.(o , digits = 2)
    return (α , γ , o )
end

