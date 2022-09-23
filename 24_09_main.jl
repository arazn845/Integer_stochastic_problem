

#########################################
# main model
#########################################

main = Model(GLPK.Optimizer)
@variable(main, ψ[1:H , 1:J], Bin)
@variable(main, 0 ≤ χ[1:H , 1:J])
@variable(main, α[1:H , 1:J , 1:Tᵈ], Bin)
@variable(main, -1000000 ≤ θ)

@variable(main, z1[1:H , 1:J , 1:Tᵈ], Bin)
@variable(main, z2[1:H , 1:J , 1:Tᵈ])

@objective(main, Min, sum( zᵖ[h] * ψ[h,j] * (Tᵈ + Tˢ) + zᵐ[h,j] * χ[h,j] 
                           for h in 1:H for j in 1:J) + θ);


#con1_s1 = @NLconstraint(main, demand_met[j in 1:J , t in 1:Tᵈ], 
#                          sum(  α[h,j,t] *(ψ[h,j] + χ[h,j]) for h in 1:H ) ≥ d[j,t] );

#################
# LINEARIZING
#################
@constraint(main, linear1[h = 1:H , j = 1:J , t =1:Tᵈ], z1[h,j,t] ≤ α[h,j,t] )
@constraint(main, linear2[h = 1:H , j = 1:J , t =1:Tᵈ], z1[h,j,t] ≤ ψ[h,j] )
@constraint(main, linear3[h = 1:H , j = 1:J , t =1:Tᵈ], z1[h,j,t] ≥ α[h,j,t] + ψ[h,j] -1 )

@constraint(main, linear4[h = 1:H , j = 1:J , t =1:Tᵈ], z2[h,j,t] ≤ α[h,j,t] )
@constraint(main, linear5[h = 1:H , j = 1:J , t =1:Tᵈ], z2[h,j,t] ≤ χ[h,j] )
@constraint(main, linear6[h = 1:H , j = 1:J , t =1:Tᵈ], z2[h,j,t] ≥ α[h,j,t] + χ[h,j] -1 )

con1_s1 = @constraint(main, demand_met[ j in 1:J , t in 1:Tᵈ ], 
    sum(z1[h,j,t] + z2[h,j,t] for h in 1:H ) ≥ d[j,t] );

con2_s1 = @constraint(main, permanent_pool[h in 1:H , j in 1:J], ψ[h,j] ≤ s[h,j]);

con3_s1 = @constraint(main, allocation_recruitment_training[h in 1:H , j in 1:J , t in 1:Tᵈ],
                      α[h,j,t] ≤ ψ[h,j] + 10 * χ[h,j] );

con4_s1 = @constraint(main, no_more_than_one_station[h in 1:H , t in 1:Tᵈ],
              sum(α[h,j,t] for j in 1:J) ≤ 1  )

con4_s1 = @constraint(main, single_skilling_at_recruitment[h in 1:H , j in 1:J], 
    χ[h,j] ≤ ( 1 - ψ[h,j] ) * χᵐᵃˣ[h,j] );

con5_s1 = @constraint(main, training_recruited[ h in 1:H ], 
  
