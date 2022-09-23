using JuMP
using NLPModelsJuMP, NLPModels
using GLPK
using Distributions
using Random
using DataFrames
using CSV
using Printf

function print_iteration(k, args...)
    f(x) = Printf.@sprintf("%12.4e", x)
    println(lpad(k, 9), " ", join(f.(args), " "))
    return
end

###############################################################
# PARAMETERS
###############################################################
H = 4
J = 3
Tᵈ = 1
Tˢ = 1
Ξ = 1

s = [1.0  0.0  0.0;  0.0  1.0  0.0; 0.0  0.0  1.0; 1.0  0.0  0.0]

ψ = [ 1 0 0 ; 0 1 0 ; 0 0 1 ; 1 0 0 ];
χ = zeros(H,J);
d = [
    0.7  0.8  0.8  0.6  0.9  1.0  1.0  0.9  0.0  1.0
    0.2  0.2  1.0  0.9  0.7  0.7  0.4  1.0  0.3  0.7
    0.5  0.6  0.5  0.1  0.4  0.6  0.4  0.0  0.5  0.6
    0.5  1.4  0.3  1.0  0.8  0.8  0.2  0.8  0.6  0.6
    0.8  0.2  1.3  0.6  0.3  0.4  0.6  0.8  0.4  0.8
    0.3  0.4  0.9  0.2  0.6  1.0  0.4  0.3  0.4  0.6
]


Random.seed!(1234)

μ = [1.2, 0.8, 0.6]

Σ = [0.5 0.0 0.0;
     0.0 0.5 0.0;
     0.0 0.0 0.5
    ]

dξ = reshape(rand(MvNormal(μ,Σ), Ξ ), (J,Tˢ,Ξ) )
dξ = round.(dξ , digits = 2)

χᵐᵃˣ = fill(1, H ,J)


zᵖ = [7; 9; 11;7 ] * 1000
zᶜ = zᵖ *2

zᵐ = fill(0, H, J)

for h=1
        zᵐ[h , 2] = (1/2) * zᵖ[h]
        zᵐ[h , 3] = (3/4) * zᵖ[h]
end

for h=2
        zᵐ[h , 3] = (1/2) * zᵖ[h]
        zᵐ[h , 1] = (3/4) * zᵖ[h]
end

for h=3
        zᵐ[h , 1] = (1/2) * zᵖ[h]
        zᵐ[h , 2] = (3/4) * zᵖ[h]
end

for h=4
        zᵐ[h , 2] = (1/2) * zᵖ[h]
        zᵐ[h , 3] = (3/4) * zᵖ[h]
end

;

#####################################
#####################################
#####################################
skill = Model(GLPK.Optimizer)

@variable(skill, ψ[1:H , 1:J], Bin)
@variable(skill, 0 ≤ χ[1:H , 1:J] ≤ 1)
@variable(skill, α[1:H , 1:J , 1:Tᵈ], Bin)
@variable(skill, z1[1:H , 1:J , 1:Tᵈ], Bin)
@variable(skill, z2[1:H , 1:J , 1:Tᵈ])

####################
# stage 2 variables
####################

@variable(skill, 0 ≤ γξ[1:J , 1:Tˢ , 1:Ξ])
@variable(skill, αξ[1:H , 1:J , 1:Tˢ , 1:Ξ], Bin)
@variable(skill, z3[1:H , 1:J , 1:Tˢ , 1:Ξ], Bin)
@variable(skill, z4[1:H , 1:J , 1:Tˢ , 1:Ξ])

@objective(skill, Min, sum( zᵖ[h] * ψ[h,j] * (Tᵈ + Tˢ) for h in 1:H for j in 1:J ) +
                            sum(zᵐ[h,j] * χ[h,j] for h in 1:H for j in 1:J ) +
                         (1/Ξ) * sum(zᶜ[j] * γξ[j,t,ξ] for j in 1:J for t in 1:Tˢ for ξ in 1:Ξ)
        )

#####################
# stage 1 constraints
#####################
@constraint(skill, linear1[h = 1:H , j = 1:J , t =1:Tᵈ], z1[h,j,t] ≤ α[h,j,t] )
@constraint(skill, linear2[h = 1:H , j = 1:J , t =1:Tᵈ], z1[h,j,t] ≤ ψ[h,j] )
@constraint(skill, linear3[h = 1:H , j = 1:J , t =1:Tᵈ], z1[h,j,t] ≥ α[h,j,t] + ψ[h,j] -1 )

@constraint(skill, linear4[h = 1:H , j = 1:J , t =1:Tᵈ], z2[h,j,t] ≤ α[h,j,t] )
@constraint(skill, linear5[h = 1:H , j = 1:J , t =1:Tᵈ], z2[h,j,t] ≤ χ[h,j] )
@constraint(skill, linear6[h = 1:H , j = 1:J , t =1:Tᵈ], z2[h,j,t] ≥ α[h,j,t] + χ[h,j] -1 )

@constraint(skill, demand_met[ j in 1:J , t in 1:Tᵈ ], sum(z1[h,j,t] + z2[h,j,t] for h in 1:H ) ≥ d[j,t] )
@constraint(skill, permanent_pool[h in 1:H , j in 1:J], ψ[h,j] ≤ s[h,j] );
@constraint(skill, allocation_recruitment_training[h in 1:H , j in 1:J , t in 1:Tᵈ], α[h,j,t] ≤ ψ[h,j] + 10 * χ[h,j] )
@constraint(skill, no_more_than_one_station[ h in 1:H , t in 1:Tᵈ ], sum(α[h,j,t] for j in 1:J) ≤ 1 );
@constraint(skill, single_skilling_at_recruitment[h in 1:H , j in 1:J], χ[h,j] ≤ ( 1 - ψ[h,j] ) * χᵐᵃˣ[h,j] );
@constraint(skill, training_recruited[ h in 1:H ], sum(χ[h,j] for j in 1:J ) ≤ sum(ψ[h,j] for j in 1:J) * J  );

###################################################
# Stage 2
###################################################


@constraint(skill, linear7[h = 1:H , j = 1:J , t =1:Tˢ , ξ = 1:Ξ], z3[h,j,t,ξ] ≤ αξ[h,j,t,ξ] )
@constraint(skill, linear8[h = 1:H , j = 1:J , t =1:Tˢ , ξ = 1:Ξ], z3[h,j,t,ξ] ≤ ψ[h,j] )
@constraint(skill, linear9[h = 1:H , j = 1:J , t =1:Tˢ , ξ = 1:Ξ], z3[h,j,t,ξ] ≥ αξ[h,j,t,ξ] + ψ[h,j] -1 )

@constraint(skill, linear10[h = 1:H , j = 1:J , t =1:Tˢ , ξ = 1:Ξ], z4[h,j,t,ξ] ≤ αξ[h,j,t,ξ] )
@constraint(skill, linear11[h = 1:H , j = 1:J , t =1:Tˢ , ξ = 1:Ξ], z4[h,j,t,ξ] ≤ χ[h,j] )
@constraint(skill, linear12[h = 1:H , j = 1:J , t =1:Tˢ , ξ = 1:Ξ], z4[h,j,t,ξ] ≥ αξ[h,j,t,ξ] + χ[h,j] -1 )

@constraint(skill, stage2_demand_met[j = 1:J , t =1:Tˢ , ξ = 1:Ξ], sum(z3[h,j,t,ξ] for h in 1:H ) + sum(z4[h,j,t,ξ] for h in 1:H) + γξ[j,t,ξ] ≥ dξ[j,t,ξ] )
@constraint(skill, stage2_permanent_allocability[h in 1:H , j in 1:J , t in 1:Tˢ , ξ = 1:Ξ], αξ[h,j,t,ξ] ≤ ψ[h,j] + 10 * χ[h,j] )
@constraint(skill, stage2_no_more_than_one_station[ h in 1:H , t in 1:Tˢ , ξ = 1:Ξ], sum(αξ[h,j,t,ξ] for j in 1:J) ≤ 1 )

#####################################
# optimize
#####################################
optimize!(skill);

println("objective value : ", objective_value(skill))



df_ψ = DataFrame(value.(ψ), :auto)
@show df_ψ 

df_χ = DataFrame(value.(χ), :auto)
@show df_χ

println("total number of recruited workforce : ", sum(value.(ψ)) )

println("total amount of secondary skills : ", sum(round.(value.(χ), digits = 2) ) )
