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

d = [
    0.7  0.8  0.8  0.6  0.9  1.0  1.0  0.9  0.0  1.0
    0.2  0.2  1.0  0.9  0.7  0.7  0.4  1.0  0.3  0.7
    0.5  0.6  0.5  0.1  0.4  0.6  0.4  0.0  0.5  0.6
    0.5  1.4  0.3  1.0  0.8  0.8  0.2  0.8  0.6  0.6
    0.8  0.2  1.3  0.6  0.3  0.4  0.6  0.8  0.4  0.8
    0.3  0.4  0.9  0.2  0.6  1.0  0.4  0.3  0.4  0.6
]


Random.seed!(1234)

#μ = [1.2, 1.8, 1.6]
#Σ = [0.5 0.0 0.0; 0.0 0.5 0.0; 0.0 0.0 0.5]
#dξ = reshape(rand(MvNormal(μ,Σ), Ξ ), (J,Tˢ,Ξ) )
#dξ = round.(dξ , digits = 2)

dξ = [1.2 , 1.8 , 3.2]

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
     sum(χ[h,j] for j in 1:J ) ≤ sum(ψ[h,j] for j in 1:J) * J  );


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


################################################################
# sub_dual -> λ
###########################################
# when we want to get the dual the x (first stage variables are considered as fixed)
# when we want to get the coefficents of x, x should be variable
# therefore, we need to define two models one for dual and one for coefficient
###########################################
###############################################################

function sub_dual(ψ , χ , αᵏ)
    sub = Model(GLPK.Optimizer) 
    @variable(sub,  α[1:H, 1:J, 1:Tˢ , 1:Ξ])
    @variable(sub, 0 ≤ γ[1:J , 1:Tˢ , 1:Ξ])
    @objective(sub, Min,  sum(zᶜ[j] * γ[j,t,ξ] for j in 1:J for t in 1:Tˢ for ξ in 1:Ξ ) )
    for h in 1:H
            for j in 1:J
                for t in 1:Tˢ
                    for ξ in 1:Ξ
                        if (αᵏ[h , j , 1 , 1] == 0 )
                            con = @constraint(sub, α[h,j,t,ξ] == 0)
                        else
                            con = @constraint(sub, α[h,j,t,ξ] == 1)
                        end
                    end
                end
            end
        end

    con1_S2 = @constraint(sub, demand_met[ j in 1:J , t in 1:Tˢ , ξ in 1:Ξ],
         dξ[j,t,ξ] ≤ sum(α[h,j,t,ξ] * (ψ[h,j] + χ[h,j]) for h in 1:H ) + γ[j,t,ξ]   )
    con2_s2 = @constraint(sub, permanent_allocability[h in 1:H , j in 1:J , t in 1:Tˢ , ξ in 1:Ξ],
                α[h,j,t,ξ] ≤ ψ[h,j] + 10 * χ[h,j] );
    con3_s2 = @constraint(sub, no_more_than_one_station[h in 1:H , t in 1:Tˢ , ξ in 1:Ξ],
               sum( α[h,j,t,ξ] for j in 1:J) ≤ 1 );
    optimize!(sub)

    con_equal = all_constraints(sub, AffExpr, MOI.EqualTo{Float64})
    con_less = all_constraints(sub, AffExpr, MOI.LessThan{Float64})
    λ1 = dual.(con_equal)
    λ2 = dual.(con_less)
    λ = append!(λ1 , λ2)

    no_con_equal = length(con_equal)
    no_con_less = length(con_less)
    no_all_con = no_con_equal + no_con_less;
    #@show no_con_equal
    #@show no_con_less
    #@show no_all_con;
    
    return λ
end



#############################################
# function for coefficients of  𝜓  and  𝜒
#############################################

function sub_coeff(αᵏ)  
    sub = Model(GLPK.Optimizer) 
    @variable(sub, ψ[1:H , 1:J], Bin )
    @variable(sub, 0 ≤ χ[1:H , 1:J])
    @variable(sub,  α[1:H, 1:J , 1:Tˢ, 1:Ξ])
    @variable(sub, 0 ≤ γ[1:J , 1:Tˢ , 1:Ξ])
     @objective(sub, Min,  sum(zᶜ[j] * γ[j,t,ξ] for j in 1:J for t in 1:Tˢ for ξ in 1:Ξ ) )
    for h in 1:H
            for j in 1:J
                for t in 1:Tˢ
                    for ξ in 1:Ξ
                        if αᵏ[h , j , 1,  1] == 0
                            con = @constraint(sub, α[h,j,t,ξ] == 0)
                        else
                            con = @constraint(sub, α[h,j,t,ξ] == 1)
                        end
                    end
                end
            end
        end

    con1_S2 = @NLconstraint(sub, demand_met[ j in 1:J , t in 1:Tˢ , ξ in 1:Ξ],
         dξ[j,t,ξ] ≤ sum(α[h,j,t,ξ] * (ψ[h,j] + χ[h,j]) for h in 1:H ) + γ[j,t,ξ]   )
    con2_s2 = @constraint(sub, permanent_allocability[h in 1:H , j in 1:J , t in 1:Tˢ , ξ in 1:Ξ],
                α[h,j,t,ξ] ≤ ψ[h,j] + 10 * χ[h,j] );
    con3_s2 = @constraint(sub, no_more_than_one_station[h in 1:H , t in 1:Tˢ , ξ in 1:Ξ],
               sum( α[h,j,t,ξ] for j in 1:J) ≤ 1 );

    vr = all_variables(sub)
    vr_index = [vr[i].index.value for i in 1:length(vr)]
    df = DataFrame(variable = vr , index = vr_index); 
    #@show df
    
    nlp = MathOptNLPModel(sub)
    q = zeros(nlp.meta.nvar)
    jac(nlp, q)
    A1 = jac(nlp, q)[ : , 1:24]
    return A1
end


println("        k   lower bound  upper bound   gap")
function initiate()
    for k in 1:10
        optimize!(main)
        ψᵏ = value.(ψ)
        χᵏ = value.(χ)
        x = all_variables(main)[1 : (H *J*2) ]
        lb = objective_value(main)
        sub_solve(ψᵏ , χᵏ)
        αᵏ = sub_solve(ψᵏ , χᵏ)[1]
        γᵏ = sub_solve(ψᵏ , χᵏ)[2]
        oᵏ = sub_solve(ψᵏ , χᵏ)[3]
        λᵏ = sub_dual(ψᵏ , χᵏ , αᵏ)
        A1 = sub_coeff(αᵏ)
        ub = sum( zᵖ[h] * ψᵏ[h,j] * (Tᵈ + Tˢ) + zᵐ[h,j] * χᵏ[h,j] for h in 1:H for j in 1:J) + 
             sum(zᶜ[j] * γᵏ[j,t,ξ] for j in 1:J for t in 1:Tˢ for ξ in 1:Ξ ) 
        gap = (ub - lb ) / ub
        print_iteration(k , lb , ub, gap)
        if gap < 0.000001
            println("************************************")
            println("***     we are at optimality   ***  ")
            break
        end
        #############################
        cut = @constraint(main, θ ≥ oᵏ - (λᵏ)' * A1 * (x .- value.(x) ) )
        println("we are adding the cut $(cut)")
    end
end

initiate()
