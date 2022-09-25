using JuMP, GLPK
using NLPModels, NLPModelsJuMP
using Printf

function print_iteration(k, args...)
    f(x) = Printf.@sprintf("%12.4e", x)
    println(lpad(k, 9), " ", join(f.(args), " "))
    return
end

##########################
# parameters main
c1 = [100 ; 150]
A1 = [1 1 ; -1 0 ; 0 -1]
b1 = [120 ;  -40 ; -20]

##########################
# parameters sub
c2 = [-24 -28 ; 
      -28 -32]
A2 = [6 10 ; 8 5 ; 1 0 ; 0 1] #y
A3 = [-60 0 ; 0 -80 ; 0 0 ; 0 0 ] # x
b2 = [0  0  500  100;
      0  0  300  300]

#############################
p = [0.4 0.6]


####################################
# defining the main model but don't solve it
####################################
main = Model(GLPK.Optimizer)
@variable(main, x[1:2])
@variable(main, -1000000 ≤ θ)
@objective(main, Min, c1' * x)
@constraint(main, A1 * x .≤ b1);

#################
# I define the sub problem inside functions

function sub_dual(x)
    λ = zeros(2 , 4) # 2 realization and 4 constraints
    for i in 1:2 # 2 realization
        sub = Model(GLPK.Optimizer)
        @variable(sub, 0 ≤ y[1:2])
        @objective(sub, Min, c2[i , :]' * y)
        @constraint(sub, A2 * y + A3 * x .≤ b2[i , :])
        optimize!(sub)
        all_con = all_constraints(sub, AffExpr, MOI.LessThan{Float64})
        λ[i , :] = dual.(all_con)
        λ = round.(λ , digits =2)
    end
    return λ
end

###############################
# I define a function to give me the h and T
# I don't need to solve this model 

function sub_coef()
    T = [] # coefficient of x
    h = zeros(2,4)  # RHS
    for i in 1:2
        sub = Model(GLPK.Optimizer)
        @variable(sub, x[1:2])
        @variable(sub, 0 ≤ y[1:2])
        @objective(sub, Min, c2[i , :]' *y)
        @constraint(sub, A2 * y + A3 * x .≤ b2[i , :])
        nlp = MathOptNLPModel(sub)
        s = zeros(nlp.meta.nvar)
        T = jac(nlp , s)[ : , 1:2]
        h[ i , :] = nlp.meta.ucon
    end 
    return (T , h)
end

sub_coef()

function initialize()
    for i in 1:10
        optimize!(main)
        xᵏ = value.(x)
        θᵏ = value(θ)
        λᵏ = sub_dual(xᵏ)
        Tᵏ = sub_coef()[1]
        hᵏ = sub_coef()[2]
        eᵏ = sum(p[i] * λᵏ[i,j]' * hᵏ[i,j] for i in 1:2 for j in 1:4)
        Eᵏ = sum(p[i] * λᵏ[i, :]' * Tᵏ  for i in 1:2)
        wᵛ = eᵏ - Eᵏ * xᵏ 
        print_iteration(i,    θᵏ,    wᵛ)
        if θᵏ ≥ wᵛ -10
            println("**********************************")
            println("****   we are at optimality*****")
            break
        end
        cut = @constraint(main, θ ≥ eᵏ - Eᵏ * x )
        println("add the $(cut)" )
    end
end

initialize()
