using JuMP, Gurobi, NLPModels, NLPModelsJuMP

probability = [0.4 , 0.6]
rhs = zeros(2,4)
I=2

c2 = [-24  -28 ; 
         -28 -32]
A2 = [6 10 ; 8 5 ; 1 0 ; 0 1]
h = [-60 0 ; 0 -80 ; 0 0 ; 0 0]
b2 = [0 0 500 100 ; 0 0 300 300] 



for i in 1:I
    sub_for_para = Model(Gurobi.Optimizer)
    @variable(sub_for_para , x[1:2])
    @variable(sub_for_para , y[1:2])
    @objective(sub_for_para , Min, c2[i , :]' * y)
    @constraint(sub_for_para,  A2 * y + h * x .≤ b2[i , :])
    nlp=MathOptNLPModel(sub_for_para)
    x = zeros(nlp.meta.nvar)
    rhs[i , :] = nlp.meta.ucon
end
rhs
