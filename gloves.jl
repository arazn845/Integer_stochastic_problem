using JuMP, Gurobi
gloves = Model(Gurobi.Optimizer);
@variable(gloves, 0 ≤ r)
@variable(gloves, 0 ≤ c);
@objective(gloves, Max, 5*r + 8*c)
con1 = @constraint(gloves, cutting, r + 1.5*c ≤ 900)
con2 = @constraint(gloves, finishing, 0.5*r + 0.33*c ≤300)
con3 = @constraint(gloves, packing, 0.125*r+0.25*c ≤100)
optimize!(gloves);

dual.(cutting)
