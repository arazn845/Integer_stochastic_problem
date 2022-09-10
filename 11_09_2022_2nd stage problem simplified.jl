using JuMP, GLPK

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
Tˢ = 2

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

sub = Model(GLPK.Optimizer)
@variable(sub, 0 ≤ γ[1:J , 1:Tˢ])
@variable(sub, α[1:H , 1:J , 1:Tˢ], Bin)

@objective(sub, Min,  sum(zᶜ[j] * γ[j,t] for j in 1:J for t in 1:Tˢ  ) )

con1 = @constraint(sub, stage2_demand_met[j = 1:J , t =1:Tˢ], sum( α[h,j,t]*ψ[h,j]  for h in 1:H ) + sum(α[h,j,t]*χ[h,j] for h in 1:H) + γ[j,t] ≥ d[j,t] )
con2 = @constraint(sub, stage2_permanent_allocability[h in 1:H , j in 1:J , t in 1:Tˢ], α[h,j,t] ≤ ψ[h,j] + 10 * χ[h,j] )
con3 = @constraint(sub, stage2_no_more_than_one_station[ h in 1:H , t in 1:Tˢ], sum(α[h,j,t] for j in 1:J) ≤ 1 )
;
