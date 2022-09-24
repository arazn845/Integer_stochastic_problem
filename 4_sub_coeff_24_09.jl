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
