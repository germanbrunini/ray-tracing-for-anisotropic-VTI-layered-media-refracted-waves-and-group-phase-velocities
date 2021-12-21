


function isotropic_ray(
    e     :: Array{Float64,1},
    v     :: Array{Float64,1},
    n_d   :: Int64           ,
    s     :: Array{Float64,1},
    r     :: Array{Float64,1},
    ftol  :: Float64         ,
    d_out :: Array{Float64,1})

    # #====================================================================
    # #
    # # Esta subrutina calcula el rayo en un medio isotropo, utilizando
    # # la técnica propuesta en
    # #
    # # "A rapid and accurate two-point ray tracing method in
    # # horizontally layered velocity model" de Tian y Chen  #
    # # IN: e     --> espesores de las capas entre la fuente y el receptor.
    # #     v     --> velocidad de la onda, puede ser s o p
    # #     n_int --> cantida de disc. entre s y r
    # #     s     --> coordenadas de la fuente
    # #     r     --> coordenadas del receptor
    # #
    # # OUT: d_out --> coordenadas x de las intersecciones entre
    # #                el rayo y las disc. del modelo
    # #
    # #====================================================================
    #
    # #averiguo cuantas y cuales son las capas con v=v_max, y calculo el
    # #espesor e_max equivalente

    v_max,v_max_ind = findmax(v)       ; # Max vel from mod and index
    e_max           = sum(e[v.==v_max]); # Add all H with v_max. Equivalent H

    # Horizontal distance between source and receiver
    delta = norm(s[1:2]-r[1:2])  ;

    # calculo el q inicial
    a = 0.0;
    b = 0.0;
    q = 0.0;
    for i = 1:n_d                 # eq 9a
        a = a + (v[i]/v_max)*e[i]
        if (v[i] != v_max)
            b = b + ((v[i]/v_max)*e[i])/sqrt(1.0-(v[i]/v_max)^2); # eq 10a
        end
    end
    a     = a/e_max     ;  # eq 9a

    aba_1 = (a*b)/(a-1) ; # eq 12

    if (delta < aba_1);
        q = delta/a;
    elseif (delta > aba_1)
        q = delta - b;
    end
    # hago newton-rhapson para averiguar el parámetro q
    delta_f = F(e,v,v_max,e_max,q) #delta_f inicial
    while (abs(1.0-delta_f/delta) > ftol)
        num     = F(e,v,v_max,e_max,q)-delta;
        den     = F_p(e,v,v_max,e_max,q);
        q       = q - num/den;
        delta_f = F(e,v,v_max,e_max,q)
    end

    # calculo el tiempo de viaje
    t_out = 0.0;
    for i = 1:n_d
        raiz = sqrt(1.0-((q/(v_max*sqrt(e_max*e_max + q*q)))*v[i])^2);
        t_out = t_out + e[i]/(v[i]*raiz);
    end
    # calculo las distancias parciales hasta el receptor
    # que tambien sirven para calcular las coordenadas de la intersección
    # del rayo con las discontinuidades
    for i=1:n_d-1
        d_out[i] = F(e[1:i],v[1:i],v_max,e_max,q);
    end
    return d_out,t_out
end

function F(
    e     :: Array{Float64,1},
    v     :: Array{Float64,1},
    v_max :: Float64         ,
    e_max :: Float64         ,
    q     :: Float64         )

    # !=====================================================================
    # ! Esta funcion calcula la distancia delta en función del modelo de
    # ! velocidades y del parámetro q
    # ! la función devuelve un arreglo con los valores de la distancia
    # ! horizontal recorrida por cada segmento del rayo
    # !=====================================================================
    F = 0.0 ;
    for i = 1:length(e)
        vimax = v[i]/v_max;
        F = F + (vimax)*e[i]*q/sqrt(e_max*e_max + (1.0-(vimax*vimax))*q^2);
    end
    return F;
end

function F_p(
    e     :: Array{Float64,1},
    v     :: Array{Float64,1},
    v_max :: Float64         ,
    e_max :: Float64         ,
    q     :: Float64)
    # !======================================================================
    # !Esta funcion es la derivada de la anterior, mas o menos...
    # !======================================================================
    F_p = 0.0 ;
    for i = 1:length(e)
        vimax = v[i]/v_max;
        F_p = F_p + (vimax)*e[i]/(e_max*e_max + (1.0-(vimax*vimax))*q^2)^(1.5);
    end
    F_p = F_p * e_max * e_max;
    return F_p;
end







# function isotropic_ray_g(
#     h     :: Array{Float64,1},
#     v     :: Array{Float64,1},
#     n_d   :: Int64           ,
#     s     :: Array{Float64,1},
#     r     :: Array{Float64,1},
#     ftol  :: Float64         ,
#     d_out :: Array{Float64,1})
#     # #====================================================================
#     # #
#     # # Esta subrutina calcula el rayo en un medio isotropo, utilizando
#     # # la técnica propuesta en
#     # #
#     # # "A rapid and accurate two-point ray tracing method in
#     # # horizontally layered velocity model" de Tian y Chen  #
#     # # IN: e     --> espesores de las capas entre la fuente y el receptor.
#     # #     v     --> velocidad de la onda, puede ser s o p
#     # #     n_int --> cantida de disc. entre s y r
#     # #     s     --> coordenadas de la fuente
#     # #     r     --> coordenadas del receptor
#     # #
#     # # OUT: d_out --> coordenadas x de las intersecciones entre
#     # #                el rayo y las disc. del modelo
#     # #
#     # #====================================================================
#     #
#     # #averiguo cuantas y cuales son las capas con v=v_max, y calculo el
#     # #espesor e_max equivalente

#     vₘ,vₘ_ind = findmax(v)     ; # Max vel from mod and index
#     hₘ        = sum(h[v.==vₘ]) ; # Add all H with v_max. Equivalent H

#     # Horizontal distance between source and receiver
#     Δ = norm(s[1:2]-r[1:2])  ;

#     # calculo el q inicial
#     a = 0.0; b = 0.0; q = 0.0;

#     ε  = v/vₘ ; # eq.6 auxiliars
#     ε² = ε.*ε ;

#     for i = 1:n_d                 # eq 9a
#         a = a + ε[i]*h[i];
#         if (v[i] != vₘ)
#             b = b + (ε[i]*h[i])/sqrt(1.0-ε²[i]); # eq 10a
#         end
#     end
#     a     = a/hₘ     ;  # eq 9a

#     Δc = (a*b)/(a-1.0) ; # eq 12
#     if (Δ < Δc); # eq 13
#         q = Δ/a;
#     elseif (Δ > Δc)
#         q = Δ - b;
#     end

#     # newton-rhapson para averiguar el parámetro q
#     Δf = F(h,ε,hₘ,q); # Δf inicial
#     while (abs(1.0-Δf/Δ) > ftol)
#         num = F(h,ε,hₘ,q)-Δ       ; # eq. 7
#         @show den = ∇F(h,ε,hₘ,q)  ;
#         q   = q - num/den         ; # eq.8
#         Δf  = F(h,ε,hₘ,q)         ; # evaluate in new "q".
#     end

#     # calculo el tiempo de viaje
#     t_out = 0.0;
#     for i = 1:n_d
#         raiz = sqrt(1.0-((q/(vₘ*sqrt(hₘ*hₘ + q*q)))*v[i])^2);
#         t_out = t_out + h[i]/(v[i]*raiz);
#     end
#     # calculo las distancias parciales hasta el receptor
#     # que tambien sirven para calcular las coordenadas de la intersección
#     # del rayo con las discontinuidades
#     for i=1:n_d-1
#         d_out[i] = F(h[1:i],v[1:i],hₘ,q);
#     end
#     return d_out,t_out
# end

# # !=====================================================================
# # ! Esta funcion calcula la distancia delta en función del modelo de
# # ! velocidades y del parámetro q
# # ! la función devuelve un arreglo con los valores de la distancia
# # ! horizontal recorrida por cada segmento del rayo
# # !=====================================================================

# F(
# h :: Array{Float64,1},
# ε :: Array{Float64,1},
# hₘ:: Float64         ,
# q ::Float64) = q*sum((ε.*h)./(sqrt.(hₘ*hₘ .+  (1.0.-(ε.*ε))*q^2 )));


# # !======================================================================
# # !Esta funcion es la derivada de la anterior, mas o menos...
# # !======================================================================
# ∇F(
# h :: Array{Float64,1},
# ε :: Array{Float64,1},
# hₘ:: Float64         ,
# q ::Float64) = hₘ*hₘ*sum((ε.*h)./( hₘ*hₘ .+ (1.0.-(ε.*ε))*q^2 ).^(1.5));


