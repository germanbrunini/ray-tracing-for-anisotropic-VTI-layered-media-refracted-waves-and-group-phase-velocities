function critical_angles(
    rangId :: Array{Int64,1},
    an     :: Array{Float64,2},
    vp     :: Array{Float64,1},
    vs     :: Array{Float64,1})

    # extracted from 
    # Landrø, M., & Tsvankin, I. (2007). 
    # Seismic critical-angle reflectometry: A method to characterize
    # azimuthal anisotropy? GEOPHYSICS, 72(3), D41–D50.
    # doi:10.1190/1.2437145 
    # Equation 1 for any phase, PP, SV, SH. 

    # extract only involved layers and on their correct order.  
    layₙ = length(rangId)
    ϵ,δ,γ  = an[rangId,1],an[rangId,2],an[rangId,3];
    vₚₒ    = vp[rangId]  ;
    vₛₒ    = vs[rangId]  ;

    Tₚₚ = zeros(Float64,layₙ-1);
    Tₛᵥ = copy(Tₚₚ);
    Tₛₕ = copy(Tₚₚ);
    for l = 2:layₙ
        FPP,FSV,FSH = Criticalfunctions(vₚₒ,vₛₒ,ϵ,δ,γ,l);
        # angles need to be in radians
        Tₚₚ[l-1] = rad2deg(bisection_2(FPP,0.0,pi/2.0));
        Tₛₕ[l-1] = rad2deg(bisection_2(FSH,0.0,pi/2.0));
        Tₛᵥ[l-1] = rad2deg(bisection_2(FSV,0.0,pi/2.0));

    end
    return (Tₚₚ,Tₛₕ,Tₛᵥ);
end    



function Criticalfunctions(
    vₚₒ :: Array{Float64,1},
    vₛₒ :: Array{Float64,1},
    ϵ :: Array{Float64,1},
    δ :: Array{Float64,1},
    γ :: Array{Float64,1},
    l :: Int64)

    numₚₚ(θ) = vₚₚₐ(vₚₒ[l-1],vₛₒ[l-1],θ    ,ϵ[l-1],δ[l-1],γ[l-1]);
    denₚₚ    = vₚₚₐ(vₚₒ[l  ],vₛₒ[l  ],π/2.0,ϵ[l  ],δ[l  ],γ[l  ]);
    FuPP(θ)  = numₚₚ(θ)/denₚₚ;
    FPP(θ)   = sin(θ)-FuPP(θ); # function ready for bisection

    numₛᵥ(θ) = vₛᵥₐ(vₚₒ[l-1],vₛₒ[l-1],θ    ,ϵ[l-1],δ[l-1],γ[l-1]);
    denₛᵥ    = vₛᵥₐ(vₚₒ[l  ],vₛₒ[l  ],π/2.0,ϵ[l  ],δ[l  ],γ[l  ]);
    FuSV(θ)  = numₛᵥ(θ)/denₛᵥ;
    FSV(θ)   = sin(θ)-FuSV(θ); # function ready for bisection

    numₛₕ(θ) = vₛₕₐ(vₚₒ[l-1],vₛₒ[l-1],θ    ,ϵ[l-1],δ[l-1],γ[l-1]);
    denₛₕ    = vₛₕₐ(vₚₒ[l  ],vₛₒ[l  ],π/2.0,ϵ[l  ],δ[l  ],γ[l  ]);
    FuSH(θ)  = numₛₕ(θ)/denₛₕ;
    FSH(θ)   = sin(θ)-FuSH(θ); # function ready for bisection

    return FPP,FSV,FSH;
end


function bisection_2(
    f       :: Function,
    a       :: Float64,
    b       :: Float64;
    tol     :: AbstractFloat=1e-5,
    maxiter :: Integer=100)

    fa = f(a);
    fb = f(b);

    if fa*fb > 0
        println("No critical angle  in [0,π/2]");
        return pi/2.0; # if there is no root, then there is no crit. ang, 
                   # therefore, crit ang. is set to 90 degrees (horizontal)
    end

    i = 0;
    local c;
    while b-a > tol
        i += 1
        if i >= maxiter
            println("Max iteration exceeded")
            return c;
        end
        c = (a+b)/2.0;
        fc = f(c);
        if fc == 0.0
            break;
        elseif fa*fc > 0.0
            a  = c;  # Root is in the right half of [a,b].
            fa = fc;
        else
            b = c;  # Root is in the left half of [a,b].
        end
    end
    return c;
end

function final_Amplitudes(
    θ      :: Array{Float64,1},
    αc     :: Array{Float64,1},
    diramp :: Array{Float64,2},
    noramp :: Array{Float64,2},
    gsf    :: Float64)

    out = prod((diramp.*noramp)) * gsf;    
    critcond = any(i->(i==1), θ[1:end-1].>αc);
    if(!critcond)
        return out;
    else
        return 0.0;
    end
end     


# function critical_angles(
#     rangId :: Array{Int64,1},
#     an     :: Array{Float64,2},
#     vp     :: Array{Float64,1},
#     vs     :: Array{Float64,1},
#     ρ      :: Array{Float64,1})

#     C  = Thomsen_to_Cij(an,vp,vs,ρ);
#     C = C[:,:,rangId];     
#     vp = vp[rangId];   
#     vs = vs[rangId];    
#     ρ  = ρ[rangId];     

#     layₙ = length(ρ);
#     # calculate my C matrix for all layers: lay_1 to lay_n
#     C  = Thomsen_to_Cij(an,vp,vs,ρ);
    
#     # begin calc crit. ang. for lay_1... lay_n-1
#     Tₚₚ = zeros(Float64,layₙ-1);
#     Tₛᵥ = copy(Tₚₚ);
#     Tₛₕ = copy(Tₚₚ);
#     for lay = 2:layₙ
#         vx1 = sqrt(C[1,1,lay-1]/ρ[lay-1])
#         vz1 = sqrt(C[3,3,lay-1]/ρ[lay-1])
#         vx2 = sqrt(C[1,1,lay]/ρ[lay])
#         @show Tₚₚ[lay-1] = acotd(Fcrit(vx1,vx2,vz1));
#         vx1 = sqrt(C[5,5,lay-1]/ρ[lay-1])
#         vz1 = sqrt(C[5,5,lay-1]/ρ[lay-1])
#         vx2 = sqrt(C[5,5,lay]/ρ[lay])
#         @show Tₛᵥ[lay-1] = acotd(Fcrit(vx1,vx2,vz1));
 
#         vx1 = sqrt(C[6,6,lay-1]/ρ[lay-1])
#         vz1 = sqrt(C[5,5,lay-1]/ρ[lay-1])
#         vx2 = sqrt(C[6,6,lay]/ρ[lay])
#         @show Tₛₕ[lay-1] = acotd(Fcrit(vx1,vx2,vz1));
#     end
     
#     return Tₚₚ,Tₛᵥ,Tₛₕ
# end        
#Fcrit(vx1,vx2,vz1)= sqrt(vz1/vx1)*sqrt((vx2/vx1)^2-1.0);
#Fcrit(vx1,vx2,vz1)= sqrt((vx2/vx1)^2-1.0);



# for (idx, θ) in enumerate(θᵣ);
#     println(idx," ",θ)
#     # first medium properties: incident and reflected
#     l1_pp = s1s3_PP_c(C[:,:,lay-1],θ,ρ[lay-1],"down");
#     l1_sh = s1s3_SH_c(C[:,:,lay-1],θ,ρ[lay-1],"down");
#     l1_sv = s1s3_SV_c(C[:,:,lay-1],θ,ρ[lay-1],"down");
#     # second medium properties: refracted
#     l2_pp = s1s3_PP_c(C[:,:,lay],90.0,ρ[lay],"down");
#     l2_sh = s1s3_SH_c(C[:,:,lay],90.0,ρ[lay],"down");
#     l2_sv = s1s3_SV_c(C[:,:,lay],90.0,ρ[lay],"down");

#     @show Rₚₚ,Rₚₛ,Tₚₚ[lay-1,idx],Tₚₛ = Coef_PP_c(l1_pp,l2_pp,l1_sv,l2_sv);
#     @show Rᵥₚ,Rᵥᵥ,Tᵥₚ,Tᵥᵥ[lay-1,idx] = Coef_SV_c(l1_pp,l2_pp,l1_sv,l2_sv);
#     @show Rₕₚ,Rₕₕ,Tₕₚ,Tₕₕ[lay-1,idx] = Coef_SH_c(l1_sh,l2_sh); # h:SH and v:SV


# function s1s3_PP_c(
#     C   :: Array{Float64,2},
#     θ   :: Float64,
#     ρ   :: Float64,
#     dir :: String)

#     C₁₁   = C[1,1] ;
#     C₁₃   = C[1,3] ;
#     C₅₅   = C[5,5] ;
#     C₃₃   = C[3,3] ;
#     C₃₃⁻¹ = 1.0/C₃₃;
#     C₅₅⁻¹ = 1.0/C₅₅;

#     # eq.6.126 - Auxiliars
#     cos²    = cosd(θ)*cosd(θ)         ;
#     sin²    = sind(θ)*sind(θ)         ;
#     sin²_2θ = sind(2.0*θ)*sind(2.0*θ) ;
#     term1   = ((C₃₃ - C₅₅) * cos² - (C₁₁ - C₅₅) * sin²)^2.0 ;
#     term2   =  (C₁₃ + C₅₅) * (C₁₃ + C₅₅) * sin²_2θ          ;
#     # eq.6.126 -  Equation
#     ℂ       = sqrt(term1 + term2)                             ; # eq.6.126
#     # eq. 6.125 - Auxiliars (for SV with "-C")
#     aux = (C₅₅ + C₁₁*sin² + C₃₃*cos² + ℂ);
#     # eq. 6.125 - Equation
#     v_θ   = sqrt(aux/(2.0*ρ));
#     s₁    = sind(θ)/v_θ;           # eq. 6.124(1)

#     s₁² = s₁*s₁;

#     K₁ = ρ*(C₅₅⁻¹ + C₃₃⁻¹) + C₅₅⁻¹ * (C₁₃*C₃₃⁻¹*(C₁₃+2.0*C₅₅)-C₁₁) * s₁²; # eq.6.114(1)
#     K₂ = C₃₃⁻¹*(C₁₁*s₁²-ρ)                                              ; # eq.6.114(2)
#     K₃ = s₁²-ρ*C₅₅⁻¹                                                    ; # eq.6.114(3)

#     con1 = 1.0/sqrt(2.0);
#     if dir == "down"
#         s₃ = cosd(θ)/v_θ;       # eq. 6.124(2)
#         #s₃ = con1*sqrt(K₁ - vpsqrt(K₁^2-4.0*K₂*K₃))
#     elseif dir == "up"
#         #s₃ = -abs(cosd(θ)/v_θ); # eq. 6.124(2)
#         s₃ = -con1*sqrt(K₁ - vpsqrt(K₁^2-4.0*K₂*K₃))
#     end

#     s₃² = s₃*s₃;
#     num1 = C₅₅ * s₁² + C₃₃ * s₃² - ρ;
#     num2 = C₁₁ * s₁² + C₅₅ * s₃² - ρ;
#     den  = C₁₁ * s₁² + C₃₃ * s₃² + C₅₅*(s₁² + s₃²) - 2.0*ρ;
#     β    =  vpsqrt(num1/den); # eq 6.116
#     ζ    =  vpsqrt(num2/den); # eq 6.117

#     return (C,s₁,s₃,β,ζ);
# end

# function s1s3_SH_c(
#     C   :: Array{Float64,2},
#     θ   :: Float64,
#     ρ   :: Float64,
#     dir :: String)

#     C₄₄ = C[4,4];
#     C₄₆ = C[4,6];
#     C₆₆ = C[6,6];

#     cos²    = cosd(θ)*cosd(θ)         ;
#     sin²    = sind(θ)*sind(θ)         ;
#     sin_2θ  = sind(2.0*θ)             ;
#     ρ⁻¹      = 1.0/ρ                  ;

#     int_sqrt = C₄₄*cos² + C₆₆*sin² + C₄₆*sin_2θ  ; # eq. 1.273
#     vp_θ     = sqrt(int_sqrt * ρ⁻¹)              ; # eq. 1.273
#     c² = C₄₄*C₆₆-C₄₆*C₄₆                         ; # eq. 1.262
#     s₁ = sind(θ)/vp_θ                            ; # eq. 1.272(1)
#     if dir == "down"
#         s₃ = (1.0/C₄₄)*( -C₄₆ * s₁ + sqrt(ρ*C₄₄ - c²*s₁*s₁));
#     elseif dir == "up"
#         s₃ = (1.0/C₄₄)*( -C₄₆ * s₁ - sqrt(ρ*C₄₄ - c²*s₁*s₁));
#     end
#     β = 0.0; # irrelevant, just for output completeness
#     ζ = 0.0; # irrelevant, just for output completeness
#     return (C,s₁,s₃,β,ζ);
# end

# function s1s3_SV_c(
#     C   :: Array{Float64,2},
#     θ   :: Float64         ,
#     ρ   :: Float64         ,
#     dir :: String           )

#     C₁₁   = C[1,1] ;
#     C₁₃   = C[1,3] ;
#     C₅₅   = C[5,5] ;
#     C₃₃   = C[3,3] ;
#     C₃₃⁻¹ = 1.0/C₃₃;
#     C₅₅⁻¹ = 1.0/C₅₅;

#     # eq.6.126 - Auxiliars
#     cos²    = cosd(θ)*cosd(θ)         ;
#     sin²    = sind(θ)*sind(θ)         ;
#     sin²_2θ = sind(2.0*θ)*sind(2.0*θ) ;
#     term1   = ((C₃₃ - C₅₅) * cos² - (C₁₁ - C₅₅) * sin²)^2.0 ;
#     term2   =  (C₁₃ + C₅₅) * (C₁₃ + C₅₅) * sin²_2θ          ;
#     # eq.6.126 -  Equation
#     ℂ       = sqrt(term1 + term2)                             ; # eq.6.126

#     # eq. 6.125 - Auxiliars (for SV with "-C")
#     aux = (C₅₅ + C₁₁*sin² + C₃₃*cos² - ℂ);
#     # eq. 6.125 - Equation
#     v_θ   = sqrt(aux/(2.0*ρ));
#     s₁    = sind(θ)/v_θ;           # eq. 6.124(1)

#     s₁² = s₁*s₁;

#     K₁ = ρ*(C₅₅⁻¹ + C₃₃⁻¹) + C₅₅⁻¹ * (C₁₃*C₃₃⁻¹*(C₁₃+2.0*C₅₅)-C₁₁) * s₁²; # eq.6.114(1)
#     K₂ = C₃₃⁻¹*(C₁₁*s₁²-ρ)                                              ; # eq.6.114(2)
#     K₃ = s₁²-ρ*C₅₅⁻¹                                                    ; # eq.6.114(3)

#     con1 = 1.0/sqrt(2.0)
#     if dir == "down"
#         s₃ = cosd(θ)/v_θ;       # eq. 6.124(2)
#         #s₃ = con1*sqrt(K₁ + vpsqrt(K₁^2-4.0*K₂*K₃))
#     elseif dir == "up"
#         #s₃ = -abs(cosd(θ)/v_θ); # eq. 6.124(2)
#         s₃ = -con1*sqrt(K₁ + vpsqrt(K₁^2-4.0*K₂*K₃))
#     end

#     s₃² = s₃*s₃;
#     num1 = C₅₅ * s₁² + C₃₃ * s₃² - ρ;
#     num2 = C₁₁ * s₁² + C₅₅ * s₃² - ρ;
#     den  = C₁₁ * s₁² + C₃₃ * s₃² + C₅₅*(s₁² + s₃²) - 2.0*ρ;
#     β    =   vpsqrt(num1/den); # eq 6.116
#     ζ    =  -vpsqrt(num2/den); # eq 6.117

#     return (C,s₁,s₃,β,ζ);
# end

# function Coef_PP_c(
#     l₁ᴾ :: Tuple{Array{Float64,2},Float64,Float64,Float64,Float64},
#     l₂ᴾ :: Tuple{Array{Float64,2},Float64,Float64,Float64,Float64},
#     l₁ˢ :: Tuple{Array{Float64,2},Float64,Float64,Float64,Float64},
#     l₂ˢ :: Tuple{Array{Float64,2},Float64,Float64,Float64,Float64} )

#     C₁ᴾ  = l₁ᴾ[1]; C₁ˢ  = l₁ˢ[1];
#     s₁₁ᴾ = l₁ᴾ[2]; s₁₁ˢ = l₁ˢ[2];
#     s₃₁ᴾ = l₁ᴾ[3]; s₃₁ˢ = l₁ˢ[3];
#     β₁ᴾ  = l₁ᴾ[4]; β₁ˢ  = l₁ˢ[4];
#     ζ₁ᴾ  = l₁ᴾ[5]; ζ₁ˢ  = l₁ˢ[5];

#     C₂ᴾ  = l₂ᴾ[1]; C₂ˢ  = l₂ˢ[1];
#     s₁₂ᴾ = l₂ᴾ[2]; s₁₂ˢ = l₂ˢ[2];
#     s₃₂ᴾ = l₂ᴾ[3]; s₃₂ˢ = l₂ˢ[3];
#     β₂ᴾ  = l₂ᴾ[4]; β₂ˢ  = l₂ˢ[4];
#     ζ₂ᴾ  = l₂ᴾ[5]; ζ₂ˢ  = l₂ˢ[5];


#     W₁ᴾ,X₁ᴾ,Z₁ᴾ = WXZ(C₁ᴾ,s₁₁ᴾ,s₃₁ᴾ,β₁ᴾ,ζ₁ᴾ); # eq. 6.121
#     W₂ᴾ,X₂ᴾ,Z₂ᴾ = WXZ(C₂ᴾ,s₁₂ᴾ,s₃₂ᴾ,β₂ᴾ,ζ₂ᴾ); # eq. 6.121

#     W₁ˢ,X₁ˢ,Z₁ˢ = WXZ(C₁ˢ,s₁₁ˢ,s₃₁ˢ,β₁ˢ,ζ₁ˢ); # eq. 6.121
#     W₂ˢ,X₂ˢ,Z₂ˢ = WXZ(C₂ˢ,s₁₂ˢ,s₃₂ˢ,β₂ˢ,ζ₂ˢ); # eq. 6.121

#     mat  = [β₁ᴾ  β₁ˢ -β₂ᴾ -β₂ˢ ;   # scattering mtrix eq. 6.139
#     ζ₁ᴾ  ζ₁ˢ  ζ₂ᴾ  ζ₂ˢ ;
#     Z₁ᴾ  Z₁ˢ -Z₂ᴾ -Z₂ˢ ;
#     W₁ᴾ  W₁ˢ  W₂ᴾ  W₂ˢ];
#     RHS  =  [-β₁ᴾ;ζ₁ᴾ;-Z₁ᴾ;W₁ᴾ] ;    # eq. 6.139

#     coefs = inv(mat)*RHS;
#     Rₚₚ = coefs[1];
#     Rₚₛ = coefs[2];
#     Tₚₚ = coefs[3];
#     Tₚₛ = coefs[4];
#     return Rₚₚ,Rₚₛ,Tₚₚ,Tₚₛ;
# end

# function Coef_SH_c(l₁ˢ,l₂ˢ)
#     C₁ˢ  = l₁ˢ[1];
#     s₁₁ˢ = l₁ˢ[2];
#     s₃₁ˢ = l₁ˢ[3];

#     C₂ˢ  = l₂ˢ[1];
#     s₁₂ˢ = l₂ˢ[2];
#     s₃₂ˢ = l₂ˢ[3];

#     Zᴵ = C₁ˢ[4,6]*s₁₁ˢ + C₁ˢ[4,4]*s₃₁ˢ;  # eq.1.269
#     Zᵀ = C₂ˢ[4,6]*s₁₂ˢ + C₂ˢ[4,4]*s₃₂ˢ;  # eq.1.269
#     Zᴿ = - Zᴵ                         ;  # eq.1.281

#     Rₛₛ = (Zᴵ- Zᵀ)/(Zᵀ-Zᴿ); # eq. 1.279(1)
#     Tₛₛ = 1.0 + Rₛₛ       ; # eq. 1.277
#     Rₛₚ = 0.0             ; # eq. irrelevant, NO SH reflection to P
#     Tₛₚ = 0.0             ; # eq. irrelevant, NO SH refraction to P
#     return Rₛₚ,Rₛₛ,Tₛₚ,Tₛₛ;
# end

# function Coef_SV_c(
#     l₁ᴾ :: Tuple{Array{Float64,2},Float64,Float64,Float64,Float64},
#     l₂ᴾ :: Tuple{Array{Float64,2},Float64,Float64,Float64,Float64},
#     l₁ˢ :: Tuple{Array{Float64,2},Float64,Float64,Float64,Float64},
#     l₂ˢ :: Tuple{Array{Float64,2},Float64,Float64,Float64,Float64})

#     C₁ᴾ  = l₁ᴾ[1]; C₁ˢ  = l₁ˢ[1];
#     s₁₁ᴾ = l₁ᴾ[2]; s₁₁ˢ = l₁ˢ[2];
#     s₃₁ᴾ = l₁ᴾ[3]; s₃₁ˢ = l₁ˢ[3];
#     β₁ᴾ  = l₁ᴾ[4]; β₁ˢ  = l₁ˢ[4];
#     ζ₁ᴾ  = l₁ᴾ[5]; ζ₁ˢ  = l₁ˢ[5];

#     C₂ᴾ  = l₂ᴾ[1]; C₂ˢ  = l₂ˢ[1];
#     s₁₂ᴾ = l₂ᴾ[2]; s₁₂ˢ = l₂ˢ[2];
#     s₃₂ᴾ = l₂ᴾ[3]; s₃₂ˢ = l₂ˢ[3];
#     β₂ᴾ  = l₂ᴾ[4]; β₂ˢ  = l₂ˢ[4];
#     ζ₂ᴾ  = l₂ᴾ[5]; ζ₂ˢ  = l₂ˢ[5];


#     W₁ᴾ,X₁ᴾ,Z₁ᴾ = WXZ(C₁ᴾ,s₁₁ᴾ,s₃₁ᴾ,β₁ᴾ,ζ₁ᴾ); # eq. 6.121
#     W₂ᴾ,X₂ᴾ,Z₂ᴾ = WXZ(C₂ᴾ,s₁₂ᴾ,s₃₂ᴾ,β₂ᴾ,ζ₂ᴾ); # eq. 6.121

#     W₁ˢ,X₁ˢ,Z₁ˢ = WXZ(C₁ˢ,s₁₁ˢ,s₃₁ˢ,β₁ˢ,ζ₁ˢ); # eq. 6.121
#     W₂ˢ,X₂ˢ,Z₂ˢ = WXZ(C₂ˢ,s₁₂ˢ,s₃₂ˢ,β₂ˢ,ζ₂ˢ); # eq. 6.121

#     mat  = [β₁ᴾ  β₁ˢ -β₂ᴾ -β₂ˢ ;   # scattering mtrix eq. 6.139
#     ζ₁ᴾ  ζ₁ˢ  ζ₂ᴾ  ζ₂ˢ ;
#     Z₁ᴾ  Z₁ˢ -Z₂ᴾ -Z₂ˢ ;
#     W₁ᴾ  W₁ˢ  W₂ᴾ  W₂ˢ];
#     RHS  =  [-β₁ˢ;ζ₁ˢ;-Z₁ˢ;W₁ˢ] ;    # eq. 6.140

#     coefs = inv(mat)*RHS;
#     Rₛₚ = coefs[1];
#     Rₛₛ = coefs[2];
#     Tₛₚ = coefs[3];
#     Tₛₛ = coefs[4];
#     return Rₛₚ,Rₛₛ,Tₛₚ,Tₛₛ;
# end

# function WXZ_c(
#     C  :: Array{Float64,2},
#     s₁ :: Float64,
#     s₃ :: Float64,
#     β  :: Float64,
#     ζ  :: Float64 )

#     W = ζ*C[5,5]*s₁ + β*C[5,5]*s₃; # eq. 6.121
#     X = β*C[1,1]*s₁ + ζ*C[1,3]*s₃; # eq. 6.122
#     Z = β*C[1,3]*s₁ + ζ*C[3,3]*s₃; # eq. 6.123
#     return W,X,Z;
# end


# function vpsqrt_c(arg)
#     #principal value of sqrt
#     if arg < 0.0
#         return sqrt(Complex(arg));
#     else
#         return sqrt(arg);
#     end
# end





# # function critical_angles(an,vp,vs,ρ)
# #     P  = Thomsen_to_Cij(an,vp,vs,ρ);
# #     θᵪ = zeros(Float64,length(vp)-1);
 
# #     for i = 2:length(vp)
 
# #        p¹₄₄ = P[4,4,i-1]; p²₄₄ = P[4,4,i];
# #        p¹₆₆ = P[6,6,i-1]; p²₆₆ = P[6,6,i];
# #        p¹₄₆ = P[4,6,i-1]; p²₄₆ = P[4,6,i];
# #        ρ¹   = ρ[i-1]    ; ρ²   = ρ[i]    ;  
 
# #        pⱼ² = p¹₄₄*p¹₆₆ - p¹₄₆*p¹₄₆ ; # p² eq. 6.40 Carcione
# #        pₒ² = p²₄₄*p²₆₆ - p²₄₆*p²₄₆ ;
# #        p⁻¹₄₄ = 1.0/p¹₄₄;
# #        @show arg  = (ρ¹*p¹₄₄)/(p²₄₄*ρ²)*pₒ²-pⱼ²;
# #        if (arg >= 0.0)
# #           arg2 = p⁻¹₄₄ * (-p¹₄₆* + sqrt(arg));
# #           θᵪ[i-1] =  acotd(arg2); #eq 6.79. Solution exist in very particular situations
# #        else
# #           θᵪ[i-1] = 90.0;
# #        end
# #     end
# #     @show θᵪ
# #     error()
# #     return θᵪ;
# #  end   
 