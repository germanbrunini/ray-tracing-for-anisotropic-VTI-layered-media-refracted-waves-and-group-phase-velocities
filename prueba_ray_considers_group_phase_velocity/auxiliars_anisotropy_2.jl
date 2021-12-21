function VTI_UnitPolVec(
    C   :: Array{Float64,2},
    θ   :: Float64         ,
    ρ   :: Float64         ,
    dir :: String          ,
    key :: String);

    # based on equation 6.115 of Carccione's Book
    # WARNING: This vectors are contained in the
    # ray propagation plane, not in the cartesian
    # one. So, in order to analize them in xyz, it
    # is necesarily to rotate them according to ray
    # direction.

     if key=="pp"

        nada,s₁,s₃,β,ζ = s1s3_PP(C,θ,ρ,dir)
        v  = 1.0/sqrt(s₁*s₁+s₃*s₃); # equivalent to phase velocity eq.(4.20)Ruger Book
        p  = [s₁,0.0,s₃]; # slowness vector for VTI
        n  = p*v;         # eq 2.58 Gretcka book
        C𝜌 = Cρ(C,ρ) ;
        Γ  = ΓΓ(n,C𝜌);
        vₚₚ,vₛₕ,vₛᵥ,Uₚₚ,Uₛₕ,Uₛᵥ = Christoffel(Γ);
        #U = [β,0.0,ζ] # equation 6.115 of Carccione's Book
        return n,p,Uₚₚ;

    elseif key =="sh"
        # in VTI, sh polarization is ortogonal
        # to ray propagation plane.

        nada,s₁,s₃,β,ζ = s1s3_SH(C,θ,ρ,dir);
        v = 1.0/sqrt(s₁*s₁+s₃*s₃); # equivalent to phase velocity eq.(4.20)Ruger Book
        p = [s₁,0.0,s₃]; # slowness vector for VTI
        n = p*v;         # eq 2.58 Gretcka book
        C𝜌 = Cρ(C,ρ) ;
        Γ  = ΓΓ(n,C𝜌);
        vₚₚ,vₛₕ,vₛᵥ,Uₚₚ,Uₛₕ,Uₛᵥ = Christoffel(Γ);
        return n,p,Uₛₕ;

    elseif key=="sv"

        nada,s₁,s₃,β,ζ = s1s3_SV(C,θ,ρ,dir);
        v = 1.0/sqrt(s₁*s₁+s₃*s₃); # equivalent to phase velocity eq.(4.20)Ruger Book
        p = [s₁,0.0,s₃]; # slowness vector for VTI
        n = p*v;         # eq 2.58 Gretcka book
        C𝜌 = Cρ(C,ρ) ;
        Γ  = ΓΓ(n,C𝜌);
        vₚₚ,vₛₕ,vₛᵥ,Uₚₚ,Uₛₕ,Uₛᵥ = Christoffel(Γ);
        return n,p,Uₛᵥ;
        #U = [β,0.0,ζ] # equation 6.115 of Carccione's Book

    end
end


function ΓΓ(
    l :: Array{Float64,1},
    C :: Array{Float64,2})

    # equation 1.69 Carcione Book. "wave field in Real..."
    # and 2.47 gretchka book.
    # (3X1)*(3X3)*(1X3) = (3X3)
    Γ = zeros(Float64,(3,3))

    C₁₁=C[1,1]; C₁₂=C[1,2]; C₁₃=C[1,3]; C₁₄=C[1,4]; C₁₅=C[1,5]; C₁₆=C[1,6];
    C₂₁=C₁₂   ; C₂₂=C[2,2]; C₂₃=C[2,3]; C₂₄=C[2,4]; C₂₅=C[2,5]; C₂₆=C[2,6];
    C₃₁=C₁₃   ; C₃₂=C₂₃   ; C₃₃=C[3,3]; C₃₄=C[3,4]; C₃₅=C[3,5]; C₃₆=C[3,6];
    C₄₁=C₁₄   ; C₄₂=C₂₄   ; C₄₃=C₃₄   ; C₄₄=C[4,4]; C₄₅=C[4,5]; C₄₆=C[4,6];
    C₅₁=C₁₅   ; C₅₂=C₂₅   ; C₅₃=C₃₅   ; C₅₄=C₄₅   ; C₅₅=C[5,5]; C₅₆=C[5,6];
    C₆₁=C₁₆   ; C₆₂=C₂₆   ; C₆₃=C₃₆   ; C₆₄=C₄₆   ; C₆₅=C₅₆   ; C₆₆=C[6,6];

    l₁  = l[1] ; l₂  = l[2] ; l₃  = l[3] ;
    l₁² = l₁*l₁; l₂² = l₂*l₂; l₃² = l₃*l₃;
    # equation 1.73 Carccione Book.
    Γ[1,1] = C₁₁*l₁² + C₆₆*l₂² + C₅₅*l₃² + 2.0*C₅₆    *l₂*l₃ + 2.0*C₁₅    *l₃*l₁ + 2.0*C₁₆    *l₁*l₂;
    Γ[2,2] = C₆₆*l₁² + C₂₂*l₂² + C₄₄*l₃² + 2.0*C₂₄    *l₂*l₃ + 2.0*C₄₆    *l₃*l₁ + 2.0*C₂₆    *l₁*l₂;
    Γ[3,3] = C₅₅*l₁² + C₄₄*l₂² + C₃₃*l₃² + 2.0*C₃₄    *l₂*l₃ + 2.0*C₃₅    *l₃*l₁ + 2.0*C₄₅    *l₁*l₂;
    Γ[1,2] = C₁₆*l₁² + C₂₆*l₂² + C₄₅*l₃² + (C₄₆ + C₂₅)*l₂*l₃ + (C₁₄ + C₅₆)*l₃*l₁ + (C₁₂ + C₆₆)*l₁*l₂;
    Γ[1,3] = C₁₅*l₁² + C₄₆*l₂² + C₃₅*l₃² + (C₄₅ + C₃₆)*l₂*l₃ + (C₁₃ + C₅₅)*l₃*l₁ + (C₁₄ + C₅₆)*l₁*l₂;
    Γ[2,3] = C₅₆*l₁² + C₂₄*l₂² + C₃₄*l₃² + (C₄₄ + C₂₃)*l₂*l₃ + (C₃₆ + C₄₅)*l₃*l₁ + (C₂₅ + C₄₆)*l₁*l₂;

    Γ[2,1] = copy(Γ[1,2]);
    Γ[3,1] = copy(Γ[1,3]);
    Γ[3,2] = copy(Γ[2,3]);

    return Γ;
end

function Cρ(
    C :: Array{Float64,2},
    ρ :: Float64)
    #equation 2.41 gretchka book
    return C/ρ; #
end

function Christoffel(Γ::Array{Float64,2})
    F = eigen(Γ); # solve eigenvalue problem

    # set by descending order, v₁>v₂>v₃. eq. 2.50 gretcka book.
    λ =    F.values[3:-1:1] ; # eigenvalues
    Λ = F.vectors[:,3:-1:1] ; # eigenvectors

              v = sqrt.(λ);
    vₚₚ,vₛₕ,vₛᵥ = v;

    ΔVₛ = 100.0*abs(vₛₕ-vₛᵥ)/abs(vₛᵥ);

    # if both S-velocities are similar, (happens often in anisotropy,
    # and always in isotropy)
    if ΔVₛ < 1.0
        # then, manually compute eigenvectors!
        # This is because, setting them in order can be twisted!

        Uₚₚ = Λ[:,1];         # this is obvioulsy related to biggest √λ
        Uₛₕ = [0.0,1.0,0.0];  # this is ortoghonal to [x1,x3] plane
                              # equation 6.115 of Carccione's Book
        Uₛᵥ = [Uₚₚ[2]*Uₛₕ[3]-Uₛₕ[2]*Uₚₚ[3], # this is ortoghonal to other two
                      -Uₚₚ[1]*Uₛₕ[3]+Uₛₕ[1]*Uₚₚ[3],
                              Uₚₚ[1]*Uₛₕ[2]-Uₛₕ[1]*Uₚₚ[2]];
        return vₚₚ,vₛₕ,vₛᵥ,Uₚₚ,Uₛₕ,Uₛᵥ;

    else

        Uₚₚ = Λ[:,1];         # this is obvioulsy related to biggest √λ
        Uₛₕ = Λ[:,2];  # this is ortoghonal to [x1,x3] plane
        Uₛᵥ = Λ[:,3];
        return vₚₚ,vₛₕ,vₛᵥ,Uₚₚ,Uₛₕ,Uₛᵥ;

    end
end

function R𝑧(
    ϕ   :: Float64,
    v   :: Transpose{Float64,Array{Float64,1}},
    key :: String)

    if (size(v)==(1,3))
        v = transpose(v);
    end

    if (key=="cw") #clock-wise rotation
    # rotation around simetry z-axis.
    # equation 1.59 of Carcione's Book.
    R = [cosd(ϕ) -sind(ϕ) 0.0;
         sind(ϕ)  cosd(ϕ) 0.0;
          0.0     0.0   1.0 ];
        return R*v;
    elseif (key=="ccw") #counter clock-wise rotation
        # rotation around simetry z-axis.
        # equation 1.59 of Carcione's Book.
        R = [cosd(-ϕ) -sind(-ϕ) 0.0;
            sind(-ϕ)  cosd(-ϕ) 0.0;
            0.0        0.0    1.0 ];
        return R*v;
    end
end


function R𝑧(
    ϕ   :: Float64,
    v   :: Array{Float64,1},
    key :: String)

    if (size(v)==(1,3))
        v = transpose(v);
    end

    if (key=="cw") #clock-wise rotation
    # rotation around simetry z-axis.
    # equation 1.59 of Carcione's Book.
    R = [cosd(ϕ) -sind(ϕ) 0.0;
         sind(ϕ)  cosd(ϕ) 0.0;
          0.0     0.0   1.0 ];
        return R*v;
    elseif (key=="ccw") #counter clock-wise rotation
        # rotation around simetry z-axis.
        # equation 1.59 of Carcione's Book.
        R = [cosd(-ϕ) -sind(-ϕ) 0.0;
            sind(-ϕ)  cosd(-ϕ) 0.0;
            0.0        0.0    1.0 ];
        return R*v;
    end
end

#
# function K(
#     C   :: Array{Float64,2},
#     ρ   :: Float64         ,
#     s₁                     );
#
#     # Based on Andreas Ruger book page 50. Similar (expressions equal to Carcione)
#     C₁₁   = C[1,1] ;
#     C₁₃   = C[1,3] ;
#     C₅₅   = C[5,5] ;
#     C₃₃   = C[3,3] ;
#     C₃₃⁻¹ = 1.0/C₃₃;
#     C₅₅⁻¹ = 1.0/C₅₅;
#
#     C₃₅⁻¹  = 1.0/(C₃₃*C₅₅)
#     C₁₃₅² = (C₁₃+C₅₅)*(C₁₃+C₅₅)
#     s₁² = s₁*s₁;
#
#     # eq.4.22 Andreas Ruger
#     K₁ = ρ*(C₅₅⁻¹ + C₃₃⁻¹) - (C₁₁*C₅₅⁻¹ + C₅₅*C₃₃⁻¹ - (C₃₅⁻¹)*C₁₃₅²)*s₁² ;
#     K₂ = C₃₃⁻¹*(C₁₁*s₁²-ρ)                                               ;
#     K₃ = s₁²-ρ*C₅₅⁻¹                                                     ;
#
#     return K₁,K₂,K₃
# end
#
# function q(K₁,K₂,K₃,key::String)
#     c = 1.0/sqrt(2.0);
#     # based on equation 4.21 of Andreas rugger Book
#     if key=="pp"
#         q𝛼 = c*sqrt(K₁-sqrt(K₁*K₁-4.0*K₂*K₃)); #equation
#         return q𝛼;
#     elseif key=="sv"
#         qᵦ = c*sqrt(K₁+sqrt(K₁*K₁-4.0*K₂*K₃));
#         return qᵦ
#     end
# end

#
#
#
# function UnitpolarizationVectors(traj,C,ID,ρ,vᵣ,vₛ,key::String)
#
#     # unit wavefront normal at source (eq.2.43 gretcka)
#     nₛ = (traj[2,:]-traj[1,:])'      /norm(traj[2,:]-traj[1,:])       ;
#     # unit wavefront normal at receiver (eq.2.43 gretcka)
#     nᵣ = (traj[end,:]-traj[end-1,:])'/norm(traj[end,:]-traj[end-1,:]) ;
#
#     idᵣ = ID[1]; # id of receiver layer
#     idₛ = ID[2]; # id of source layer
#
#     Cρₛ = Cρ(C[:,:,idₛ],ρ[idₛ]);
#     Cρᵣ = Cρ(C[:,:,idᵣ],ρ[idᵣ]);
#
#
#     Γₛ    = ΓΓ(nₛ,Cρₛ);
#     Γᵣ    = ΓΓ(nᵣ,Cρᵣ);
#     Vₛ,Uₛ = Christoffel(Γₛ);
#     Vᵣ,Uᵣ = Christoffel(Γᵣ);
#
#     if     key=="pp"
#         return Uₛ[:,1], Uᵣ[:,1];
#     elseif key=="sh"
#         return Uₛ[:,2], Uᵣ[:,2];
#     else
#         return Uₛ[:,3], Uᵣ[:,3];
#     end
# end
