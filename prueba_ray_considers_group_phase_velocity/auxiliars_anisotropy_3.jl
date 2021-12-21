
# auxiliars
cos²(ξ) = cos(ξ)*cos(ξ);
sin²(ξ) = sin(ξ)*sin(ξ);
sin³(ξ) = sin(ξ)*sin(ξ)*sin(ξ);
sin⁴(ξ) = sin(ξ)*sin(ξ)*sin(ξ)*sin(ξ);
# transformation relations

# eq. 1.74 Tsvankin Book
ψₚₚ(vₚₒ,vₛₒ,ψ,θ,ϵ,δ,γ) = ψ - (θ + (δ + 2.0*(ϵ-δ)*sin²(θ))*sin(2.0*θ)) ;
 # eq. 1.75 Tsvankin Book
ψₛᵥ(vₚₒ,vₛₒ,ψ,θ,ϵ,δ,γ) = ψ - (θ + 0.5*(ϵ-δ)*((vₚₒ/vₛₒ)^2)*sin(4.0*θ)) ;
# eq. 22c Weak elastic Anisotropy.
# from Seismic and acoustic velocities in Reservoir rocks Volume 2,
# Theoretical and model Studies
ψₛₕ(vₚₒ,vₛₒ,ψ,θ,ϵ,δ,γ) = atan(tan(ψ) - tan(θ)*(1.0 + 2.0*γ))               ;


# # derivatives
#
# ∇ψₚₚ(vₚₒ,vₛₒ,ψ,θ,ϵ,δ,γ) = 1.0 + 2.0*cos(2.0*θ)*(δ + 2.0*(ϵ - δ)*sin²(θ)) + 4.0*(ϵ - δ)*cos(θ)*sin(θ)*sin(2.0*θ)





###### PHASE VELOCITIES
# eq.1.61 Tsvankin Book
vₚₚₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ) = vₚₒ * (1.0 + δ * sin²(θ) + (ϵ-δ) * sin⁴(θ) )   ;
# eq.1.65 Tsvankin Book
vₛᵥₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ) = vₛₒ * (1.0 + (ϵ-δ)*((vₚₒ/vₛₒ)^2) * sin²(θ)*cos²(θ)) ;
# eq.4.16 Andreas Ruger Book
vₛₕₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ) = vₛₒ * (1.0 + γ * sin²(θ))                         ;

###### DERIVATIVES OF PHASE VELOCITIES
∇vₚₚₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ) = vₚₒ * (2.0*δ*sin(θ)*cos(θ) + 4.0*(ϵ-δ)*cos(θ)*sin³(θ) )          ;
∇vₛᵥₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ) = vₛₒ * (ϵ-δ)*((vₚₒ/vₛₒ)^2) * (2.0*sin(θ)*cos(θ))*(1.0-2.0*sin²(θ));
∇vₛₕₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ) = vₛₒ * 2.0 * γ * sin(θ)*cos(θ)                                       ;

###### SECOND DERIVATIVES OF PHASE VELOCITIES (wolfram!)
∇v²ₚₚₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ) = vₚₒ * (cos²(θ)*(2.0*δ*+12.0*(ϵ-δ)*sin²(θ)) - sin²(θ)*(2.0*δ* + 4.0*(ϵ-δ)*sin²(θ)));
∇v²ₛᵥₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ) = vₛₒ * (ϵ-δ)*((vₚₒ/vₛₒ)^2) * (4.0*sin⁴(θ)-2.0*sin²(θ)+(2.0-12.0*sin²(θ))*cos²(θ)  );
∇v²ₛₕₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ) = vₛₒ * 2.0 * γ * (cos²(θ)-sin²(θ))                                                 ;

#### GROUP VELOCITIES
#### group velocity auxiliars
∇v_vₚₚₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ) = ∇vₚₚₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ)/vₚₚₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ);
∇v_vₛᵥₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ) = ∇vₛᵥₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ)/vₛᵥₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ);

𝐕ₚₚₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ) = vₚₚₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ) * sqrt(1.0 + (∇v_vₚₚₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ))^2); # eq.170 Tsvankin Book
𝐕ₛᵥₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ) = vₛᵥₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ) * sqrt(1.0 + (∇v_vₛᵥₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ))^2); # eq.170 Tsvankin Book
𝐕ₛₕₐ(vₚₒ,vₛₒ,ψ,ϵ,δ,γ) = vₛₒ * sqrt(1.0 + 2.0*γ) / sqrt(1.0 + 2.0 * γ * cos²(ψ) )        ; # eq.169 Tsvankin Book

#### GROUP VELOCITIES DERIVATIVES
function ∇𝐕ₚₚₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ)
    V  = vₚₚₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ)   ;
    C  = ∇vₚₚₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ)  ;
    ΔV = C                       ;
    ΔC = ∇v²ₚₚₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ) ;

    num = C*ΔC + V*ΔV    ;
    den = sqrt(C^2 + V^2);
    return num/den; # by deriving eq. 1.70 in wolfram!!
end

function ∇𝐕ₛᵥₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ)
    V  = vₛᵥₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ)   ;
    C  = ∇vₛᵥₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ)  ;
    ΔV = C                       ;
    ΔC = ∇v²ₛᵥₐ(vₚₒ,vₛₒ,θ,ϵ,δ,γ) ;

    num = C*ΔC + V*ΔV    ;
    den = sqrt(C^2 + V^2);
    return num/den; # by deriving eq. 1.70 in wolfram!!
end

function ∇𝐕ₛₕₐ(vₚₒ,vₛₒ,ψ,ϵ,δ,γ)
    num = 2.0 * vₛₒ * γ * sqrt(1.0 + 2.0 * γ) * cos(ψ) * sin(ψ);
    den = (1.0 + 2.0 * γ *cos(ψ)*cos(ψ))^(1.5)                 ;
    return num/den;
end



function θ_angle(
    ψ𝐅   :: Function,
    vₚₒ  :: Float64,
    vₛₒ  :: Float64,
    ψ₀   :: Float64, # middle of interval to look for (in radians!)
    ϵ    :: Float64,
    δ    :: Float64,
    γ    :: Float64,
    porc :: Float64)

    a = ψ₀ - porc*abs(ψ₀); # looking between porc (%) difference interval
    b = ψ₀ + porc*abs(ψ₀);
    θ₀ = bisection(ψ𝐅,a,b,(vₚₒ,vₛₒ,ψ₀,ϵ,δ,γ));   # phase angle
    return θ₀
end

function bisection(
    f       :: Function,
    a       :: Float64,
    b       :: Float64,
    args    :: Tuple=();
    tol     :: AbstractFloat=1e-5,
    maxiter :: Integer=100)

    vₚₒ = args[1];
    vₛₒ = args[2];
    ψ₀  = args[3];
    ϵ   = args[4];
    δ   = args[5];
    γ   = args[6];

    fa = f(vₚₒ,vₛₒ,ψ₀,a,ϵ,δ,γ);
    fb = f(vₚₒ,vₛₒ,ψ₀,b,ϵ,δ,γ);

    if fa*fb > 0
        println("No real root in [a,b]");
        return (a+b)/2.0;
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
        fc = f(vₚₒ,vₛₒ,ψ₀,c,ϵ,δ,γ);
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
