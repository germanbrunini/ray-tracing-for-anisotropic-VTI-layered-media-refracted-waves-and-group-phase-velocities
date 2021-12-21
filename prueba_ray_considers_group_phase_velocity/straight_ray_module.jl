function straight_ray(
    vₚₒ :: Float64,
    vₛₒ :: Float64,
    an :: Array{Float64,1},
    s  :: Array{Float64,1},
    r  :: Array{Float64,1})
    # =======================================================================
    #
    # Esta rutina calcula los tiempos de arrivo de las ondas p,ssv y ssh
    # para el modelo de weak anisotropy de thompsen.
    #
    # an[1] = ϵ
    # an[2] = δ
    # an[3] = γ
    # ======================================================================

    ϵ = an[1];
    δ = an[2];
    γ = an[3];

    dist = norm(r-s);
    # Calculate ray angle (same ray for all phases, as it is a straight ray!)
    ψ₀ = atan(norm(r[1:2]-s[1:2])/abs(r[3]-s[3])); # this one is in radians!

    # For the above ray angle, calculate the phase angle (diferent phase angles!)
    θ₀ₚₚ = θ_angle(ψₚₚ,vₚₒ,vₛₒ,ψ₀,ϵ,δ,γ,0.15);   # phase angle in radians
    θ₀ₛᵥ = θ_angle(ψₛᵥ,vₚₒ,vₛₒ,ψ₀,ϵ,δ,γ,0.15);   # phase angle in radians

    # theory
    # Travel times have to be calculated througth ray trayectories
    # For this, group veocity has to be calculated. Because group
    # velocities depends on phase velocities, the phase angle (to
    # evaluate phase velocities) has to be calculated by solving
    # equations: eq. 1.74 Tsvankin Book.

    # anisotropic Ray travel time (angles in radians!)
    tₚₚ = dist/𝐕ₚₚₐ(vₚₒ,vₛₒ,θ₀ₚₚ,ϵ,δ,γ);
    tₛᵥ = dist/𝐕ₛᵥₐ(vₚₒ,vₛₒ,θ₀ₛᵥ,ϵ,δ,γ);
    tₛₕ = dist/𝐕ₛₕₐ(vₚₒ,vₛₒ,ψ₀  ,ϵ,δ,γ); # this one is on group angle! eq.1.69 Tsvankin

    return tₚₚ,tₛᵥ,tₛₕ;
end #straight_ray
