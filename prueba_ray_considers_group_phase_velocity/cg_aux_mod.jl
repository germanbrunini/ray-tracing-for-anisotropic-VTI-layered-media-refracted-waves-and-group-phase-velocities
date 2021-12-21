# ==============================================================================
function func(
  x_in :: Array{Float64,1},
  e    :: Array{Float64,1},
  s    :: Array{Float64,1},
  r    :: Array{Float64,1},
  vₚₒ  :: Array{Float64,1},
  vₛₒ  :: Array{Float64,1},
  an   :: Array{Float64,2},
  ray  :: String) # output: Tiempo de viaje and Δt intervals

  len_e = size(e,1)              ;
  xin1  = size(x_in,1)           ;
  x     = zeros(Float64,xin1 + 2);

  x[1         ] = 0.0                ;
  x[2:xin1 + 1] = x_in               ;
  x[xin1   + 2] = norm(r[1:2]-s[1:2]);

  func    = 0.0                      ;
  mem_fun = zeros(Float64,(len_e,1)) ;
  for i = 2:len_e + 1
    Δx   = abs(x[i-1]-x[i]); # radial distance between ray-breaking points
    ψ   = atan(Δx/e[i-1])  ; # in radians!
    vₐˣ = an_vel(ψ,vₚₒ[i-1],vₛₒ[i-1],an[i-1,:],ray)*sin(ψ); # radial velocity!
    func         = func + Δx / vₐˣ; # sumes up all Δt
    mem_fun[i-1] = Δx / vₐˣ       ; # stores each Δt
  end
  return func,mem_fun;
end #function func

# ==============================================================================
function dfunc(
  x_in :: Array{Float64,1},
  e    :: Array{Float64,1},
  s    :: Array{Float64,1},
  r    :: Array{Float64,1},
  vp   :: Array{Float64,1},
  vs   :: Array{Float64,1},
  an   :: Array{Float64,2},
  ray  :: String)

  # !Para calcular la derivada de un punto x[i] se necesita también el punto
  # x(i+1) y el punto x[i-1], por lo tanto se van a calcular las derivadas
  # desde i=2 a i=N-1, siendo x(N) el receptor; x[1] es la fuente.

  xin1  = size(x_in,1)           ;
  x     = zeros(Float64,xin1 + 2);
  dfunc = zeros(Float64,xin1    );

  x[1         ] = 0.0      ;
  x[2:xin1 + 1] = x_in     ;
  x[xin1   + 2] = norm(r[1:2]-s[1:2]);

  for i = 2:size(x,1)-1
    dx1   = x[i]   - x[i-1]              ;
    dx2   = x[i+1] - x[i]                ;
    phi1  = atan(dx1 / e[i-1])           ;
    phi2  = atan(dx2 / e[i])             ;
    dphi1 = e[i-1]   / (e[i-1]*e[i-1] + dx1*dx1);
    dphi2 = -e[i]    / (e[i]  *e[i]   + dx2*dx2);
    v1    = an_vel(phi1, vp[i-1], vs[i-1], an[i-1,:], ray);
    v2    = an_vel(phi2, vp[i]  , vs[i]  , an[i,:]  , ray);
    sp1   = sin(phi1);
    sp2   = sin(phi2);
    dv1   = an_dvel(phi1, vp[i-1], vs[i-1], an[i-1,:], ray);
    dv2   = an_dvel(phi2, vp[i]  , vs[i]  , an[i,:]  , ray);

    aux1  = cos(phi1)*v1 + sp1*dv1;
    aux2  = cos(phi2)*v2 + sp2*dv2;

    den1 = (sp1*v1)^2;
    den2 = (sp2*v2)^2;

    num1 = sp1*v1 - dx1*dphi1*aux1;
    num2 = sp2*v2 + dx2*dphi2*aux2;
    dfunc[i-1] = (num1/den1)-(num2/den2);

  end
  return dfunc;
end #function dfunc

# ==============================================================================

function an_vel(
  ψₒ  :: Float64, # ray angle
  vₚₒ :: Float64,
  vₛₒ :: Float64,
  an  :: Array{Float64,1},
  ray :: String)

  an_vel = 0.0;
  ϵ = an[1];
  δ = an[2];
  γ = an[3];

  if (ray=="vp")
    # For the above ray angle, calculate the phase angle (diferent phase angles!)
    θ₀ₚₚ = θ_angle(ψₚₚ,vₚₒ,vₛₒ,ψₒ,ϵ,δ,γ,0.30) ;  # phase angle in radians
    vₐ   = 𝐕ₚₚₐ(vₚₒ,vₛₒ,θ₀ₚₚ,ϵ,δ,γ)         ;
  elseif (ray=="sv")
    θ₀ₛᵥ = θ_angle(ψₛᵥ,vₚₒ,vₛₒ,ψₒ,ϵ,δ,γ,0.30); # phase angle in radians
    vₐ   = 𝐕ₛᵥₐ(vₚₒ,vₛₒ,θ₀ₛᵥ,ϵ,δ,γ)          ;
  elseif (ray=="sh")
    vₐ   = 𝐕ₛₕₐ(vₚₒ,vₛₒ,ψₒ,ϵ,δ,γ)            ;
  else
    @warn("input wrong key word argument in an_vel")
  end
  return vₐ;
end # function an_vel

# ==============================================================================

function an_dvel(
  ψₒ  :: Float64, # ray angle
  vₚₒ :: Float64,
  vₛₒ :: Float64,
  an  :: Array{Float64,1},
  ray :: String)

  an_dvel = 0.0;
  ϵ = an[1];
  δ = an[2];
  γ = an[3];

  if (ray=="vp")
    θ₀ₚₚ = θ_angle(ψₚₚ,vₚₒ,vₛₒ,ψₒ,ϵ,δ,γ,0.30) ;  # phase angle in radians
    Δvₐ  = ∇𝐕ₚₚₐ(vₚₒ,vₛₒ,θ₀ₚₚ,ϵ,δ,γ)          ;
  elseif (ray=="sv")
    θ₀ₛᵥ = θ_angle(ψₛᵥ,vₚₒ,vₛₒ,ψₒ,ϵ,δ,γ,0.30) ;  # phase angle in radians
    Δvₐ  = ∇𝐕ₛᵥₐ(vₚₒ,vₛₒ,θ₀ₛᵥ,ϵ,δ,γ)          ;
  elseif (ray=="sh")
    Δvₐ = ∇𝐕ₛₕₐ(vₚₒ,vₛₒ,ψₒ,ϵ,δ,γ)             ;
  end
  return Δvₐ;
end # function an_dvel
