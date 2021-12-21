function synth_generator(
	gdata :: Tuple{Float64,Int64,Int64,Float64,Float64},
    tₜ    :: Float64, # anisotropic travel time
	ℒ     :: Float64, # Geometrical Spreading due to last segment
	ϕ     :: Float64, # azimuth for rotation around symetry axis "z=x3"
    vᵣₑ   :: Float64, # anisotropic phase velocity on receiver
    vₛₒ   :: Float64, # anisotropic phase velocity on source
    ρᵣₑ   :: Float64, # density on receiver
    ρₛₒ   :: Float64, # density on source
    M     :: Array{Float64,2}, # moment tensor
	𝐏ₛ    :: Array{Float64,1}, # phase slowness vector
	𝐔ᵣ    :: Array{Float64,1}, # Unit phase polarization vector at receiver
	𝐔ₛ    :: Array{Float64,1}, # Unit phase polarization vector at source
    phase :: String)::Array{Float64,2}; # which phase to be calculated "pp"sh"sv"

	Δₜ  = gdata[1];  # time sampling
    ns  = gdata[2];  # number of samples in trace
    nr  = gdata[3];  # number rof receivers, in this case, one at a time
    t₀ₛ = gdata[4];  # initial source time (usually set to 0.0)
    f₀  = gdata[5];  # source frequency

	# we use Aki richards way (1/v^3), so only use direction of P (|P|=1)
	# Gretchka book uses (1/v^2) and |P|=1/v.
	𝐏ₛ = transpose(𝐏ₛ/norm(𝐏ₛ)); # make it row vector
	𝐔ᵣ = transpose(𝐔ᵣ)        ;
	𝐔ₛ = transpose(𝐔ₛ)        ;

	# now, "ap" is on the propagation plane [x₁,x₂,x₃], which contains
	# the symetry axis z=x₃. In fact, ap exists in [x₁,x₃] due to
	# VTI anysitropy (page 322 of Gretchka book). It needs to be
	# rotated towards the xyz(north east down) system.

	# we need a clockwise "cw" rotation. This is because we need to
	# see the vectors move back (counter clockwise), and therefore,
	# we need to advance the system clockwise (put the vectors behind!).

	𝐏ₛ = R𝑧(ϕ,𝐏ₛ,"cw");
	𝐔ᵣ = R𝑧(ϕ,𝐔ᵣ,"cw");
	𝐔ₛ = R𝑧(ϕ,𝐔ₛ,"cw");

	# Microseismic Monitoring. Vladimir Gretchka. 2017. Pag.109.
	# Altought the results expresed in equation 3.80 is obtained for homogeneos
	# isotropic media, it turns out to be valid universally, that is, for arbitrary
	# anisotropy and heterogeneity (e.g. Cerveny 2001, Chapter 2; also chapetr4)

    # values are multiplied by a constant, for numerical stability
	# (I don't like seismograms to have small numbers) but, this
	# value can be change as desired. (set to 1.0, for example)

	scale = 1e3      ;
	#𝐏ₛ  = 𝐏ₛ  * scale; # PD: V is divided -> P is multiplied!
	vᵣₑ = vᵣₑ / scale;
	vₛₒ = vₛₒ / scale;
    ρᵣₑ = ρᵣₑ / scale;
    ρₛₒ = ρₛₒ / scale;

    vₛₒ²  = vₛₒ*vₛₒ             ; # eq. 4.93 Aki-Rich & eq.4.108 Gretchka Book
	rw    = Ricker(dt=Δₜ, f0=f₀);
    nrw   = length(rw)          ;
    rw    = reshape(rw,(1,nrw)) ;
    nf    = 4 * nextpow(2,ns)   ;
    dw    = 2.0 * pi/(nf*Δₜ)    ;
    rwpad = hcat(rw,zeros(nf-nrw)');
    RW    = fft(rwpad)             ;

    U     = zeros(Float64, (ns, 3*nr));

    if phase=="pp"
        @inbounds for j = 1:nr
			# eq.4.93                Aki-Richards and
			# eq 3.73 / 3.84 / 4.108 Gretchka Book
			cp = 1.0/(4.0 * pi * sqrt(ρₛₒ*ρᵣₑ * vₛₒ*vᵣₑ) * ℒ * vₛₒ²);

			ℛᵈ  = 𝐔ₛ' * M * 𝐏ₛ ; # eq.6.23 Gretchka Book.
			URᵈ = 𝐔ᵣ * ℛᵈ      ; # eq.6.22 Gretchka Book.
			ap  = cp * URᵈ     ;

            Uj  = zeros(Complex{Float64}, nf, 3);
            Ujp = zeros(Complex{Float64}, nf, 3);
            nf2 = floor(Int, nf/2);
            for k = 2:nf2
                w    = (k-1)*dw;
                imw  = im*w;
                imwR = imw*RW[k];
                tsp  = t₀ₛ + tₜ;
                Ujp[k, :]     = ap*exp(-imw*tsp);
                Uj[k, :]      = imwR * (Ujp[k,:]);
                Uj[nf+2-k, :] = conj(Uj[k, :]);
            end
            uj = real(ifft(Uj, 1))[1:ns, :];
            jx = 3*(j-1) + 1;
            jy = jx + 1     ;
            jz = jx + 2     ;
            U[:, jx:jz] = uj;
    	end
		return U;
    end

    if phase=="sh"
        @inbounds for j = 1:nr
			# eq.4.93                Aki-Richards and
			# eq 3.73 / 3.84 / 4.108 Gretchka Book
			cs = 1.0/(4.0 * pi * sqrt(ρₛₒ*ρᵣₑ * vₛₒ*vᵣₑ) * ℒ * vₛₒ²);

			ℛᵈ  = 𝐔ₛ' * M * 𝐏ₛ ; # eq.6.23 Gretchka Book.
			URᵈ = 𝐔ᵣ * ℛᵈ      ; # eq.6.22 Gretchka Book.
			as  = cs * URᵈ     ;

            Uj  = zeros(Complex{Float64}, nf, 3);
            Ujs = zeros(Complex{Float64}, nf, 3);
            nf2 = floor(Int, nf/2);
            for k = 2:nf2
                w    = (k-1)*dw;
                imw  = im*w;
                imwR = imw*RW[k];
                tss  = t₀ₛ + tₜ;
                Ujs[k, :] = as*exp(-imw*tss);
                Uj[k, :]  = imwR * (Ujs[k,:]);
                Uj[nf+2-k, :] = conj(Uj[k, :]);
            end
            uj = real(ifft(Uj, 1))[1:ns, :];
            jx = 3*(j-1) + 1;
			jy = jx + 1     ;
            jz = jx + 2;
            U[:, jx:jz] = uj;
		end
		return U;
    end

    if phase=="sv"
        @inbounds for j = 1:nr
			# eq.4.93                Aki-Richards and
			# eq 3.73 / 3.84 / 4.108 Gretchka Book
			cs = 1.0/(4.0 * pi * sqrt(ρₛₒ*ρᵣₑ * vₛₒ*vᵣₑ) * ℒ * vₛₒ²);

			ℛᵈ  = 𝐔ₛ' * M * 𝐏ₛ ; # eq.6.23 Gretchka Book.
			URᵈ = 𝐔ᵣ * ℛᵈ      ; # eq.6.22 Gretchka Book.
			as  = cs * URᵈ     ;

            Uj  = zeros(Complex{Float64}, nf, 3);
            Ujs = zeros(Complex{Float64}, nf, 3);
            nf2 = floor(Int, nf/2);
            for k = 2:nf2
                w    = (k-1)*dw;
                imw  = im*w;
                imwR = imw*RW[k];
                tss  = t₀ₛ + tₜ;
                Ujs[k, :] = as*exp(-imw*tss);
                Uj[k, :]  = imwR * (Ujs[k,:]);
                Uj[nf+2-k, :] = conj(Uj[k, :]);
            end
            uj = real(ifft(Uj, 1))[1:ns, :];
            jx = 3*(j-1) + 1;
			jy = jx + 1     ;
            jz = jx + 2;
            U[:, jx:jz] = uj;
		end
		return U;
    end

end

"""
Ricker(; <keyword arguments>)
Create a Ricker wavelet.
# Keyword arguments
* `dt::Real=0.002`: sampling interval in secs.
* `f0::Real=20.0`: central frequency in Hz.
# Examples
```julia
julia> w = Ricker(); plot(w);
julia> w = Ricker(dt=0.004, f0=20); plot(w);
```
# Reference
Sheriff, Robert, 2002, Encyclopedic Dictionary of Applied Geophysics, fourth
ed.: Society of Exploration Geophysicists. Geophysical Reference Series No. 13.
"""

function Ricker(;dt::Real=0.002, f0::Real=20.0)::Array{Float64,1}
    nw = 2.0/(f0*dt);
    nc = floor(Int, nw/2);
    t  = dt*collect(-nc:1:nc);
    p = [f0];

    # writedlm("dRick",dRick.(t,p));
    # writedlm("Rick",Rick.(t,p));
    return Rick.(t,p);
end

# Ricker
Rick(t,p)  = (1.0-2.0*(pi*p[1]*t)^2)*exp(-(pi*p[1]*t)^2);
