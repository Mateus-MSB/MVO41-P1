include("rotacoes.jl")

# constantes
Re = 6_378 # km
μ = 398_600 # km3 / s2

# Estrutura
struct Orbita
    a::Float64
    e::Float64
    i::Float64
    ω::Float64
    Ω::Float64
    ν::Float64
end

# Bases
i = [1,0,0]
j = [0,1,0]
k = [0,0,1]

using LinearAlgebra
unit(vetor) = vetor/norm(vetor) # vetor unitario

# Metodos
function get_elements_from_rv(R::Vector,V::Vector)
    r = norm(R)
    v = norm(V)
    
    h = cross(R,V)
    ϵ = v^2/2 - μ/r
    B = cross(V,h) - μ/r * R
    e_vec = B/μ
    N = cross(k,h)

    a = -μ/2ϵ
    e = norm(e_vec)
    I = acosd(dot(k,unit(h)))

    ω = acosd(dot(unit(e_vec),unit(N)))
    if(dot(k,e_vec)<=0)
        ω = 360 - ω
    end

    Ω = acosd(dot(i,unit(N)))
    if(dot(j,N)<=0)
        Ω = 360 - Ω
    end
    
    ν = acosd(dot(unit(R),unit(e_vec)))
    if(dot(R,V)<=0)
        ν = 360 - ν
    end

    orbit = Orbita(a,e,I,ω,Ω,ν)
    show_orbit_elements(orbit)
    
    return orbit
end

function get_perifocal_from_elements(a,e,ν)
    p = a*(1-e^2)
    r = p/(1+e*cosd(ν))

    R = [r*cosd(ν);r*sind(ν);0]
    V = [-sqrt(μ/p)*sind(ν);sqrt(μ/p)*(e+cosd(ν));0]
    
    println("Rₚ: ",R)
    println("Vₚ: ",V,"\n")
    return R, V
end

function get_eci_from_perifocal(R_pf::Vector,V_pf::Vector,i,ω,Ω)
    R_eci = Rz(-Ω) * Rx(-i) * Rz(-ω) * R_pf 
    V_eci = Rz(-Ω) * Rx(-i) * Rz(-ω) * V_pf 

    println("Rₓ: ",R_eci)
    println("Vₓ: ",V_eci,"\n")

    return R_eci, V_eci
end 

function get_eci_from_elements(a,e,ν,i,ω,Ω)
    (R_pf,V_pf) = get_perifocal_from_elements(a,e,ν)

    (R_eci,V_eci) = get_eci_from_perifocal(R_pf,V_pf,i,ω,Ω);

    return R_eci, V_eci
end

function get_RaDec_from_r(R::Vector)
    r = norm(R)
    l,m,n = R/r
    
    δ = asind(n)
    α = acosd(l/cosd(δ))
    if(m<=0)
        α = 360 - α
    end

    println("α: ",α)
    println("δ: ",δ,"\n")

    return α, δ
end

# Print
function show_orbit_elements(orbit::Orbita)
    println("a: ",orbit.a)
    println("e: ",orbit.e)
    println("i: ",orbit.i)
    println("ω: ",orbit.ω)
    println("Ω: ",orbit.Ω)
    println("ν: ",orbit.ν)
end

