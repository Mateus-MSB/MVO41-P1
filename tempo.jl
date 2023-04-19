using Roots

# constantes
Re = 6_378 # km
μ = 398_600 # km3 / s2

"""
    Compute time in hours
"""
function get_time(a::Float64, e::Float64, ν::Float64, τ::Float64=0.0)

    E = 2*atan(tand(ν/2)*sqrt((1-e)/(1+e)))
    M = E - e*sin(E)
    
    n = sqrt(μ/a^3)
    t = τ + M/n
    return t/3600
end

"""
    Compute true anomaly in degrees
"""
function get_anomaly(t::Float64,a::Float64,e::Float64)
    n = sqrt(μ/a^3)
    M = n*t

    F(E) = M - E + e*sin(E)
    if(M < π)
        E0 = M + e/2
        E = find_zero(F,E0)
    else
        E0 = M - e/2
        E = find_zero(F,E0)
    end
    println(E)
    ν = 2 * atand(tan(E/2)*sqrt((1+e)/(1-e)))
    println(ν)
    
    if(ν<0)
        ν = 360 + ν
    end
    return ν
end