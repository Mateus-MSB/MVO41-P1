using Roots

# constantes
Re = 6_378 # km
μ = 398_600 # km3 / s2

"""
    Compute time in hours
"""
function get_time(a::Real, e::Real, ν::Real, τ::Real=0.0)
    a = abs(a)
    
    if (e<0)
        error("e CANNOT BE NEGATIVE")

    elseif (e<1)
        E = 2*atan(tand(ν/2)*sqrt((1-e)/(1+e)))
        M = E - e*sin(E)
        
        n = sqrt(μ/a^3)
        t = τ + M/n
        return t/3600

    elseif (e==1)
        p = a*(1-e^2)
        D = sqrt(p)*tand(ν/2)
        t = τ + 1/(2*sqrt(μ)) * (p*D + D^3/3)
        return t/3600

    else
        F = 2*atanh(tand(ν/2)*sqrt((e-1)/(e+1)))
        Mh = e*sinh(F) - F
        
        n = sqrt(μ/a^3)
        t = τ + Mh/n
        return t/3600
    end
end

"""
    Compute true anomaly in degrees.
    Inform time in hours.
    If the orbit is parabolic, then the second argument if h, not a.
"""
function get_anomaly(t::Real,a::Real,e::Real)
    t = 3600t # hours to seconds
    a = abs(a)

    if (e<0)
        error("e CANNOT BE NEGATIVE")

    elseif (e<1)
        n = sqrt(μ/a^3)
        M = n*t

        F(E) = M - E + e*sin(E)
        if(M < π) # Algoritmo 3.1 do Curtis
            E0 = M + e/2
            E = find_zero(F,E0)
        else
            E0 = M - e/2
            E = find_zero(F,E0)
        end
        ν = 2 * atand(tan(E/2)*sqrt((1+e)/(1-e)))
        
        if (ν<0)
            ν = 360 + ν
        end

        return ν

    elseif (e==1)
        Mp = μ^2 * t / a^3
        ν = 2 * atand((3*Mp + sqrt(9Mp^2+1))^(1/3) - (3*Mp + sqrt(9Mp^2+1))^(-1/3))
        return ν

    else
        n = sqrt(μ/a^3)
        Mh = n*t
        @show Mh

        G(F) = e*sinh(F) - F - Mh

        F = find_zero(G,0)
        @show F
        ν = 2 * atand(tanh(F/2)*sqrt((e+1)/(e-1)))
        
        @show ν
        if (ν<0)
            ν = 360 + ν
        end
        return ν 
    end
end
