## Includes
include("orbitas.jl")
include("tempo.jl")

# constantes
Re = 6_378 # km
μ = 398_600 # km3 / s2

## Example 4.1
r = [-5368, -1784, 3691]
RA, Dec = get_RaDec_from_r(r);


## Example 4.3
r = [-6045,-3490,2500]
v = [-3.457,6.618,2.533]

orbita43 = get_elements_from_rv(r,v);

## Example 4.7
h = 80_000
e = 1.4
a = h^2/(μ*(1-e^2))

i = 30
Ω = 40
ω = 60
ν = 30

R_perifocal, V_perifocal = get_perifocal_from_elements(a,e,ν);
R_eci, V_eci = get_eci_from_perifocal(R_perifocal,V_perifocal,i,ω,Ω);
R_eci2, V_eci2 = get_eci_from_elements(a,e,ν,i,ω,Ω);

## Example 3.1
rp = 9600
ra = 21000

a = 0.5*(rp+ra)
e = (ra-rp)/(rp+ra)
ν = 120.0
t = get_time(a,e,ν)

## Example 3.2
t_new = 3.0*3600 # s

ν_new = get_anomaly(t_new,a,e)