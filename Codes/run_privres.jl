include("main_privres.jl")

sr = SOEres(Nb = 9, Na = 7, Nz = 5, Euler = true);
vfi!(sr)

save("SOEres.jld", "sr", sr)