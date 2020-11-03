include("main_privres.jl")

sr = SOEres(Nb = 9, Na = 7, Nz = 5, Euler = true);
vfi!(sr, verbose=true)

save("SOEres.jld", "sr", sr)