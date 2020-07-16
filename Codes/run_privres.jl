include("main_privres.jl")

sr = SOEres();

vfi!(sr)

save("SOEres.jld", "sr", sr)