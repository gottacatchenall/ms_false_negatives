using Plots
using EcologicalNetworks
using StatsBase
using Statistics: std

include("_common.jl")

function add_falsenegatives(truenetwork, p_fn)
    mat = adjacency(truenetwork)
    for (i,el) in enumerate(mat)
        if el == 1 && rand() < p_fn
            mat[i] = 0 
        elseif el == 0 && rand() < p_fn
            mat[i] = 1
        end
    end
    return UnipartiteNetwork(mat)
end

function errors(f;
    S = 100,
    fnr = 0:0.02:1,
    nreps= 30
)

    err = zeros(length(fnr))
    sd_err = zeros(length(fnr))
    for (fni,p_fn) in enumerate(fnr)
        errorlist = zeros(nreps)

        for rep in 1:nreps

            truenetwork = flexiblelinksmodel(S)
            observednet = add_falsenegatives(truenetwork, p_fn)

            errorlist[rep] = sqrt((f(truenetwork) - f(observednet) )^2)
        end
        err[fni] = mean(filter(x -> !isnan(x), errorlist))
        sd_err[fni] = std(filter(x -> !isnan(x), errorlist))
    end

    return fnr, err, sd_err
end

meandegcent(x) = mean(collect(values(centrality_degree(x))))
meankatzcent(x) = mean(collect(values(centrality_katz(x))))
meanclosenessctr(x) = mean(collect(values(centrality_closeness(x))))

fnr, connect_err, connect_sd = errors(connectance, nreps=2000)
fnr, degcent_err, degcent_sd = errors(meandegcent, nreps=2000)

fnr, specrad_err, specrad_sd = errors(ρ, nreps=300) 


fnt = font(20, "Roboto")

plot(dpi=300, fontfamily=fnt, size=(500,500), frame=:box, ylim=(0,1),xlim=(0,1))
plot!(fnr, connect_err, label="", ribbon=connect_sd, c=:dodgerblue4, fa=0.2, )
plot!(fnr, connect_err, label="", ribbon=2connect_sd, c=:dodgerblue4, fa=0.1)
scatter!(fnr, connect_err, mc=:white, msc=:dodgerblue4, msw=1.5, ms=4,label="connectance")

plot!(fnr, specrad_err, label="", ribbon=specrad_sd, c=:gray50, fa=0.2)
plot!(fnr, specrad_err, label="", ribbon=2specrad_sd, c=:gray50, fa=0.1)
scatter!(fnr, specrad_err, mc=:white, msw=1.5, ms=4, msc=:gray50,label="spectral radius")


plot!(fnr, degcent_err, label="", ribbon=degcent_sd, c=:teal, fa=0.2)
plot!(fnr, degcent_err, label="", ribbon=2degcent_sd, c=:teal, fa=0.1)
scatter!(fnr, degcent_err, mc=:white, msw=1.5, ms=4, msc=:teal,label="degree centrality")
xaxis!("false negative rate")
yaxis!("mean-squared-error in property estimate")



savefig("props_specrad.png")