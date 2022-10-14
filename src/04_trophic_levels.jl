using Makie, CairoMakie
using EcologicalNetworks
using StatsBase
using Statistics: std
using ProgressMeter
using DataFrames
using Colors
using CSV
using Random
Random.seed!(42)

include("flexiblelinks.jl")

function add_falsenegatives(truenetwork, p_fn)
    mat = adjacency(truenetwork)
    for (i,el) in enumerate(mat)
        if el == 1 && rand() < p_fn
            mat[i] = 0 
        elseif el == 1 
            mat[i] = 1
        end
    end
    return size(truenetwork,1) == size(truenetwork,2) ? UnipartiteNetwork(mat) : BipartiteNetwork(mat)
end

function make_perfect_nested(T,B)
    mat = zeros(Bool, T, B)
    
    for i in 1:T
        for j in 1:B
            if i > j
                mat[i,j] = 1
            end 
        end
    end
    BipartiteNetwork(mat)
end 


S = 150
nreps = 50

trophic_means, trophic_vars = Float32[], Float32[]
fnrs = Float32[]

df = DataFrame(trophic_means=Float32[], 
    trophic_vars=Float32[], index=[], fnr=Float32[], numspecies=Float32[])

NS = [100, 250, 500, 1000]
for ns in NS
    @info ns
    truenetwork = flexiblelinksmodel(ns)
    for (ind,p_fn) in enumerate(0:0.05:1)
        for i in 1:50
            observednet = add_falsenegatives(truenetwork, p_fn)
            vals = collect(values(trophic_level(observednet)))
            push!(df.fnr,p_fn )
            push!(df.trophic_means, mean(vals))
            push!(df.trophic_vars, var(vals))
            push!(df.numspecies, ns)
            push!(df.index, ind)
        end 
    end 
end 

fig = Figure(resolution=(700,700))

ax1 = Axis(fig[1,1], 
    xticks=0:0.1:1,
    xlabel="Added FNR",
    ylabel="Mean trophic level",
    xminorgridvisible=true,
    xminorticks=0:0.05:1,
    yminorticksvisible=true,
    yminorticks=0:0.5:1,
    xlabelsize=22,
    ylabelsize=22,
)
limits!(ax1, 0,1,1,8)
df

findall(!isfinite, df.trophic_means)

cols = [
    colorant"seagreen4",
    colorant"dodgerblue3",
    colorant"mediumpurple4",
    colorant"indianred3",
]

function get_mns()
    allmns, allsds = [], []
    for (i,ns) in enumerate(NS)
        thisdf = filter(r->r.numspecies == ns, df)
        #scatterlines!(ax1, thisdf.fnr, thisdf.trophic_means)
        mns = Float32[]
        sds = Float32[]
        fnrs = 0:0.05:1
        for (i,p_fn) in enumerate(fnrs)
            thisfnr = filter(r->r.index == i, thisdf)
            push!(mns, mean(filter(!isnan, thisfnr.trophic_means)))
            push!(sds, mean(sqrt.(thisfnr.trophic_vars)))
        end 
        push!(allmns, mns)
        push!(allsds, sds)
    end
    allmns, allsds
end 

mns, sds = get_mns()
for (i,ns) in enumerate(NS)
    fnrs = 0:0.05:1
    c = cols[i]    
    band!(ax1, 
        fnrs,  
        mns[i] .+ sds[i],  
        mns[i] .- sds[i],  
        transparency=true,
        color=RGBA(c.r,c.g,c.b,0.3))  
end 

for (i,ns) in enumerate(NS)
    scatterlines!(ax1, 
        fnrs,
        mns[i],
        markercolor=:white, 
        color=cols[i],
        strokewidth=2.5, 
        markersize=12,
        strokecolor=cols[i],
        label=ns)
end 


els = [PolyElement(color=c) for c in cols]
Legend(fig[1,2],els, string.(NS),  "Species Richness")


fig

CSV.write(joinpath("src", "artifacts", "trophic_levels.csv"), df)