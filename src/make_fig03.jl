using Makie, CairoMakie
using DataFrames
using Colors
using StatsBase
using CSV

CairoMakie.activate!(; px_per_unit=2)

model_data = CSV.read(joinpath("src", "artifacts", "model_comparison.csv"), DataFrame)
trophic_data = CSV.read(joinpath("src", "artifacts", "trophic_levels.csv"), DataFrame)

fig = Figure(resolution=(900,900))
ga = GridLayout(fig[1,1])

ax1 = Axis(ga[1,1],
    xticks=0:0.1:1,
    yticks=0:0.1:1,
    xminorticks=0:0.05:1,
    xminorgridvisible=true,
    xminorticksvisible=true,
    yminorticks=0:0.05:1,
    yminorgridvisible=true,
    yminorticksvisible=true,
    xlabel="Added FNR",
    ylabel="ROC-AUC",
    xlabelsize=22,
    ylabelsize=22
)
limits!(ax1, 0, 0.95, 0,1)


grey, red, orange, green = ([colorant"grey35", colorant"red", colorant"orange", colorant"green"])

addalpha(x,alpha) = RGBA(x.r,x.g,x.b, alpha)

grey, red, orange, green = addalpha.([grey, red, orange, green], 0.3) 

poly!(ax1, 
    Rect(0, 0, 1.0, 0.5),
    color = (:grey, 0.3),
    transparency = true
)
text!(ax1, 0.63, 0.46, text="Worse than random", textsize=12, justification=:center)


poly!(ax1, 
    Rect(0, 0.5, 1.0, 0.25),
    color = red,
    transparency = true
)
text!(ax1, 0.68, 0.52, color=:firebrick, textsize=12, text="Close to random")

poly!(ax1, 
    Rect(0, 0.75, 1.0, 0.15),
    color = orange,
    transparency = true
)
text!(ax1, 0.72, 0.85, color=:salmon2, textsize=12, text="Fair Classifier")

poly!(ax1, 
    Rect(0, 0.9, 1.0, 0.1),
    color = green,
    transparency = true
)
text!(ax1, 0.65, 0.94, color=:darkgreen, textsize=12, text="Excellent Classifier")

models = unique(model_data.model)

colors = Colors.JULIA_LOGO_COLORS
marks = [:circle, :utriangle, :dtriangle, :rect]
for (i,m) in enumerate(models)
    thisdf = filter(r->r.model==m, model_data)
    scatterlines!(ax1, thisdf.added_fnr, thisdf.rocauc,  linewidth=2, color=colors[i],strokewidth=1.5, strokecolor=colors[i], markercolor=:white, markersize=15, marker = marks[i])
end

gb = GridLayout(fig[1,2])
ax2 = Axis(
    gb[1,1],
    xticks=0:0.1:1,
    yticks=0:0.1:1,
    xminorticks=0:0.05:1,
    xminorgridvisible=true,
    xminorticksvisible=true,
    yminorticks=0:0.05:1,
    yminorgridvisible=true,
    yminorticksvisible=true,
    xlabel="Added FNR",
    ylabel="PR-AUC",
    xlabelsize=22,
    ylabelsize=22
)
limits!(ax2, 0, 0.95, 0,1)


grey, red, orange, green = ([colorant"grey35", colorant"red", colorant"orange", colorant"green"])

addalpha(x,alpha) = RGBA(x.r,x.g,x.b, alpha)

grey, red, orange, green = addalpha.([grey, red, orange, green], 0.3) 

poly!(ax2, 
    Rect(0, 0, 1.0, 0.5),
    color = (:grey, 0.3),
    transparency = true
)
text!(ax2, 0.63, 0.46, text="Worse than random",  textsize=12,justification=:center)


poly!(ax2, 
    Rect(0, 0.5, 1.0, 0.25),
    color = red,
    transparency = true
)
text!(ax2, 0.68, 0.52, color=:firebrick,  textsize=12,text="Close to random")

poly!(ax2, 
    Rect(0, 0.75, 1.0, 0.15),
    color = orange,
    transparency = true
)
text!(ax2, 0.72, 0.85, color=:salmon2,  textsize=12, text="Fair Predictor")

poly!(ax2, 
    Rect(0, 0.9, 1.0, 0.1),
    color = green,
    transparency = true
)
text!(ax2,0.65, 0.94, color=:darkgreen,  textsize=12, text="Excellent Predictor")

models = unique(model_data.model)

colors = Colors.JULIA_LOGO_COLORS
marks = [:circle, :utriangle, :dtriangle, :rect]
for (i,m) in enumerate(models)
    thisdf = filter(r->r.model==m, model_data)
    scatterlines!(ax2, thisdf.added_fnr, thisdf.prauc,  linewidth=2, color=colors[i],strokewidth=1.5, strokecolor=colors[i], markercolor=:white, markersize=15, marker = marks[i])
end


gc = GridLayout(fig[2,1])

ax3 = Axis(gc[1,1], 
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
limits!(ax3, 0,1,1,8)

cols = [
    colorant"seagreen4",
    colorant"dodgerblue3",
    colorant"mediumpurple4",
    colorant"indianred3",
]

function get_mns()
    NS = [100, 250, 500, 1000]
    allmns, allsds = [], []
    for (i,ns) in enumerate(NS)
        thisdf = filter(r->r.numspecies == ns, trophic_data)
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

NS = [100, 250, 500, 1000]
mns, sds = get_mns()
for (i,ns) in enumerate(NS)
    fnrs = 0:0.05:1
    c = cols[i]    
    band!(ax3, 
        fnrs,  
        mns[i] .+ sds[i],  
        mns[i] .- sds[i],  
        transparency=true,
        color=RGBA(c.r,c.g,c.b,0.3))  
end 

fnrs = 0:0.05:1

for (i,ns) in enumerate(NS)
    scatterlines!(ax3, 
        fnrs,
        mns[i],
        markercolor=:white, 
        color=cols[i],
        strokewidth=1.5, 
        markersize=9,
        strokecolor=cols[i],
        label=ns)
end 



gd = GridLayout(fig[2,2])
els = [PolyElement(color=c) for c in cols]


modelels = [MarkerElement(marker=marks[i],strokecolor=colors[i], markercolor=:white, strokewidth=1.5, markersize=12) for i in 1:4]

Legend(gd[1,1],modelels, models,  "Models")
Legend(gd[2,1],els, string.(NS),  "Species Richness")


colsize!(fig.layout, 1, Relative(0.5))
rowsize!(fig.layout, 1, Relative(0.5))
colsize!(fig.layout, 2, Relative(0.5))
rowsize!(fig.layout, 2, Relative(0.5))

fig

for (label, layout) in zip(["A", "B", "C"], [ga, gb, gc])
    Label(layout[1, 1, TopLeft()], label,
        textsize = 22,
        font = "TeX Gyre Heros Bold",
        padding = (10, 10, 10, 10),
        halign = :right)
end

fig |> save(joinpath("figures", "fig3.png"))