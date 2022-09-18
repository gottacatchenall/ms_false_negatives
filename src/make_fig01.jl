using Makie, CairoMakie
using Distributions
using ColorSchemes
using Formatting
using MathTeXEngine: latexstring
using Colors
using DataFrames, CSV

#=
    Top left panel:

    Geometric distribution
=#

CairoMakie.activate!(; px_per_unit=2)

pdf(x, fnr) = 1.0 - fnr^x
xs = 1:15
fnrs = 0

f = Figure(resolution = (1300, 1200))
ga = GridLayout(f[1,1])

ax1 = Axis(ga[1,1],
    xticks = 1:15,
    yticks = 0:0.2:1,
    yminorticks = 0:0.1:1,
    yminorticksvisible = true,
    yminorgridvisible = true,
    xlabel = "Number of observed negatives", 
    ylabel = "Probability interaction is a true-negative",
    xlabelsize=22,
    ylabelsize=22)

for (i,fnr) in enumerate([0.1,0.3,0.5,0.7,0.9])
    scatterlines!(ax1, xs, pdf.(xs, fnr),
        strokewidth = 2,
        linewidth = 4,
        markersize = 15, 
        color = ColorSchemes.Paired_7[i],
        strokecolor = ColorSchemes.Paired_7[i], 
        markercolor = :white,
        label = "$fnr")
end
limits!(ax1, (1,15), (0,1))
axislegend("False Negative Rate", position=:rb)

#= 
    Top right panel:

    Expected observations of null  

=#

function maketex(x)
    #prev = ["10^$(Int32((log10(xi))))" for xi in x]
    #latexstring.(prev)

    string.("10^",Int32.(log10.(x)))
end

gb = GridLayout(f[1,2])

ax2 = Axis(gb[1,1],
    xticks=0:0.1:1,
    yticks=[10^i for i in 0:7],
    xminorticks = 0:0.05:0.5,
    xminorticksvisible = true,
    xminorgridvisible = true,
    ytickformat=maketex,
    xlabel="Relative abundance of particular species",
    ylabel="Expected number of observations of all species",
    yscale=log10,
    xlabelsize=22,
    ylabelsize=22,
    yticksize=28)

limits!(ax2, 0, 0.5, 1, 10^7)
function expected_samples(numobs)
    relativeabund = Float64[]
    obs = Float64[]
    for relabd in vcat(0.001:0.001:0.5)
        push!(relativeabund, relabd)
        push!(obs, numobs/relabd)
    end
    relativeabund, obs
end

for (i,numobs) in enumerate([1,10,100,1000,10000])
    relabd, req = expected_samples(numobs)
    lines!(ax2,relabd, req,
        strokewidth = 2,
        linewidth = 4,
        markersize = 12, 
        color = ColorSchemes.tableau_sunset_sunrise[i],
        strokecolor = ColorSchemes.tableau_sunset_sunrise[i], 
    #    markercolor = :white,
        label="$numobs"
    )

end
axislegend("Target obs. of particular species",orientation = :horizontal, position=:rt)



gc = GridLayout(f[2,1])

ax3 = Axis(gc[1,1],
    xlabel="Total number of observed interactions across all species",
    ylabel="Realised false negative rate",
    xticks=0:250:1500,
    xminorgridvisible=true,
    xminorticks=0:100:1500,
    yticks=0:0.2:1,
    yminorticksvisible=true,
    yminorticks=0:0.1:1,
    xlabelsize=22,
    ylabelsize=22
)
limits!(ax3, 0,1500, 0, 1)

simulated_df = CSV.read(joinpath("..", "artifacts", "1c_simulated_fnr.csv"), DataFrame)

richnesses = unique(simulated_df.richness)
#cs1 = ColorScheme(range(colorant"dodgerblue", colorant"cyan4",length=3))
cs1= ColorSchemes.cmyk[[1,4,8]]
for (i,r) in enumerate(richnesses)
    thisdf = filter(x->x.richness == r, simulated_df)

    c = cs1[i]

    band!(ax3, 
        thisdf.sampling_effort,  
        thisdf.mean_fnr + 2thisdf.sd_fnr,  
        thisdf.mean_fnr- 2thisdf.sd_fnr,
        transparency=true,
        color=RGBA(c.r,c.g,c.b,0.3))  
    
    band!(ax3, 
        thisdf.sampling_effort,  
        thisdf.mean_fnr + thisdf.sd_fnr,  
        thisdf.mean_fnr- thisdf.sd_fnr,
        transparency=true,
        color=RGBA(c.r,c.g,c.b,0.3))  
    scatterlines!(ax3, 
        thisdf.sampling_effort,
        thisdf.mean_fnr,
        strokewidth = 2,
        linewidth = 4,
        markersize = 10,
        label="$r",
        color=c,
        markercolor=:white,
        strokecolor=cs1[i])    
end
axislegend("Richness", position=:rt)


gd = GridLayout(f[2,2])

ax4 = Axis(gd[1,1],
    xlabel="Total number of observed interactions across all species",
    ylabel="Realised false negative rate",
    xticks=0:250:1500,
    xminorgridvisible=true,
    xminorticks=0:100:1500,
    yticks=0:0.2:1,
    yminorticksvisible=true,
    yminorticks=0:0.1:1,
    xlabelsize=22,
    ylabelsize=22,
)

limits!(ax4, 0,1500, 0,1)


mangal_df = CSV.read(joinpath("..", "artifacts", "1d_mangal_fnrs.csv"), DataFrame)

maxrich = max(unique(mangal_df.richness)...)
cs2 = ColorScheme(range(colorant"dodgerblue", colorant"cyan4",length=maxrich))


I = unique(mangal_df.index)

for i in I
    thisdf = filter(x->x.index == i, mangal_df)
    c =  cs2[thisdf.richness[begin]]
    lines!(ax4,
        thisdf.sampling_effort,
        thisdf.mean_fnr,
        color=RGBA(c.r,c.g,c.b,0.3))
end

Colorbar(gd[1,2], label="Richness", colormap=cs2[1:end], limits=(0,maxrich), ticks=0:100:700)


for (label, layout) in zip(["A", "B", "C", "D"], [ga,gb,gc,gd])
    Label(layout[1, 1, TopLeft()], label,
        textsize = 22,
        font = "TeX Gyre Heros Bold",
        padding = (10, 10, 10, 10),
        halign = :right)
end

f


f |> save(joinpath("..", "figures","fig1.png"))