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
    xticks=0:0.1:0.5,
    xminorticks = 0:0.05:0.5,
    yticks=0:0.1:0.5,
    yminorticks = 0:0.05:0.5,
    xminorticksvisible = true,
    xminorgridvisible = true,
    xlabel=L"P(A)",
    ylabel=L"P(B)",
    xlabelsize=22,
    ylabelsize=22,
    yticksize=28)


p_a = 0.01:0.01:0.5
p_b = 0.01:0.01:0.5
expectednumobs(target_negatives, p_a, p_b) = target_negatives/(p_a*p_b)
exp_obs = zeros(length(p_a), length(p_b))
targ = 10
for I in CartesianIndices(size(exp_obs))
    exp_obs[I] = expectednumobs(targ, p_a[I[1]], p_a[I[2]])
end

function maketicks_1b(x)
    x
    #latexstring.(["10^{$(i)}" for i in x])
end


cmap_1b = parse.(Colorant, [
    "#fa6d3f",
    "#f69e43",
    "#e4ba4f",
    "#a4b167",
    "#31934d",
    "#3e9c8a",
    "#19caa7",
    #"#3b606c",
    "#37414c"])

limits!(ax2,0.01,0.5,0.01,0.5)
co = contourf!(
    ax2, 
    p_a, 
    p_b, 
    log10.(exp_obs), 
    colorlimits=1:5, 
   # colormap = ColorSchemes.BuGn_9,
    colormap = ColorSchemes.devon,
    levels=range(1,5, length = 10)
)
Colorbar(gb[1,2], co, ticksize=20, label="Expected number of observations of all individuals to see A and B together 10 times",ticks=[i for i in 1:0.5:5], tickformat=maketicks_1b)





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

simulated_df = CSV.read(joinpath("src", "artifacts", "1c_simulated_fnr.csv"), DataFrame)

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


mangal_df = CSV.read(joinpath("src", "artifacts", "1d_mangal_fnrs.csv"), DataFrame)

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

#f |> save(joinpath("figures","fig1.png"))