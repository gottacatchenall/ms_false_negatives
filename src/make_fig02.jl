using Makie, CairoMakie
using Distributions
using ColorSchemes
using Formatting
using MathTeXEngine: latexstring
using Colors
using DataFrames, CSV
using HypothesisTests
using StatsBase

CairoMakie.activate!(; px_per_unit=3)

diffs = CSV.read(joinpath("src", "artifacts", "joint_marginal_diffs.csv"), DataFrame)

metadata = CSV.read(joinpath("src", "artifacts", "metadata.csv"), DataFrame)

titles = Dict(
    "kolpelke_et_al_2017" => "Kopelke et al. (2017)",
    "hadfield_2014" => "Hadfield et al. (2014)",
    "havens_1992" => "Havens (1992)", 
    "ponisio_2017" => "Ponisio et al. (2017)",
    "RMBL_pollination" => "CaraDonna et al. (2015)",        
    "closs_1994" => "Closs & Lake (1994)",
    "nz_stream_foodweb" => "Townsend & Thompson (1995)";
)
 

typeofnet = Dict(
    "kolpelke_et_al_2017" => :foodweb,
    "hadfield_2014" => :hostparasite,
    "havens_1992" => :foodweb, 
    "ponisio_2017" => :plantpollinator,
    "RMBL_pollination" => :plantpollinator,     
    "closs_1994" => :foodweb,
    "nz_stream_foodweb" => :foodweb,
)


colmap = Dict(
    :foodweb=>:steelblue3,
    :plantpollinator=>:mediumseagreen,
    :hostparasite=>:mediumpurple4,
)




fig = Figure(resolution=(1200,600))
i = 1
indlist = [(1,1), (1,2), (1,3), (1,4), 
           (2,1), (2,2), (2,3)]

xpos = [0.0025, 0.000075,0.00035,0.0008, 0.00025, 0.0006,0.0005]
        
xticks = [0:0.002:0.004, 0:0.00005:0.0001, 0.0:0.00025:0.0005, 0.0:0.0005:0.001,-0.005:0.0005:0.0005, 0:0.0005:0.001, -0.001:0.001:0.001]

for (n,t) in titles
    thisdf = filter(r->r.dataset==n,diffs)
    lo, hi = min(thisdf.diff...),max(thisdf.diff...)
    mn = mean(thisdf.diff)
    sd = std(thisdf.diff)
    ylab = indlist[i][2] == 1 ? "Density" : ""

    xformatter(tickvec) = begin
        string.(tickvec)
    end
    ax = Axis(fig[indlist[i]...],xticks=xticks[i],xtickformat=xformatter,ylabel=ylab,xlabel=L"P(AB) - P(A)P(B)", title=titles[n])
    


    xlims!(ax, lo-sd, hi+sd)
    scalefactor = (filter(r->r.dataset==n,metadata).metaweb_richness[1])
    hist!(ax, thisdf.diff, offset=i,normalization=:pdf,color=(colmap[typeofnet[n]], 0.7))
    vlines!(ax, [0.], linewidth=3, color = :red, linestyle=:dash)
    density!(ax, thisdf.diff,  strokecolor =(:slategrey, 0.7), strokewidth = 3, color = (:slategray, 0.1),)


    ymin, ymax = Makie.getylimits(ax)
    xmin, xmax = Makie.getxlimits(ax)

    p = pvalue(OneSampleTTest(thisdf.diff))
    pexp = p != 0 ? Int64(floor(log10(p))+1) : 0

    str = "P < 10^{$pexp}"    

    lstr =  p == 0 ? L"P \approx 0" : latexstring(str)
    text!(ax, Point(xpos[i], ymax-(ymax/5)), text=lstr, textsize=17)
    i += 1

end


fw = PolyElement(color = (:steelblue, 0.7),)
pp = PolyElement(color = (:mediumseagreen, 0.7),)
hp = PolyElement(color = (:mediumpurple4, 0.7),)

Legend(fig[2,4], width=220, [fw,pp,hp], ["Food Web", "Plant-Pollinator", "Host-Parasite"])

fig

fig |> save("fig3.png")