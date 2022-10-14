using CairoMakie, Distributions
using MathTeXEngine
using MathTeXEngine: latexstring
CairoMakie.activate!(; px_per_unit=2)


p_a = 0.01:0.01:0.5
p_b = 0.01:0.01:0.5


exp_obs = zeros(length(p_a), length(p_b))

expectednumobs(target_negatives, p_a, p_b) = target_negatives/(p_a*p_b)


# p_a = 0.1, p_b = 0.2
# how many of all species to see this pair 5 times?
# 0.1  *  0.2

# to see A and B once, you need an expected 1/(p_a*p_b) exp_obs

targ = 10
for I in CartesianIndices(size(exp_obs))
    exp_obs[I] = expectednumobs(targ, p_a[I[1]], p_a[I[2]])
end

function maketicks(x)
    latexstring.(["10^{$(i)}" for i in x])
end

fig = Figure(resolution=(800,700))
ax = Axis(fig[1,1], xlabel=L"P(A)", ylabel=L"P(B)", xlabelsize=22,
ylabelsize=22)
limits!(ax,0.01,0.5,0.01,0.5)
co = contourf!(ax, p_a, p_b, log10.(exp_obs), colorlimits=1:5, levels=range(1,5, length = 10))
Colorbar(fig[1,2], co, ticksize=20, label="Expected number of observations of all individuals to see A and B together 10 times",ticks=[i for i in 1:0.5:5], tickformat=maketicks)




fig

fig |> save("new_1b.png")
