using EcologicalNetworks
using Distributions
using ProgressMeter
using CairoMakie

include(joinpath("src", "utils.jl"))

ids = map(i -> i.ID, filter(i -> contains(i.Reference, "Hadfield"), web_of_life()))
B = convert.(BipartiteNetwork, web_of_life.(ids))
M = reduce(∪, B)

function get_metadata()
    C = zeros(Int64, size(M)) 
    N = zeros(Int64, size(M))  

    # We then move through the metaweb step by step
    for (i,si) in enumerate(species(M; dims=1)), (j,sj) in enumerate(species(M; dims=2))
        for b in B
            if si in species(b)
                if sj in species(b)
                    C[i,j] += 1

                    i_b_idx, j_b_idx = findall(x->x==si, b.T)[begin], findall(x->x==sj, b.B)[begin]
                    if b[i_b_idx, j_b_idx] == 1
                        N[i,j] += 1
                    end
                end
            end
        end
    end
    return C,N
end 

C,N = get_metadata()

probabilistic_network = BipartiteProbabilisticNetwork(rand(size(M)...)) # Prob net object from EcologicalNetworks
RAs = [rand(LogNormal(0,1), i) for i in size(probabilistic_network)]
RAs = [r ./ sum(r) for r in RAs] 

relabd = RelativeAbundance([M.T, M.B], RAs)
cooc = CooccurrenceMatrix([M.T,M.B], C)
ints = InteractionsMatrix([M.T,M.B], N)

exp_fpr = Exponential(0.02)
cooccurrence_sample = InteractionUncertaintySampler.sample(CooccurrenceSampler(relabd), 1000)
mean_fnr = FNR(cooccurrence_sample, C, M)


using DataFrames, CSV
df = CSV.read("modelfit.csv", DataFrame)

# predicted prob metaweb
P = zeros(size(M))
for r in eachrow(df)
    i = findfirst(i->M.T[i]==string(r.speciesI), 1:length(M.T))
    j = findfirst(i->M.B[i]==string(r.speciesJ), 1:length(M.B))
    P[i,j] = r.prediction
end 

function resample(p, fnr, fpr, num_particles=100, num_samples=200) 
    samples = []
    for _ in 1:num_samples
        s = 0
        for _ in 1:num_particles
            observed = rand(Bernoulli(p))
            if !observed
                is_fn = rand(Bernoulli(fnr))

                s += is_fn 
            elseif observed
                is_fp = rand(Bernoulli(fpr))

                s += !is_fp
            end
        end 
        push!(samples, s/num_particles)
    end 
    samples
end

normal_approximation(pᵢⱼ, fnr, fpr, num_particles) = begin
    p_adj = pᵢⱼ*(1-fpr) + (1-pᵢⱼ)*fnr
    σ = sqrt(p_adj*(1-p_adj))/sqrt(num_particles)
    if sqrt(p_adj*(1-p_adj)) < 0.
        @info p_adj
        @info "sigma <= 0"
        return Normal(p_adj, 0.)
    end
    Normal(p_adj, σ)
end 

fpr = 0.02
num_particles = 150
num_samples = 1000

ipr = InteractionProbabilityResampler(fnr_prior, Exponential(mean_fpr))

hosts = ["Lemmus sibiricus", "Lagurus lagurus","Talpa altaica"]
parasites = ["Xenopsylla cheopis", "Coptopsylla olgae", "Stenoponia ivanovi"]
h_j, p_i = [2,9,17], [1,160,90]

f = Figure(resolution=(1100,1600))
f[1,1] = g = GridLayout()
par_colors = [:dodgerblue, :red, :forestgreen]
commonnames = ["Little ground squirrel", "European snow vole", "Wood mouse"]

for axnum in 1:3
    ax1 = Axis(
        g[axnum,2],
        xticks=0:0.1:1,
        xminorticksvisible = true,
        xminorticks=0:0.05:1,
        xminorgridvisible=true,
        yticklabelsvisible = false,
        yticksvisible = false,
        xlabel = axnum == 3 ? L"P(i \leftrightarrow j)" : "",
        xlabelsize=32,
        titlesize=22,
        subtitlesize=16,
        titlealign=:left,
        subtitle=M.B[h_j[axnum]],
        title=commonnames[axnum]
        )
    hideydecorations!(ax1, ticks = false)
    xlims!(ax1, 0,1)
    for i in 1:length(p_i)
        hlines!(ax1,[-5i], linewidth=0.4,color=(par_colors[i], 0.9))
    end

    gauss(i, mu, var) = 1.0/(sqrt(var)*2π) * exp(-1*((i-mu)/(2*sqrt(var))^2))
    x = 0:0.01:1
    for (i, p) in reverse(collect(enumerate(p_i)))
        pij = P[p,h_j[axnum]]
        exp_pstar = min(pij*(2-mean_fnr)*(1-mean_fpr), 1)

#        limiting_dist = Normal(exp_pstar,
#        sqrt(exp_pstar*(1-exp_pstar))/sqrt(num_particles))

        # TODO 
        # correct num_particles based on obs co-occ

        limiting_dist = normal_approximation(pij, mean_fnr, mean_fpr, num_particles)


        y =  [pdf(limiting_dist, i) for i in x]
        lines!(ax1, x, -5i.+y,linewidth=2,color = (par_colors[i], 1))
    end

    for (i, p) in reverse(collect(enumerate(p_i)))
        pij = P[p,h_j[axnum]]
        vector = resample(pij, mean_fnr, mean_fpr, num_particles, num_samples)
        hist!(ax1,vector, normalization=:pdf, offset = -5i, color = (par_colors[i], 0.5),
            bandwidth = 0.2)
    end

    for (i, p) in reverse(collect(enumerate(p_i)))
        #vlines!(ax1,P[p,h_j[axnum]], linewidth=3,color = (par_colors[i], 0.5))
        linesegments!(ax1, [(P[p,h_j[axnum]], -5i), (P[p,h_j[axnum]], -5i+15)], linestyle=:dash,linewidth=3,color = (par_colors[i], 1))
    end

end 

f

leg = f[1,1]


p1 = PolyElement(color = (par_colors[1], 0.7),)
p2 = PolyElement(color = (par_colors[2], 0.7),)
p3 = PolyElement(color = (par_colors[3], 0.7),)


paranames = [M.T[i] for i in p_i]

Legend(f[1,1][1,1], labelsize=22, width=400, [p1,p2,p3],paranames, "Parasite")

colgap!(g, 200)
f


f |> save(joinpath("resampled_uncert.png"))

#= old 

 plts = []
for i in p_i
    plt = plot()
    for j in h_j
    
        histogram!(plt, [tmp[x][i,j] for x in 1:length(tmp)], label="$i,$j", bins=20)
        vline!(plt, label="", [P[i,j]], lw=5)
    end
    push!(plts, plt)
end

plot(plts...)

i = 1, j = 2
i = 1, j = 9
i = 1, j = 17
i = 160, j = 17
i = 162, j = 17

j = 17
tmpvec = [mean([tmp[x][i,j] for x in 1:length(tmp)]) for i in 1:length(M.T)]
scatter(1:length(tmpvec), tmpvec)

M.T[i] =#