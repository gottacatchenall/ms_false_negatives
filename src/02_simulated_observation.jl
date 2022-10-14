using EcologicalNetworks
using StatsBase: mean, std
using Distributions: LogNormal, Categorical
using ColorSchemes, Colors
using Plots
using Measures
using Random
Random.seed!(42)

include(joinpath("utils", "flexiblelinks.jl"))
include(joinpath("utils", "all_mangal_networks.jl"))

function observe(numobservations, relativeabundances)
    counts = zeros(length(relativeabundances))
    for i in 1:numobservations
        ind = rand(Categorical(relativeabundances))
        counts[ind] += 1
    end
    return counts
end

function countfp(A, observations)
    metaweb = adjacency(A)
    S = richness(A)
    tp, fn = 0,0
    for i in 1:S, j in 1:S
        if metaweb[i,j] == 1
            if observations[i] > 0 && observations[j] > 0
                tp += 1
            elseif observations[i] > 0 || observations[j] > 0
                fn += 1
            end
        end
    end
    return tp, fn
end

function getabundances(A)
    S = richness(A)
    # By trophic levels: 
    # Z = 2
    # abundances = Z.^[trophdict["s$i"]-1 for i in 1:S]
    # abundance_dist = abundances ./ sum(abundances)

    # By lognormal dist:
    abundances = rand(LogNormal(),S)
    return  abundances ./ sum(abundances)
end

function runreplicate(numobservations, A)
    relativeabundances = getabundances(A)
    observations = observe(numobservations, relativeabundances)
    tp, fn = countfp(A, observations)

    fnr = fn/(fn+tp)
    return fnr
end

function buildplot(A;  
    numreplicates = 200,
    samplingeffort = vcat(1,25, 50:50:1500))


    means = zeros(length(samplingeffort))
    sds = zeros(length(samplingeffort))

    for (s, numobs) in enumerate(samplingeffort)
        thesefrs = zeros(numreplicates)
        for r in 1:numreplicates
            # sample metaweb if as is a function for generating,
            # otherwise it is an empirical web so just use it 
            metaweb = typeof(A) <: Function ? A(1) : A             
            thesefrs[r] = runreplicate(numobs, metaweb)
        end
        means[s] = mean(thesefrs);
        sds[s] = std(thesefrs);
    end

    means, sds
end



samp = vcat(1,25, 50:50:1500)
numreplicates = 200
mean250, sd250 = buildplot(flexiblelinksmodel(250), numreplicates=numreplicates)
mean100, sd100 = buildplot(flexiblelinksmodel(100), numreplicates=numreplicates)
mean50, sd50 = buildplot(flexiblelinksmodel(50), numreplicates=numreplicates)

richnesses = vcat([fill(x, length(samp)) for x in [50,100,250]]...)
means = vcat(mean50, mean100, mean250)
sds = vcat(sd50, sd100, sd250)

sampeffort = [samp..., samp..., samp...]

df = DataFrame(richness=richnesses, sampling_effort=sampeffort, mean_fnr=means, sd_fnr=sds)
CSV.write(joinpath("src", "artifacts", "1c_simulated_fnr.csv"), df)

#= 
    Mangal data


=#

querymangal()


writemangaldata()

fw, para, mutu, misc = mangaldata()
means_per_fw = []
rich = [richness(f) for f in fw]

@showprogress for thisfw in fw
    means, sds = buildplot(thisfw)
    push!(means_per_fw, means)
end

samp = vcat(1,25, 50:50:1500)
means = vcat(means_per_fw[sortperm(rich)]...)

ind = 1:length(rich)

xs = [fill(x,length(samp)) for x in rich[sortperm(rich)]]
richnesses = vcat(xs...)

sampeffort = vcat([samp for i in 1:length(rich)]...)

indices = vcat([fill(ind[i], length(samp)) for i in 1:length(rich)]...)

df = DataFrame(index=indices,richness=richnesses, sampling_effort=sampeffort, mean_fnr=means)

CSV.write(joinpath("src", "artifacts", "1d_mangal_fnrs.csv"), df)

