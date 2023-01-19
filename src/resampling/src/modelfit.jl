using EcologicalNetworks
using Plots
using Random
using EcologicalNetworks
using MultivariateStats
using ProgressMeter
using Statistics
using DataFrames
using MLJ
using StatsBase
using CSV

include("computemeasures.jl")

ids = map(i -> i.ID, filter(i -> contains(i.Reference, "Hadfield"), web_of_life()))
B = convert.(BipartiteNetwork, web_of_life.(ids))
M = reduce(∪, B)

function load_features()


    N = convert(UnipartiteNetwork, M)
    K = EcologicalNetworks.mirror(N)

    pc = fit(PPCA, Float64.(Array(K.edges)))
    pr = MultivariateStats.transform(pc, Float64.(Array(K.edges)))

    nf = 15

    cooc = zeros(Bool, prod(size(M)))                  # H/P co-occurrence?
    labels = zeros(Bool, prod(size(M)))                # H/P interaction?
    features = zeros(Float64, (2 * nf, prod(size(M))))   # H/P latent traits

    # We then move through the metaweb step by step
    cursor = 0
    for i in species(M; dims=1), j in species(M; dims=2)
        cursor += 1
        # Interaction in the metaweb?
        labels[cursor] = M[i, j]
        # Values in the PCA space
        p_i = findfirst(i .== species(N))
        p_j = findfirst(j .== species(N))
        features[1:nf, cursor] .= pr[1:nf, p_i]
        features[(nf + 1):end, cursor] .= pr[1:nf, p_j]
        # Co-occurrence?
        for b in B
            if i in species(b)
                if j in species(b)
                    cooc[cursor] = true
                end
            end
        end
    end

    # We only work on the species pairs that DO co-occur
    x = Float64.(copy(features[:, :])) # Latent traits
    y = labels

    x,y
end

function feature_dataframe(featuresmatrix, labels)
    total_nf, nlabels = size(featuresmatrix)
    nf_per_species = Int64(total_nf/2)

    df = DataFrame(label=labels)
    for i in 1:nf_per_species
        df[!, Symbol("X$i")] = zeros(Float64, nlabels)
    end
    for i in 1:nf_per_species
        df[!, Symbol("Y$i")] = zeros(Float64, nlabels)
    end

    for r in 1:nlabels
        df[r,2:end] .= featuresmatrix[:,r]
    end

    df
end


featdf = feature_dataframe(load_features()...)
y, X = unpack(featdf, ==(:label); rng = 123)
y = coerce(y, Multiclass{2})

# models(matching(X,y))


BRT = @load EvoTreeClassifier pkg = EvoTrees
RandomForest = @load RandomForestClassifier pkg = DecisionTree verbosity = 0
DecisionTree = @load DecisionTreeClassifier pkg = DecisionTree verbosity = 0
Logistic = @load LogisticClassifier pkg = ScikitLearn



# Ensemble model with different true/falses when negative sampling
# 
function get_balanced_test_train_sets(X,y, n=200; testratio=0.2)
    num_obs = length(y)

    int_y = [i.ref == 2 for i in y]

    training_size = Int64(floor(num_obs*(1-testratio)))
    Itrain = sample(1:num_obs, training_size, replace=false)
    Itest = (1:num_obs)[map(i->i ∉ Itrain, 1:num_obs)]
    train_positives = findall(int_y[Itrain])
    possible_negatives = Itrain[[!int_y[i]  for i in Itrain]]

    batches = []
    for i in 1:n
        train_negatives = sample(possible_negatives, length(train_positives), replace=false)
        balanced_Itrain = shuffle!([train_negatives..., train_positives...])
        push!(batches, balanced_Itrain)
    end 
    
    Itrain, Itest, batches
end 


#=Itrain, Itest, batches = get_balanced_test_train_sets(X,y, 100)

machs = Dict()

@showprogress for (i,batch_Itrain) in enumerate(batches)
    mach = machine(RandomForest(), X,y)
    fit!(mach, rows=batch_Itrain, verbosity=0)
    machs[Symbol("learner_$i")] = mach
end 


preds = [[i.prob_given_ref[2] for i in MLJ.predict(machs[k], X[Itest,:])] for k in keys(machs)]

ensemble_predict = [mean([p[i] for p in preds]) for i in eachindex(preds[1])]
true_y = [i.ref == 2 for i in y[Itest]]
computemeasures(true_y,ensemble_predict)

    

MLJ.fit!(RandomForest(),X,y, resampling=(Itrain,Itest), measures=[AUC])

histogram(ensemble_predict[findall(x->x>0,ensemble_predict)])
=#


#include("_computemeasures.jl")


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


#=

    uniqueXs = unique([x[1] for x in nonzeroC])
    uniqueYs = unique([x[2] for x in nonzeroC])
    speciesT, speciesB = M.T[uniqueXs], M.B[uniqueYs]

    finalC, finalN = zeros(length(speciesT), length(speciesB)), zeros(length(speciesT), length(speciesB))

    @info length(uniqueXs), length(uniqueYs)

    for I in nonzeroC
        sp1, sp2 = species(M,dims=1)[I[1]], species(M,dims=2)[I[2]]
        i,j = findfirst(x->x==sp1, speciesT), findfirst(x->x==sp2, speciesB)
        
        finalC[i,j] = C[I[1],I[2]]
        finalN[i,j] = N[I[1], I[2]]
    end

    return finalC, finalN =#
    return C, N
end


function get_predictions()
    C, N = get_metadata()
    nonzeroC = findall(x->x>0, C)


    mach = machine(BRT(), X,y)
    evaluate!(mach, resampling=CV(), measure=[auc], verbosity=1)
    MLJ.fit!(mach)

    F = reshape([x.prob_given_ref[2] for x in MLJ.predict(mach)], size(C))

    RA_T, RA_B = zeros(length(species(M, dims=1))), zeros(length(species(M, dims=2)))

    RA_T = [sum(M[i,:]) for i in 1:length(species(M,dims=1))]
    RA_B = [sum(M[:,i]) for i in 1:length(species(M,dims=2))]

    return C,N,F, (RA_T./sum(RA_T),RA_B./sum(RA_B))
end 

