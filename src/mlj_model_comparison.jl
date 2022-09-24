using MLJ
using MLJFlux
using EcologicalNetworks
using DataFrames
using MultivariateStats

include("_computemeasures.jl")

function load_features()
    # Get the data from the Web of Life database, and merge them into a metaweb M
    ids = map(i -> i.ID, filter(i -> contains(i.Reference, "Hadfield"), web_of_life()))
    B = convert.(BipartiteNetwork, web_of_life.(ids))
    M = reduce(∪, B)

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
    kept = findall(cooc)
    x = Float32.(copy(features[:, kept])) # Latent traits
    y = labels[kept]     # Interaction bit

    x,y
end

function feature_dataframe(featuresmatrix, labels)
    total_nf, nlabels = size(featuresmatrix)
    nf_per_species = Int32(total_nf/2)

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

function split(X, y, prop = 0.8)
    training_size = convert(Int64, floor(nrow(X) * 0.8))
    train = sort(sample(1:nrow(X), training_size; replace = false))
    test = filter(i -> !(i in train), 1:nrow(X))
    train, test
end

featdf = feature_dataframe(load_features()...)

y, X = unpack(featdf, ==(:label); rng = 123)
y = coerce(y, Multiclass{2})


Tree = @load EvoTreeClassifier pkg = EvoTrees
RandomForest = @load RandomForestClassifier pkg = DecisionTree verbosity = 0
DecisionTree = @load DecisionTreeClassifier pkg = DecisionTree verbosity = 0
Ada = @load AdaBoostClassifier pkg = ScikitLearn verbosity = 0
Logistic = @load LogisticClassifier pkg = ScikitLearn
KNN = @load KNNClassifier 


# Does amount of bagging matter? 
brt = EnsembleModel(model = Tree(), n = 100)
rf = EnsembleModel(model = RandomForest(), n = 100)
dt = EnsembleModel(model = DecisionTree(), n = 100)
ada = EnsembleModel(model = Ada(), n = 100)
logit = Logistic()
knn = EnsembleModel(model = KNN(), n = 100)

MLJ.evaluate(rf, X,y, resampling=CV(shuffle=true), measures=[kappa])



mach = machine(dt, X, y)
fit!(mach, rows = Itrain)
pred = MLJ.predict(mach, rows = Itest)
prediction = pred.prob_given_ref[2]
obs = Bool.(y[Itest] .== true)
prediction, computemeasures(obs, prediction)