using MLJ
using MLJFlux
using EcologicalNetworks
using DataFrames
using MultivariateStats
using Plots

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


function add_false_negatives(y, Itrain, fnr)
    y_with_fn = similar(y)
    y_with_fn = zeros(Bool, length(y))

    for i in Itrain
        if y[i]==true && rand() < (1-fnr)
            y_with_fn[i] = 1
        else
            y_with_fn[i] = 0
        end
    end

    boolvec = [yi == true for yi in y]
    creal, cobs = sum(boolvec[Itrain])/length(boolvec[Itrain]), sum(y_with_fn[Itrain])/length(y_with_fn[Itrain])
    @info "Took C=$creal, added FNs to make it C=$cobs"
    coerce(y_with_fn, Multiclass{2})
end

function test(model, X, y, Itrain, Itest; fnr=0.1)

    y_with_fn = add_false_negatives(y, Itrain, fnr)

    mach = machine(model, X, y_with_fn)
    fit!(mach, rows = Itrain)
    pred = MLJ.predict(mach, rows = Itest)
    prediction = pred.prob_given_ref[2]
    obs = Bool.(y[Itest] .== true)
    prediction, computemeasures(obs, prediction)
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
brt = EnsembleModel(model = Tree(), n = 200)
rf = EnsembleModel(model = RandomForest(), n = 200)
dt = EnsembleModel(model = DecisionTree(), n = 200)
ada = EnsembleModel(model = Ada(), n = 100)
logit = Logistic()
knn = EnsembleModel(model = KNN(), n = 100)

MLJ.evaluate(rf, X,y, resampling=CV(shuffle=true), measures=[kappa])



modellist = [brt, rf, dt, logit]

df = DataFrame(model=[], added_fnr=[], rocauc=[], prauc=[])


for fnr in 0.0:0.05:0.95
    @info "FNR: $fnr"
    Itrain, Itest = split(X,y)
    for m in modellist 
        thispred, thismetrics = test(m, X,y, Itrain, Itest; fnr=fnr)

        push!(df.model, m)
        push!(df.added_fnr, fnr)
        push!(df.prauc, thismetrics[:prauc])
        push!(df.rocauc, thismetrics[:rocauc])

    end 
end


namedict= Dict(
    brt => "BRT",
    rf => "RF",
    dt => "DT",
    logit => "Logistic",
    knn => "KNN"
)


 using CSV
CSV.write("model_comparison.csv", df)

plt = plot([0, 1, 1, 0], [0, 0, 0.5, 0.5], xminorticks=2,xguidefontsize=18, yguidefontsize=18, xtickfontsize=14, ytickfontsize=14, legend_font_pointsize=16,dpi=200, size=(1200, 900),aspectratio=1, st=:shape, c=:grey, frame=:box, xticks=0:0.1:0.9, yticks=0:0.1:1, alpha=0.2, lw=0.0,  legend=:bottomleft, lab="")

annotate!(0.5, 0.48, text("Worse than random", :black, :center, 12))
plot!([0, 1, 1, 0], [0.5, 0.5, 0.75, 0.75], st=:shape, c=:red, alpha=0.2, lw=0.0, lab="")
annotate!(0.5, 0.53, text("Close to random", :red, :center, 12))
plot!([0, 1, 1, 0], [0.75, 0.75, 0.9, 0.9], st=:shape, c=:orange, alpha=0.2, lw=0.0, lab="")
annotate!(0.5, 0.88, text("Fair classifier", :orange, :center, 12))
plot!([0, 1, 1, 0], [0.9, 0.9, 1.0, 1.0], st=:shape, c=:green, alpha=0.2, lw=0.0, lab="")
annotate!(0.5, 0.98, text("Excellent classifier", :green, :center, 12))
xaxis!("Added False Negative Rate", (0, 0.95))

yaxis!("ROC-AUC", (0.,1))

cols = [:dodgerblue, :seagreen4, :firebrick3, :mediumpurple3, :grey]
shapes = [:circle, :square, :diamond, :utriangle, :dtriangle]
msizes = [8, 8, 8, 8, 8]

for (i,m) in enumerate(modellist) 
    thisdf = filter(r->r.model==m, df)
    plot!(thisdf.added_fnr, thisdf.rocauc, lw=1.5, c=cols[i], label="")
    scatter!(thisdf.added_fnr, thisdf.rocauc, marker=shapes[i],ms=msizes[i], msw=1.5,msc=cols[i], mc=:white, label=namedict[m])
end

plt



plt2 = plot([0, 1, 1, 0], [0, 0, 0.5, 0.5], xguidefontsize=18, yguidefontsize=18, xtickfontsize=14, ytickfontsize=14,legend_font_pointsize=16,dpi=200, size=(1200, 900),aspectratio=1, st=:shape, c=:grey, frame=:box, xticks=0:0.1:0.9, yticks=0:0.1:1, alpha=0.2, lw=0.0,  legend=:bottomleft, lab="")
annotate!(0.5, 0.48, text("Worse than random", :black, :center, 12))
plot!([0, 1, 1, 0], [0.5, 0.5, 0.75, 0.75], st=:shape, c=:red, alpha=0.2, lw=0.0, lab="")
annotate!(0.5, 0.53, text("Close to random", :red, :center, 12))
plot!([0, 1, 1, 0], [0.75, 0.75, 0.9, 0.9], st=:shape, c=:orange, alpha=0.2, lw=0.0, lab="")
annotate!(0.5, 0.88, text("Fair predictor", :orange, :center, 12))
plot!([0, 1, 1, 0], [0.9, 0.9, 1.0, 1.0], st=:shape, c=:green, alpha=0.2, lw=0.0, lab="")
annotate!(0.5, 0.98, text("Excellent predictor", :green, :center, 12))
xaxis!("Added False Negative Rate", (0, 0.95))

for (i,m) in enumerate(modellist) 
    thisdf = filter(r->r.model==m, df)
    plot!(thisdf.added_fnr, thisdf.prauc, lw=1.5, c=cols[i], label="")
    scatter!(thisdf.added_fnr, thisdf.prauc, marker=shapes[i],ms=msizes[i], msw=1.5,msc=cols[i], mc=:white, label=namedict[m])
end
yaxis!("PR-AUC", (0, 1.))


plt2
using Plots.PlotMeasures: mm
plot(plt, plt2, dpi=200,size=(1500, 800), margin=10mm)
savefig("fig3.png")