
struct CooccurrenceMatrix{S,F}
    species::Vector{S}
    matrix::Matrix{F}
end

struct InteractionsMatrix{S,F}
    species::Vector{S}
    matrix::Matrix{F}
end

struct RelativeAbundance{S,F}
    species::Vector{S}
    relativeabundances::Vector{F}
end

num_partitions(ra::RelativeAbundance) = length(ra.species)
partition(ra::RelativeAbundance, i) = 
    RelativeAbundance([ra.species[i]], [ra.relativeabundances[i]])

function sample(ra::RelativeAbundance, n)
    K = num_partitions(ra)
    K == 1 ? ra.species[begin][rand(Categorical(ra.relativeabundances[begin]), n)] : sample.([partition(ra, i) for i in 1:K],n) 
end 


# sample from relative abundance
struct CooccurrenceSampler
    relativeabundance::RelativeAbundance
end
function sample(cs::CooccurrenceSampler, n)
    species = cs.relativeabundance.species
    cooc_sample = sample(cs.relativeabundance, n)
    mat = length(species) > 1 ? zeros(Bool,size(species[1])[1], size(species[2])[1]) : zeros(Bool,size(species[1]), size(species[1]))

    if length(species) == 1
        for (i,s) in enumerate(species[1]), (j,t) in enumerate(species[1])
            mat[i,j] = s ∈ cooc_sample && t ∈ cooc_sample[1]
        end
    else
        for (i,s) in enumerate(species[1]), (j,t) in enumerate(species[2])
            mat[i,j] = s ∈ cooc_sample[1] && t ∈ cooc_sample[2]
        end
    end

    CooccurrenceMatrix(cs.relativeabundance.species, mat)
end 


# compute FNR stats from Coocc sample
function FNR(cooccurrence_sample, C, M)
    cooc = findall(I->C[I] > 0, CartesianIndices(adjacency(M)))
    sampled, real = cooccurrence_sample.matrix, adjacency(M)

    fn = length(findall(I->real[I] && !sampled[I], CartesianIndices(sampled)))
    fn/prod(size(real[cooc]))
end 


