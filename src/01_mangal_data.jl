using Mangal
using EcologicalNetworks
using ProgressMeter
using Plots
using DataFrames
using CSV
using BSON
using StatsBase
using Random
Random.seed!(42)

sets = Dict(
    "kolpelke_et_al_2017" => "Kopelke et al. (2017)",
    "hadfield_2014" => "Hadfield et al (2014)",
    "havens_1992" => "Havens (1992)", 
    "ponisio_2017" => "Ponisio et al. (2017)",
    "RMBL_pollination" => "CaraDonna et al. (2017)",        
    "closs_1994" => "Closs & Lake (1994)"
)



# I = findall(x->x>5, count.(MangalNetwork, datasets()))
data = Dict([s => dataset(s) for (s,n) in sets])
allnets = Dict([n=>networks(d) for (n,d) in data])


# Estimate how long htis takes, but it might involve IO? 
# Just save a BSON is the best option 
nets = Dict([n=>convert.(UnipartiteNetwork, d) for (n,d) in allnets])


# HUGE NOTE
# Each network in RMBL needs to be aggregated before it is used to build the
# metaweb because they contain sepeare 2-node networks for each interaction 
# This should be possible by unioning all the interactions, but it might take a while


# Convert each list of networks to one with the same taxonomic backbone 
function convert_taxa(net)
    name, As = net

    thesenets = []

    namesdf = CSV.read(joinpath("src", "artifacts", "species_translations", "$name.csv"), DataFrame)
    for (netnum,A) in enumerate(As)
        spnames = String[]
        Ind = []
        thisnet_df = filter(r->r.network_number==netnum, namesdf)            

        num_dropped = length(findall(r-> ismissing(r.DROP) || r.DROP==true, eachrow(thisnet_df)))
        filter!(r-> ismissing(r.DROP) || r.DROP == false, thisnet_df)

        for (i,r) in enumerate(eachrow(thisnet_df))
            mang_taxa  = !(ismissing(r.mangal_taxon_name)) ? r.mangal_taxon_name : nothing 
            mang_species = !(ismissing(r.mangal_species_name)) ? r.mangal_species_name : nothing

            if ismissing(r.DROP) || r.DROP == false
                i = findall(n -> n.name == mang_species || mang_taxa == (!ismissing(n.taxon) ? n.taxon.name : ""), species(A))
                push!(Ind, i[1])
                push!(spnames, r.CONSENSUS_NAME)
            end 
        end 


        rows_with_duplicated_consensus_names = findall(nonunique(thisnet_df,:CONSENSUS_NAME)) 
            
        dups = unique([findall(x->x.CONSENSUS_NAME == n, eachrow(thisnet_df)) for n in thisnet_df.CONSENSUS_NAME[rows_with_duplicated_consensus_names]])

        
        function findnodes(A, species_consensus)
            matching_df_rows = filter(r->r.CONSENSUS_NAME == species_consensus, thisnet_df)  
            matching_nodes = []
            for s in species(A)
                mang_taxa  = !ismissing(s.taxon) ? s.taxon.name : nothing 
                mang_species = s.name 
                
                for r in eachrow(matching_df_rows)
                    df_taxa  = !(ismissing(r.mangal_taxon_name)) ? r.mangal_taxon_name : nothing 
                    df_species = !(ismissing(r.mangal_species_name)) ? r.mangal_species_name : nothing
        
                    if (mang_taxa == df_taxa || mang_species == df_species)
                        push!(matching_nodes, s)
                    end
                end
            end
            matching_nodes
        end     

        if length(dups) > 0
            
            # Find all nodes in A that correspond to a unique species name
            spnames = unique(spnames)

            # Now return a dict that is species consensus name => list of nodes 
            d = Dict()
            for sp in spnames
                correspondingnodes = findnodes(A, sp)
                merge!(d, Dict(sp=>correspondingnodes))
            end

            adj = zeros(Bool,length(spnames),length(spnames))

            findkey(x) = begin
                for (k,v) in d
                    x ∈ v && return k
                end
            end

            for (i,sp) in enumerate(spnames)
                this_species_nodes = d[sp]
                for n in this_species_nodes
                    left, right =  A[:,n], A[n,:]
                    if length(left) > 0
                        lefttargs = findkey.(left)
                        for targ in lefttargs  
                            if !isnothing(targ)
                                targindex = findall(x->x==targ, spnames)[1]
                                adj[i,targindex] = 1
                            end
                        end
                    end 

                    if length(right) > 0 
                        righttargs = findkey.(right)
                        for targ in righttargs  
                            if !isnothing(targ)
                                targindex = findall(x->x==targ, spnames)[1]
                                adj[targindex,i] =1  
                            end
                        end
                    end 
                end
            end
            nonzero_species_indices = []
            for (i,sp) in enumerate(spnames)
                if sum(adj[i,:]) > 0 || sum(adj[:,i]) > 0
                    push!(nonzero_species_indices, i)
                end
            end 
            thisnet = UnipartiteNetwork(adj[nonzero_species_indices,nonzero_species_indices], spnames[nonzero_species_indices])
            push!(thesenets, thisnet)    

        else
            thisnet = UnipartiteNetwork(adjacency(A)[Ind,Ind], spnames)
            push!(thesenets, thisnet)    
        end


      #  @info species(A)
      #  @info Ind
      #  @info adjacency(A)
    end
    thesenets
end 

unipartite_nets = Dict()
@time for x in nets
    @info x[1] 
    As = convert_taxa(x)
    merge!(unipartite_nets, Dict(x[1]=>As))
end 



merge!(unipartite_nets, Dict("nz_stream_foodweb"=>nz_stream_foodweb()))


final_nets = Dict()
for (k,v) in unipartite_nets
    try
        b = convert.(BipartiteNetwork, v) 
        merge!(final_nets, Dict(k=>b))
    catch
        merge!(final_nets, Dict(k=>v))
    end
end

final_nets

plot([heatmap(adjacency(reduce(∪,v)), colorbar=:none, title=n) for (n,v) in unipartite_nets]..., layout=(2,4), size=(1200, 900))


Base.in(spname::String,net::NT) where {NT <: UnipartiteNetwork } = spname in net.S
Base.in(spname::String, net::NT) where {NT <: BipartiteNetwork } = spname in net.T || spname in net.B

function checkproductunipartite(M, marginalcount, jointcount)
    metaweb = adjacency(M)
    diffs = []
    joints = []
    margs = []

    S = richness(M)

    jointsum = sum(jointcount)
    margsum = sum(marginalcount)

    for i in 1:S
        P_i = marginalcount[i]/margsum
        for j in 1:S
            if i > j && metaweb[i,j] == 1
                P_j = marginalcount[j]/ margsum
                P_ij = jointcount[i,j] / jointsum # joint probab of seeing i and j
                thisdiff = P_ij - P_i*P_j
                push!(diffs, thisdiff)
                push!(joints, P_ij)
                push!(margs, P_i*P_j)
            end 
        end
    end
    return diffs, joints, margs
end

function checkproductbipartite(M, hostcounts, paracounts, jointcount)
    metaweb = adjacency(M)
    diffs = []
    joints = []
    margs = []
    numpara, numhost = size(metaweb)

    for i in 1:numpara
        P_i = paracounts[i]/sum(paracounts)
        for j in 1:numhost
            if metaweb[i,j] == 1
                P_j = hostcounts[j] / sum(hostcounts)
                P_ij = jointcount[i,j] / sum(jointcount)
                thisdiff = P_ij - P_i*P_j
                push!(diffs, thisdiff)
                push!(joints, P_ij)
                push!(margs, P_i*P_j)
            end 
        end
    end


    return diffs, joints, margs
end


function unipartitesplit(As)
    M = reduce(∪, As)
    speciespool = M.S
    S = richness(M)
    marginalcount = zeros(S)
    jointcount = zeros(S,S)
    for A in As 
        for (i,speciesi) in enumerate(speciespool)
            for (j, speciesj) in enumerate(speciespool)
                if i > j
                    marginalcount[i] +=  speciesi in A  
                    marginalcount[j] +=  speciesj in A
                    jointcount[i,j] += (speciesi in A && speciesj in A)
                end 
            end
        end
    end

    return checkproductunipartite(M,marginalcount, jointcount)
end

function bipartitesplit(Bs)
    M = convert(BipartiteNetwork,reduce(∪, Bs))
    
    hostpool = M.B
    parapool = M.T
   
    host_marginalcount = zeros(length(hostpool))
    para_marginalcount = zeros(length(parapool))
    jointcount = zeros(length(parapool),length(hostpool), )

    @showprogress for (i,para) in enumerate(parapool)
        for (j,host) in enumerate(hostpool)
            for B in Bs 
                    para_marginalcount[i] +=  para in B
                    host_marginalcount[j] +=  host in B  
                    jointcount[i,j] += (host in B && para in B)
            end
        end
    end

    return checkproductbipartite(M, host_marginalcount, para_marginalcount, jointcount)
end


function makejointmargplots(As)
    if typeof(As[begin]) <: UnipartiteNetwork
        diffs, joints, margs = unipartitesplit(As);
        return (diffs, joints, margs)
    elseif typeof(As[begin]) <: BipartiteNetwork
        diffs, joints, margs = bipartitesplit(As);
        return (diffs, joints, margs)
    end
end 

diffs = Dict()
for (n,v) in final_nets
    diff, joints, margs = makejointmargplots(v)
    merge!(diffs, Dict(n=>diff))
end


plot([histogram(d, title=n, bins=30) for (n,d) in diffs] ..., layout=(2,4),size=(1200,900))
 

#=
    Info or write the information for each network 
=#

# Write diffs 

diffs_df = DataFrame(dataset=[], diff=[])
for (n,v) in final_nets
    thesediffs = diffs[n]
    for d in thesediffs
    push!(diffs_df.diff, d)
    push!(diffs_df.dataset,n)
    end
end

CSV.write(joinpath("src", "artifacts", "joint_marginal_diffs.csv"), diffs_df)

titles = Dict(
    "kolpelke_et_al_2017" => "Kopelke et al. (2017)",
    "hadfield_2014" => "Hadfield et al (2014)",
    "havens_1992" => "Havens (1992)", 
    "ponisio_2017" => "Ponisio et al. (2017)",
    "RMBL_pollination" => "CaraDonna et al. (2017)",        
    "closs_1994" => "Closs & Lake (1994)",
    "nz_stream_foodweb" => "Townsend & Thompson (2000)";
)




# Beta div
all_βwn = Dict()
all_βos = Dict()

for (n,v) in final_nets
    βwn_vec = []
    βos_vec = []

    for i in 1:length(v), j in i+1:length(v)
        push!(βwn_vec, KGL01(βwn(v[i], v[j])))
        push!(βos_vec, KGL01(βos(v[i], v[j])))
    end 

    merge!(all_βos, Dict(n=>βos_vec))
    merge!(all_βwn, Dict(n=>βwn_vec))
end


betaplts = []
for (n,_) in all_βwn
    alpha = 100. / length(all_βwn[n])^0.85
    thisplt = scatter(all_βwn[n], all_βos[n], xlabel="βwn", ylabel="βos",frame=:box, msw=0, title=titles[n], ma=alpha, label="", xlim=(1,2), ylim=(1,2), aspectratio=1)
    plot!([1,2],[1,2], lc=:grey, linestyle=:dash, label="")
    push!(betaplts, thisplt)
end
plot(betaplts..., size=(900,900), dpi=200)

β_df = DataFrame(dataset=[], βwn=[],βos=[])
for (n,v) in all_βwn
    these_βwn = all_βwn[n]
    these_βos = all_βos[n]

    for i in 1:length(these_βos)
        push!(β_df.βos, these_βos[i])
        push!(β_df.βwn, these_βwn[i])
        push!(β_df.dataset,n)
    end
end

CSV.write(joinpath("src", "artifacts", "beta_diversity.csv"), β_df)



# Write metadata 
metadata_df = DataFrame(
    dataset=[], 
    num_networks=[],
    metaweb_richness=[],
    mean_local_richness=[],
    metaweb_connectance=[],
    mean_local_connectance=[],
    mean_βos = [],
    mean_βwn = []
)

for (n,v) in final_nets
    push!(metadata_df.dataset, n)
    push!(metadata_df.num_networks, length(v))

    mw = reduce(∪, v)
    push!(metadata_df.metaweb_richness, richness(mw))
    push!(metadata_df.metaweb_connectance, connectance(mw))

    push!(metadata_df.mean_local_richness, mean([richness(a) for a in v ]))
    push!(metadata_df.mean_local_connectance, mean(filter(!isnan, [connectance(a) for a in v ])))

    meanβos = mean(filter(!isnan, all_βos[n]))
    meanβwn = mean(filter(!isnan, all_βwn[n]))

    push!(metadata_df.mean_βos, meanβos)
    push!(metadata_df.mean_βwn, meanβwn)
end


metadata_df

CSV.write(joinpath("src", "artifacts", "metadata.csv"), metadata_df)


# old code

# overload some functions 
