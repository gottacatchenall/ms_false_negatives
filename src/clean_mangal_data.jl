using Mangal
using EcologicalNetworks
using ProgressMeter
using Plots
using DataFrames
using CSV
using BSON
datasets()


sets = Dict(
    #"hawkins_goeden_1984" => "Hawkins & Goeden (1984)",      # not enough samples
    "kolpelke_et_al_2017" => "Kopelke et al. (2017)",
    "hadfield_2014" => "Hafield et al (2014)",
    "havens_1992" => "Havens (1992)", 
    "ponisio_2017" => "Ponisio et al. (2017)",
    "RMBL_pollination" => "CaraDonna et al. (2015)",        
    "closs_1994" => "Closs & Lake (1994)",
    #"parker_huryn_2006" => "Parker & Huryn (2006)",         
    #"ricciardi_2010" => "Ricciardi & MacIsaac (2010)",       # not many species 
        #"fryer_1959" => "Fryer (1959)",                         # Not usable
    #"primack_1983" => "Primack (1983)",                     # Only 3 loctations, toss
    #"ruzicka_2012" => "Ruzicka et al. (2012)",              # This isn't usable
)



# I = findall(x->x>5, count.(MangalNetwork, datasets()))
data = Dict([s => dataset(s) for (s,n) in sets])
allnets = Dict([n=>networks(d) for (n,d) in data])


# Estimate how long htis takes, but it might involve IO? 
# Just save a BSON is the best option 
nets = Dict([n=>convert.(UnipartiteNetwork, d) for (n,d) in allnets])

bson("mangal_networks.bson")


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


        # For species that are duplicated (meaning they have diff mangal
        # taxa/species names but the same CONSENSUS_NAME because they are
        # the same species), we need to combine their rows in the adjacency
        # matrix. 


        # A[consensus,:] = A[:S1,:] .+ A[:S2,:] 
        # and 
        # A[:,consensus] = A[:,:S1] .+ A[:,:S2]


        rows_with_duplicated_consensus_names = findall(nonunique(thisnet_df,:CONSENSUS_NAME)) 
        
        # now find the first row each duplicated_row corresponds to 
        
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

species(unipartite_nets["RMBL_pollination"][1])
adjacency(unipartite_nets["RMBL_pollination"][1])

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

plot([spy(adjacency(reduce(∪,v)), title=n) for (n,v) in unipartite_nets]..., layout=(2,4), size=(1200, 900))

diffs = Dict()
for (n,v) in unipartite_nets
    diff, joints, margs = makejointmargplots(v)
    @info n, length(diff)
    merge!(diffs, Dict(n=>diff))
end

diffs["parker_huryn_2006"]




plot([histogram(d, title=n, bins=30) for (n,d) in diffs] ..., layout=(2,4),size=(1200,900))
 





# old code

# overload some functions 

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
            @info j
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