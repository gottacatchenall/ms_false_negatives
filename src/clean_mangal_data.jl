using Mangal
using EcologicalNetworks
using ProgressMeter
using Plots
using DataFrames
using CSV

datasets()


sets = Dict(
    "hawkins_goeden_1984" => "Hawkins & Goeden (1984)",     
    "kolpelke_et_al_2017" => "Kopelke et al. (2017)",
    "hadfield_2014" => "Hafield et al (2014)",
    "havens_1992" => "Havens (1992)", 
    "ricciardi_2010" => "Ricciardi & MacIsaac (2010)",
    "ponisio_2017" => "Ponisio et al. (2017)",
    "ruzicka_2012" => "Ruzicka et al. (2012)",
    "RMBL_pollination" => "CaraDonna et al. (2015)", # this is not fully connected and thus breaks 
    "fryer_1959" => "Fryer (1959)",
    "closs_1994" => "Closs & Lake (1994)",
    "primack_1983" => "Primack (1983)",
    "parker_huryn_2006" => "Parker & Huryn (2006)"
)



# I = findall(x->x>5, count.(MangalNetwork, datasets()))
datanames = collect(keys(sets))
data = [dataset(d) for d in datanames]

allnets = [networks(d) for d in data]

allnets

nets = [convert.(UnipartiteNetwork, d) for d in allnets]

splist = vcat([i for i in species.(nets[1])]...)
unique(splist)


nets[1]

# Put them in a CSV for name cleaning 


for i in eachindex(nets)
    df = DataFrame(network_number=[], mangal_taxon_name=[], mangal_species_name=[],mangal_node_id=[])
    for (i,thisnet) in enumerate(nets[i])
        thisspecies = species(thisnet)
        for s in thisspecies
            push!(df.network_number, i)
            push!(df.mangal_species_name, s.name)
            push!(df.mangal_taxon_name, ismissing(s.taxon) ? missing : s.taxon.name)
            push!(df.mangal_node_id, s.id)
        end
    end 

    CSV.write("$(datanames[i]).csv", df)
end



newnets = []

for netvec in nets
    adjmats = adjacency.(netvec)
    thesenets = UnipartiteNetwork[]
    for i in eachindex(netvec)
        names = string.([ismissing(s.taxon) ? s.name : s.taxon.id for s in species(netvec[i])])
        push!(thesenets, UnipartiteNetwork(adjmats[i], names))
    end
    
    push!(newnets, thesenets)
end


newnets

push!(newnets, UnipartiteNetwork.(adjacency.(nz_stream_foodweb())))

newnets

pltdata = makejointmargplots.(newnets)

filter!((x)->!isempty(x[1]), pltdata)

plt = plot(aspectratio=1)
diffs = [(y .- x) for (x,y) in pltdata]

plot(histogram.(diffs)...)
 

function makehists(diffs)
    hists = []
    for (i,diff) in enumerate(diffs)
        mi,mx = min(diff...), max(diff...)
        if mx > mi
            @info typeof(diff)
            binsize = (mx-mi)/100
            plot()
            push!(hists, histogram(diff, bins=mi:binsize:mx, xlim=(min(mi,0), mx)))
        else
            @info "There is only one value for diffs at index $i"
        end
    end
    hists 
end

plot(makehists(diffs)..., size=(1200, 800))





metawebs = []
for (i,thesenets) in enumerate(newnets)
    try 
        a = convert.(BipartiteNetwork, thesenets)
        push!(metawebs, reduce(∪, convert.(BipartiteNetwork, thesenets)))
    catch 
        @info "net looks unipartite..."
        @info typeof(thesenets)
        push!(metawebs, reduce(∪, thesenets))
    end
end


metawebs

UnipartiteNetwork.(uni_nets)

plot(heatmap.(adjacency.(metawebs))..., colorbar=:none)





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
        margA, margB, jointAB = unipartitesplit(As);
        return (margA, margB, jointAB)
    elseif typeof(As[begin]) <: BipartiteNetwork
        margA, margB, jointAB = bipartitesplit(As);
        return (margA, margB, jointAB)
    end
end 