"""
    This code is primanry adapted from 
    https://github.com/PoisotLab/JuliaEcoNetworksValuePack/blob/v0.0.1/paper/main.jl
    in order to query Mangal and save id's of networks based on network type 
    (foodweb,parasitism,mutalism,misc) in files in the artifacts directory.  
"""


using Mangal
using DataFrames
using ProgressMeter
using SparseArrays
using CSV

mangal_networks = nothing 


function querymangal()
    global mangal_networks

    if isnothing(mangal_networks)
        numnetworks = count(MangalNetwork)
        mangal_networks = DataFrame(zeros(Int32, numnetworks,7), [:id, :S, :L, :pred, :herb, :mutu, :para])
        
        count_per_page = 100
        number_of_pages = convert(Int, ceil(numnetworks/count_per_page))
        cursor = 1
        @showprogress for page in 1:number_of_pages
            networks_in_page = Mangal.networks("count" => count_per_page, "page" => page-1)
            for current_network in networks_in_page
                S = count(MangalNode, current_network)
                L = count(MangalInteraction, current_network)
                pred = count(MangalInteraction, current_network, "type" => "predation")
                herb = count(MangalInteraction, current_network, "type" => "herbivory")
                mutu = count(MangalInteraction, current_network, "type" => "mutualism")
                para = count(MangalInteraction, current_network, "type" => "parasitism")
                mangal_networks[cursor,:] .= (current_network.id, S, L, pred, herb, mutu, para)
                cursor = cursor + 1
            end
        end
    end
    CSV.write(joinpath("src", "artifacts", "allmangal", "mangal.csv"), mangal_networks)

    # remove networks with less than 5 links
    mangal_networks2 = mangal_networks[mangal_networks.L .> 4, :]

    # number of interactions in food webs (sum of predation and herbivory interactions)
    mangal_networks.foodweb = mangal_networks.pred .+ mangal_networks.herb

    # number of other types of interactions
    mangal_networks.other = mangal_networks.L .- mangal_networks.foodweb .- mangal_networks.mutu .- mangal_networks.para

    # remove networks with less than 5 links
    mangal_networks2 = mangal_networks[mangal_networks.L .> 4, :]

    number_interactions_max_type = maximum.(eachrow(mangal_networks2[:,[:mutu, :para,:foodweb,:other]]))
    foodwebs = mangal_networks2[mangal_networks2[!, :foodweb] .== number_interactions_max_type ,:]
    parasitism_webs = mangal_networks2[mangal_networks2[!, :para] .== number_interactions_max_type ,:]
    mutualism_webs = mangal_networks2[mangal_networks2[!, :mutu] .== number_interactions_max_type ,:]
    other_webs = mangal_networks2[mangal_networks2[!, :other] .== number_interactions_max_type ,:]


    CSV.write(joinpath("src", "artifacts", "allmangal", "foodwebs.csv"), foodwebs)
    CSV.write(joinpath("src", "artifacts", "allmangal", "parasite.csv"), parasitism_webs)
    CSV.write(joinpath("src", "artifacts", "allmangal", "mutualist.csv"), mutualism_webs)
    CSV.write(joinpath("src", "artifacts", "allmangal", "misc.csv"), other_webs)

end


function writemangaldata()
    fw_df = CSV.read(joinpath("src", "artifacts", "allmangal", "foodwebs.csv"),DataFrame)
    para_df = CSV.read(joinpath("src", "artifacts", "allmangal", "parasite.csv"), DataFrame)
    mutu_df = CSV.read(joinpath("src", "artifacts", "allmangal", "mutualist.csv"), DataFrame)
    misc_df = CSV.read(joinpath("src", "artifacts", "allmangal", "misc.csv"), DataFrame)
    
    fw = network.(fw_df.id)
    para = network.(para_df.id)
    mutu = network.(mutu_df.id)
    misc =  network.(misc_df.id)

    writeedgelists(fw, "foodweb")
    writeedgelists(para, "parasite")
    writeedgelists(mutu, "mutualist")
    writeedgelists(misc, "misc")

end


function EcologicalNetworks.convert(::Type{EcologicalNetworks.UnipartiteNetwork}, interac::Vector{MangalInteraction})

    all_object_nodes = MangalNode[]

    for i in interac
        append!(all_object_nodes, [i.from, i.to])
    end

    object_nodes = unique(all_object_nodes)
    S = length(object_nodes)
    A = zeros(Bool, (S,S))
    N = EcologicalNetworks.UnipartiteNetwork(A, object_nodes)
    for i in interac
        if !ismissing(i.strength) && i.strength != 0
            N[i.from, i.to] = true
        end
    end

    return N
end

function writeedgelists(thewebs, typepath)
    for w in thewebs
        @info w 
        thismat = convert(UnipartiteNetwork,w)
        pathname = joinpath("src", "artifacts", "allmangal", typepath, "$(w.name).csv" ) 
        Is,Js,trash = findnz(thismat.edges)
        CSV.write(pathname,  DataFrame([Is, Js], [:i,:j]))
    end
end

function read_edgelist(filepath, type)
    df = CSV.read(filepath, DataFrame)
    nrow(df) > 0 || return nothing

    sz = max(max(df.i...), max(df.j...))
    A = zeros(sz,sz)    

    
    
    for r in 1:nrow(df)
        i = df[r, :i]
        j = df[r, :j]
        A[i,j] = 1
    end
    return type(Bool.(A))
end


function read_edgelists(dir, type::Type{T}) where {T<: AbstractEcologicalNetwork}
    filenames = filter(x->endswith(x, ".csv"), readdir(dir))
    nets = []
    for file in filenames 
        mat = read_edgelist(joinpath(dir, file), type)
        !isnothing(mat) && push!(nets, mat)
    end
    return nets
end


function mangaldata()
    fw = read_edgelists(joinpath("src", "artifacts", "allmangal", "foodweb"), UnipartiteNetwork)
    para = read_edgelists(joinpath("src", "artifacts", "allmangal", "parasite"), BipartiteNetwork)
    mutu = read_edgelists(joinpath("src", "artifacts", "allmangal", "mutualist"),  BipartiteNetwork)
    misc = read_edgelists(joinpath("src", "artifacts", "allmangal", "misc"), UnipartiteNetwork)

    fw, para, mutu, misc
end
