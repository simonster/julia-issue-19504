module Sessions
using JSON
export Paradigm, Session, ChannelMetadata, ClusterMetadata, sessions,
       singlenev, nevstem, hascrash, hasrestart, isincluded, channels,
       channelmua, clusters, clustersua, filtereddatafile, metadatafile,
       filtereddata, events, mua, sua, matfile, nevfile, duration, firingrate,
       stationarity

const RECORDINGS_DIR = dirname(@__FILE__)
const DENOISED_DIR = "/mnt/data/patty_recordings_denoised"
const DATA_DIR = "/mnt/data/patty_recordings_data"
const IT_REGIONS = Set([:IT_VP])

immutable Paradigm
    nev::Vector{String}
    mat::Vector{String}
    paradigm_number::Vector{Int}
end

immutable Session
    date::String
    channels::Dict{Symbol,Vector{Int}}
    channel_status::Dict{Int,Symbol}
    paradigms::Dict{Symbol,Paradigm}
end

immutable ChannelMetadata
    date::String
    ch::Int
    region::Symbol
    status::Symbol
end

metadatafile(dir::String) = joinpath(RECORDINGS_DIR, dir, "metadata.json")

function Session(dir::String)
    metadata = JSON.parse(readstring(metadatafile(dir)))
    @assert metadata["date"] == dir
    channels = Dict{Symbol,Vector{Int}}([(Symbol(k), convert(Vector{Int}, v)::Vector{Int}) for (k, v) in metadata["channels"]])
    if haskey(metadata, "channel_status")
        channel_status = Dict{Int,Symbol}([(parse(Int, k), Symbol(v)) for (k, v) in metadata["channel_status"]])
    else
        channel_status = Dict{Int,Symbol}()
    end
    paradigms = Dict{Symbol,Paradigm}()
    for (k, v) in metadata["paradigms"]
        nev = v["nev"]
        !isa(nev, Vector) && (nev = [nev])
        mat = v["mat"]
        !isa(mat, Vector) && (mat = [mat])
        paradigm_number = v["paradigm_number"]
        !isa(paradigm_number, Vector) && (paradigm_number = [paradigm_number])
        paradigms[Symbol(k)] = Paradigm(nev, mat, paradigm_number)
    end
    Session(dir, channels, channel_status, paradigms)::Session
end

function sessions(includedonly::Bool=true)
    [Session("2016-03-31")]
end

function nevstem(s::Session, paradigm::Paradigm)
    nev = paradigm.nev
    length(nev) > 1 && error("multiple nev files for session $s.date")
    return nev[1][1:rsearch(nev[1], '.')-1]
end
nevstem(s::Session, paradigm::Symbol) = nevstem(s, s.paradigms[paradigm])

function nevstem(s::Session)
    nev = first(s.paradigms)[2].nev
    length(nev) > 1 && throw(ArgumentError("multiple nev files"))e
    nev = nev[1]
    for p in values(s.paradigms)
        if length(p.nev) > 1 || p.nev[1] != nev
            throw(ArgumentError("multiple nev files"))
        end
    end
    return nev[1:rsearch(nev, '.')-1]
end

function channels(s::Session, region::Symbol, chtype::Symbol=:all)
    !haskey(s.channels, region) && return Int[]
    chtype == :all && return [ChannelMetadata(s, ch, region) for ch in s.channels[region]]
    ret = ChannelMetadata[]
    for ch in s.channels[region]
        if !haskey(s.channel_status, ch)
            warn("no channel_status entry for $(s.date) ch$ch")
        elseif s.channel_status[ch] == chtype
            push!(ret, ChannelMetadata(s.date, ch, region, chtype))
        end
    end
    ret
end

function metadatafile(s::Session, paradigm::Symbol, ch::ChannelMetadata)
    joinpath(dirname(@__FILE__, "metadata.jld"))
end

function filtereddata(s::Session, paradigm::Symbol, ch::ChannelMetadata)
    if ch.region == :IT_VP
        stem = nevstem(s, paradigm)
        datadir = joinpath(DENOISED_DIR, string(stem, "_it"))
        f = joinpath(datadir, "$(stem)_$(suffix)_ica.dat")
        fh = open(f)
        try
            data = Mmap.mmap(fh, Matrix{Int16}, (length(s.channels[:IT_VP]), filesize(fh)÷2÷length(s.channels[:IT_VP])))
            return data[findfirst(s.channels[:IT_VP], ch.ch), :]
        finally
            close(fh)
        end
    else
        fh = open(filtereddatafile(s, paradigm, ch))
        try
            return Mmap.mmap(fh, Vector{Int16}, filesize(fh)÷2)
        finally
            close(fh)
        end
    end
end

function ChannelMetadata(s::Session, ch::Int, region::Symbol)
    if haskey(s.channel_status, ch)
        ChannelMetadata(s.date, ch, region, s.channel_status[ch])
    else
        ChannelMetadata(s.date, ch, region, :unknown)
    end
end
end