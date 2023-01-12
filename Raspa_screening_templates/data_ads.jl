#!/bin/env julia

using CSV
using DataFrames

for file in readdir(@__DIR__, join=true)
    if startswith(basename(file), "DATAloading")
        rm(file)
    end
end

function getfirst(f, l)
    for x in l
        f(x) && return x
    end
    nothing
end

function seekgrep!(elines, T, pattern; stopat=nothing, antiword=nothing)
    U = Union{Nothing, pattern isa Regex ? T : eltype(elines)}
    for l in elines
        if pattern isa Regex
            m = match(pattern, l)
            if m isa RegexMatch
                x = only(m.captures)
                return Some{U}(T <: AbstractString ? x : parse(T, x))
            end
        elseif occursin(pattern, l)
            return Some{U}(l)
        end
        !(antiword isa Nothing) && occursin(antiword, l) && return Some{U}(nothing)
        !(stopat isa Nothing) && occursin(stopat, l) && return nothing
    end
    nothing
end

function allocurrences(T, file, pattern; stopat=nothing, antiword=nothing, cyclekeep=[true])
    ret = T[]
    elines = eachline(file)
    n = length(cyclekeep)
    k = count(cyclekeep)
    j = 1
    x = seekgrep!(elines, T, pattern; stopat, antiword)
    while !(x isa Nothing)
        val = something(x)
        if val isa Nothing
            resize!(ret, max(0, length(ret) - k))
        else
            cyclekeep[mod1(j, n)] && push!(ret, val)
            j += 1
        end
        x = seekgrep!(elines, T, pattern; stopat, antiword)
    end
    return ret
end

const DATAloading = joinpath(@__DIR__, "DATAloading")
function makeDATAloading(previous)
    @assert !isfile(DATAloading)
    prepends = Dict{Tuple{Int,String},SubString{String}}()
    if !isnothing(previous)
        old = joinpath(previous, "DATAloading")
        @assert isfile(old)
        for l in eachline(old)
            splits = split(l, ','; limit=3)
            prepends[(parse(Int,splits[1]),splits[2])] = splits[3]
        end
    end

    for system in readdir(joinpath(@__DIR__, "Output"), join=false)
        files = readdir(joinpath(@__DIR__, "Output", system), join=true)
        for file in files
            TEMP = something(seekgrep!(eachline(file), Int, r"^External temperature:\s*([^\s]+)"))
            splits = split(basename(file), '_')
            bbase = join(@view(splits[2:end-1]), '_')
            pressure = parse(Float64, first(splitext(splits[end])))
            movietemp = joinpath(@__DIR__, "Movies", system)
            moviefile = getfirst(readdir(movietemp, join=false)) do name
                s = split(name, '_')
                s[1] == "Movie" || return false
                s[end] == "frameworks.pdb" || return false
                join(@view(s[2:end-2]), '_') == bbase || return false
                parse(Float64, s[end-1]) == pressure || return false
                return true
            end
            species = allocurrences(String, file, r"Component [0-9]+ \[.*\] \(([^\s]+)"; stopat="Framework Status")
            cyclekeep = species .== "Adsorbate"
            loadings = allocurrences(Float64, file, r"\s*absolute adsorption\:.*\s+([^\s]+)\s+\(avg.\s+[^\s]+\) \[mol\/uc\]"; antiword="ENERGY DRIFT", cyclekeep)
            open(DATAloading, "a") do f
                key = (TEMP, basename(file))
                print(f, TEMP, ',', key[2], ',')
                pre = get(prepends, key, nothing)
                if !isnothing(pre)
                    print(f, pre, ',')
                end
                join(f, loadings, ',')
                println(f)
            end
        end
    end
end
makeDATAloading(isempty(ARGS) ? nothing : first(ARGS))

const gd = groupby(CSV.read(DATAloading, DataFrame; header=false), 1)
const nums = combine(gd, 1=>first; renamecols=false)[!, 1]
for (df, T) in zip(gd, nums)
    titles = [first(splitext(last(split(x, '_'; limit=2)))) for x in @view(df[:,2])]
    perm = rename!(permutedims(@view df[:, 3:end]), titles)
    CSV.write(joinpath(@__DIR__, "DATAloading_$T.csv"), perm)
end
rm(DATAloading)
