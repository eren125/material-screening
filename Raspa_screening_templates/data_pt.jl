#!/bin/env julia

using CSV
using DataFrames

for file in readdir(@__DIR__, join=true)
    if startswith(basename(file), "DATAenergy")
        rm(file)
    end
end

function getfirst(f, l)
    for x in l
        f(x) && return x
    end
    nothing
end

function seekgrep!(elines, T, pattern)
    for l in elines
        if pattern isa Regex
            m = match(pattern, l)
            m isa RegexMatch && return parse(T, only(m.captures))
        elseif occursin(pattern, l)
            return l
        end
    end
    nothing
end

function allocurrences(T, file, pattern)
    ret = T[]
    elines = eachline(file)
    x = seekgrep!(elines, T, pattern)
    while !(x isa Nothing)
        push!(ret, x)
        x = seekgrep!(elines, T, pattern)
    end
    return ret
end

function _moviefinder(bbase, pressure, suff)
    suffix::String = suff*".pdb"
    k::Int = count(==('_'), suff)
    function (name)
        s = split(name, '_')
        s[1] == "Movie" || return false
        endswith(name, suffix) || return false
        join(@view(s[2:end-2-k]), '_') == bbase || return false
        parse(Float64, s[end-1-k]) == pressure || return false
        return true
    end
end

const DATAenergy = joinpath(@__DIR__, "DATAenergy")
function makeDATAenergy(previous, target="total potential")
    @assert !isfile(DATAenergy)
    prepends = Dict{Tuple{Int,String},SubString{String}}()
    if !isnothing(previous)
        old = joinpath(previous, "DATAenergy")
        @assert isfile(old)
        for l in eachline(old)
            splits = split(l, ','; limit=3)
            prepends[(parse(Int,splits[1]),splits[2])] = splits[3]
        end
    end

    ptperf = Tuple{String,Vector{Float64}}[]
    ptperf_dict = Dict{String,Int}()
    tab = target == "total potential" ? "" : "	"
    for temp in readdir(joinpath(@__DIR__, "Output"), join=false)
        files = readdir(joinpath(@__DIR__, "Output", temp), join=true)
        TEMP = seekgrep!(eachline(first(files)), Int, r"^External temperature:\s*([^\s]+)")
        for system in files
            splits = split(basename(system), '_')
            bbase = join(@view(splits[2:end-1]), '_')
            pressure = parse(Float64, first(splitext(splits[end])))
            energies = allocurrences(Float64, system, Regex("^$(tab)Current $target energy:\\s*([^\\s]+)"))
            num = length(energies)
            movietemp = joinpath(@__DIR__, "Movies", temp)
            if isdir(movietemp)
                moviefile = getfirst(_moviefinder(bbase, pressure, "frameworks"), readdir(movietemp, join=false))
                if isnothing(moviefile)
                    moviefile = getfirst(_moviefinder(bbase, pressure, "allcomponents"), readdir(movietemp, join=false))
                end
                if isnothing(moviefile)
                    printstyled("MOVIE for ", bbase, '_', pressure, " NOT FOUND\n"; color=:blue)
                else
                    num = count(l -> @view(l[1:5]) == "MODEL", eachline(joinpath(movietemp, moviefile)))
                end
            end
            ptocc = allocurrences(Float64, system, r"^System \[[0-9]+\]\<\-\>\[[0-9]+\].*\(([^\s]+)\s*\[\%\]\)$")
            ptperfkey = join(@view(split(bbase, '_')[1:3]), '_')
            idptperf = get!(ptperf_dict, ptperfkey, length(ptperf)+1)
            !isempty(ptocc) && if idptperf > length(ptperf)
                push!(ptperf, (ptperfkey, [last(ptocc)]))
            else
                push!(last(ptperf[idptperf]), last(ptocc))
            end
            # !isempty(ptocc) && push!(thisptperf, last(ptocc))
            open(DATAenergy, "a") do f
                key = (TEMP, basename(system))
                print(f, TEMP, ',', key[2], ',')
                pre = get(prepends, key, nothing)
                if !isnothing(pre)
                    print(f, pre, ',')
                end
                join(f, @view(energies[end-num+1:end]), ',')
                println(f)
            end
        end
        # push!(ptperf, thisptperf)
    end
    !isempty(ptperf) && open(joinpath(@__DIR__, "DATAptperf"), "w") do f
        for (bbase, ptperfs) in ptperf
            print(f, bbase, ": ")
            join(f, ptperfs, ',')
            println(f)
        end
        # for i in 1:length(first(ptperf))
        #     join(f, (ptperf[j][i] for j in eachindex(ptperf)), ',')
        #     println(f)
        # end
    end
end

if isempty(ARGS)
    makeDATAenergy(nothing)
elseif occursin('=', ARGS[1])
    _s = split(ARGS[1], '=')
    @assert _s[1] == "target"
    @assert length(_s) == 2
    makeDATAenergy(nothing, _s[2])
else
    makeDATAenergy(ARGS[1])
end

const gd = groupby(CSV.read(DATAenergy, DataFrame; header=false), 1)
const nums = combine(gd, 1=>first; renamecols=false)[!, 1]
for (df, T) in zip(gd, nums)
    titles = [first(splitext(last(split(x, '_'; limit=2)))) for x in @view(df[:,2])]
    perm = rename!(permutedims(@view df[:, 3:end]), titles)
    CSV.write(joinpath(@__DIR__, "DATAenergy_$T.csv"), perm)
end
rm(DATAenergy)
