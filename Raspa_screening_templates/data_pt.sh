#!/bin/env julia

using CSV
using DataFrames

for file in readdir(@__DIR__, join=true)
    if startswith(basename(file), "DATA")
        rm(file)
    end
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

const DATA = joinpath(@__DIR__, "DATA")
@assert !isfile(DATA)

for temp in readdir(joinpath(@__DIR__, "Output"), join=true)
    files = readdir(temp, join=true)
    TEMP = seekgrep!(eachline(first(files)), Int, r"^External temperature:\s*([^\s]+)")
    for system in files
        energies = allocurrences(Float64, system, r"^Current total potential energy:\s*([^\s]+)")
        open(DATA, "a") do f
            print(f, TEMP, ',', basename(system), ',')
            join(f, energies, ',')
            println(f)
        end
    end
end

const gd = groupby(DataFrame(CSV.File(DATA; header=false)), 1)
const nums = combine(gd, 1=>first; renamecols=false)[!, 1]
for (df, T) in zip(gd, nums)
    titles = [split(x, '_')[end-3] for x in @view df[:,2]]
    perm = rename!(permutedims(@view df[:, 3:end]), titles)
    CSV.write("DATA_$T.csv", perm)
end
