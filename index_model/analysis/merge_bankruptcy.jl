#####################################################################
#####################################################################
### Merge the per scenario bankruptcy tables into one tidy CSV.
###
### run_bankruptcy.jl is submitted as one HPC job per scenario, so
### each job writes its own part table into manuscript_index/parts/.
### This combines them, in sweep order, into
### manuscript_index/index_bankruptcy.csv.
###
### Run wherever the parts are - on the cluster after the jobs finish,
### or locally after `make sync_index_analysis`:
###
###     julia --project=. index_model/analysis/merge_bankruptcy.jl
###     julia --project=. index_model/analysis/merge_bankruptcy.jl \
###           --parts=DIR --out=PATH
###
### Deliberately standalone: no model, no ValueFunctionIterations, no
### Plots, so it starts in under a second and can be rerun freely while
### jobs are still landing.  It is also why the CSV reader below is
### written out rather than imported.
### Jack H. Buckner, Oregon State University, 2026
### Generated with Claude Code.
#####################################################################
#####################################################################

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_PARTS = joinpath(ROOT, "manuscript_index", "parts")
const DEFAULT_OUT   = joinpath(ROOT, "manuscript_index", "index_bankruptcy.csv")

# Sweep order for the merged table.  Arm A first, then arm B, each down the
# precision ladder, so the CSV reads the way the two arms are plotted.  A scenario
# whose arm is not one of these sorts last rather than being dropped.
const ARM_ORDER = Dict("iso_premium" => 1, "full_coverage" => 2)


#####################################################################
### A minimal CSV reader
#####################################################################
#
# Matches what setup.jl's write_csv produces: fields are quoted only when they
# contain a comma, a quote or a newline, and an embedded quote is doubled.
# Newlines inside quoted fields are not produced by write_csv and so are not
# handled here; a field containing one would be an error, not a silent misparse.

function split_csv_line(line::AbstractString)
    fields = String[]
    buf = IOBuffer()
    inquote = false
    i = firstindex(line)
    while i <= lastindex(line)
        c = line[i]
        if inquote
            if c == '"'
                # a doubled quote is a literal quote; a single one closes the field
                if i < lastindex(line) && line[nextind(line, i)] == '"'
                    write(buf, '"')
                    i = nextind(line, i)
                else
                    inquote = false
                end
            else
                write(buf, c)
            end
        elseif c == '"'
            inquote = true
        elseif c == ','
            push!(fields, String(take!(buf)))
        else
            write(buf, c)
        end
        i = nextind(line, i)
    end
    inquote && throw(ArgumentError("unterminated quoted field in: $line"))
    push!(fields, String(take!(buf)))
    return fields
end


"""
    read_part(path)

Read one part table.  Returns `(header, rows)` with `rows` a vector of the raw
data lines, unparsed: the merge only needs to reorder whole rows, and reparsing
the numbers would risk changing them.
"""
function read_part(path::AbstractString)
    lines = filter(!isempty, readlines(path))
    isempty(lines) && return (header = nothing, rows = String[])
    return (header = lines[1], rows = lines[2:end])
end


# Sort key for a data row, read off its own fields so the merge needs no
# knowledge of the sweep beyond the column names.
function row_key(line, cols)
    f = split_csv_line(line)
    get_col(name, default) = begin
        i = findfirst(==(name), cols)
        (i === nothing || i > length(f)) ? default : f[i]
    end
    arm = get_col("arm", "")
    pr  = something(tryparse(Float64, get_col("pr", "")), -Inf)
    # descending precision within an arm: the ladder runs 1.00 down to 0.33
    return (get(ARM_ORDER, arm, length(ARM_ORDER) + 1), -pr,
            get_col("cost_tag", ""), get_col("scenario", ""))
end


#####################################################################
### Main
#####################################################################

function parse_args(args)
    parts, out = DEFAULT_PARTS, DEFAULT_OUT
    for a in args
        if startswith(a, "--parts=");  parts = String(a[9:end])
        elseif startswith(a, "--out="); out = String(a[7:end])
        else
            throw(ArgumentError("unknown option $a; expected --parts=DIR or --out=PATH"))
        end
    end
    return (parts = parts, out = out)
end


function main(args)
    opts = parse_args(args)

    isdir(opts.parts) || begin
        println("no parts directory at ", opts.parts)
        println("run `make run_index_bankruptcy` on the cluster first, then ",
                "`make sync_index_analysis` locally")
        return 1
    end

    files = sort([joinpath(opts.parts, f) for f in readdir(opts.parts)
                  if endswith(f, ".csv")])
    isempty(files) && begin
        println("no .csv part tables in ", opts.parts)
        return 1
    end

    header = nothing
    rows   = String[]
    used   = String[]
    empty_files = String[]

    for f in files
        part = read_part(f)
        # A job that died before writing anything leaves an empty or header only
        # file.  Reported, not merged, and not fatal: the point of the fan out is
        # that one failure does not cost the rest.
        if part.header === nothing || isempty(part.rows)
            push!(empty_files, basename(f))
            continue
        end
        if header === nothing
            header = part.header
        elseif part.header != header
            throw(ErrorException(
                "column mismatch: $(basename(f)) has a different header from " *
                "$(basename(first(used))).  The parts were written by different " *
                "versions of run_bankruptcy.jl; re-run the sweep rather than " *
                "merging them."))
        end
        append!(rows, part.rows)
        push!(used, f)
    end

    if header === nothing
        println("every part table in ", opts.parts, " is empty; nothing to merge")
        return 1
    end

    cols = split_csv_line(header)
    sort!(rows, by = line -> row_key(line, cols))

    mkpath(dirname(abspath(opts.out)))
    open(opts.out, "w") do io
        println(io, header)
        foreach(r -> println(io, r), rows)
    end

    scen_i = findfirst(==("scenario"), cols)
    scenarios = scen_i === nothing ? String[] :
        unique([split_csv_line(r)[scen_i] for r in rows])

    println("merged ", length(used), " part tables -> ", opts.out)
    println("  ", length(rows), " rows, ", length(scenarios), " scenarios")
    isempty(empty_files) || println("  skipped (empty, the job probably failed): ",
                                    join(empty_files, ", "))
    return 0
end


if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
