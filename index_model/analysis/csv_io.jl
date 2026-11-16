#####################################################################
#####################################################################
### Reading and writing the tidy tables the analysis produces.
###
### Written by hand rather than with CSV.jl: neither CSV nor DataFrames
### is in Project.toml (analysis/table_4.jl uses both and so does not
### run under --project=.), and these tables are flat rows of numbers
### and short strings.  Keeping the index model free of the extra
### dependencies means the same environment that solves also analyses.
###
### Standalone on purpose - no model, no ValueFunctionIterations, no
### Plots - so merge_bankruptcy.jl and run_composite_figure.jl can read
### tables back in under a second, without paying to load the solver
### they never call.  analysis/setup.jl includes it too, so every
### analysis script gets write_csv transparently.
### Jack H. Buckner, Oregon State University, 2026
### Generated with Claude Code.
#####################################################################
#####################################################################


#####################################################################
### Writing
#####################################################################

csv_field(x::AbstractString) =
    any(c -> c in (',', '"', '\n', '\r'), x) ? '"' * replace(x, '"' => "\"\"") * '"' : x
csv_field(x::Symbol)  = csv_field(string(x))
csv_field(x::Bool)    = string(x)
csv_field(x::Real)    = isnan(x) ? "NA" : string(x)   # NA reads as missing in R
csv_field(x::Nothing) = "NA"
csv_field(x)          = csv_field(string(x))

"""
    write_csv(path, rows; columns = keys(first(rows)))

Write a vector of NamedTuples as a tidy CSV, creating the directory if needed.
Every row must carry `columns`; a row missing one is an error rather than a
silently blank cell.
"""
function write_csv(path::AbstractString, rows::AbstractVector; columns = nothing)
    isempty(rows) && throw(ArgumentError("nothing to write to $path"))
    cols = something(columns, collect(keys(first(rows))))
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, join(csv_field.(cols), ","))
        for r in rows
            missing_cols = [c for c in cols if !haskey(r, c)]
            isempty(missing_cols) || throw(ArgumentError(
                "row is missing " * join(string.(missing_cols), ", ")))
            println(io, join([csv_field(r[c]) for c in cols], ","))
        end
    end
    return path
end


#####################################################################
### Reading
#####################################################################

"""
    split_csv_line(line) -> Vector{String}

Split one CSV record.  Matches what `write_csv` produces: fields are quoted only
when they contain a comma, a quote or a newline, and an embedded quote is
doubled.  Newlines inside quoted fields are not produced by `write_csv` and so
are not handled; a field containing one is an error, not a silent misparse.
"""
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


# Recover the value a field was written from.  `NA` is what write_csv emits for
# NaN, so it has to come back as NaN rather than as the string "NA" - otherwise a
# degenerate lambda_ratio silently turns a numeric column into a string one and
# every downstream comparison against it fails in a confusing way.
function parse_csv_value(s::AbstractString)
    s == "NA"    && return NaN
    s == "true"  && return true
    s == "false" && return false
    isempty(s)   && return ""
    v = tryparse(Int, s)
    v === nothing || return v
    f = tryparse(Float64, s)
    f === nothing || return f
    return String(s)
end


"""
    read_csv(path) -> Vector{NamedTuple}

Read a table written by `write_csv` back into the same shape it went out as:
one NamedTuple per row, with numeric fields parsed and `NA` restored to `NaN`.

Column names become symbols, so a row is indexed as `r.P10` or `r[:P10]` exactly
like the rows the analysis scripts build.  A ragged row - more or fewer fields
than the header - is an error rather than a silently shifted record.
"""
function read_csv(path::AbstractString)
    isfile(path) || throw(ArgumentError("no table at $path"))
    lines = filter(!isempty, readlines(path))
    isempty(lines) && throw(ArgumentError("$path is empty"))

    cols = Tuple(Symbol.(split_csv_line(lines[1])))
    rows = NamedTuple[]
    for (i, line) in enumerate(lines[2:end])
        f = split_csv_line(line)
        length(f) == length(cols) || throw(ArgumentError(
            "$path line $(i + 1) has $(length(f)) fields, expected $(length(cols))"))
        push!(rows, NamedTuple{cols}(Tuple(parse_csv_value.(f))))
    end
    return rows
end
