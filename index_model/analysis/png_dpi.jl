#####################################################################
#####################################################################
### Stamp the physical resolution into a PNG.
###
### GR rasterises at whatever `dpi` Plots asks for - a 9 x 2.8 inch
### figure at dpi = 600 really is written as 5400 x 1680 pixels - but
### it records no PHYSICAL size in the file.  A PNG without a pHYs
### chunk is read as 72 dpi by `sips`, by Preview, by Word, and by the
### submission systems that check "is this at least 300 dpi", so a
### perfectly sharp figure gets reported as 72 dpi and bounced.
###
### This writes the pHYs chunk so the file says what it is.  It does
### not resample anything: the pixels are untouched, only the metadata
### that says how large they are meant to be.
###
### Standalone and dependency free - a PNG chunk is a length, a type, a
### payload and a CRC.  The one wrinkle is that PNG uses CRC-32 (IEEE),
### while Julia's CRC32c stdlib provides the Castagnoli variant, which
### is a different polynomial and would produce a file every reader
### rejects.  So the table is built here.
### Jack H. Buckner, Oregon State University, 2026
### Generated with Claude Code.
#####################################################################
#####################################################################

const PNG_SIGNATURE = UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]

# CRC-32 (IEEE 802.3), reflected, polynomial 0xEDB88320.
const CRC32_TABLE = let table = zeros(UInt32, 256)
    for n in 0:255
        c = UInt32(n)
        for _ in 1:8
            c = (c & 0x00000001) != 0 ? (0xedb88320 ⊻ (c >> 1)) : (c >> 1)
        end
        table[n + 1] = c
    end
    table
end

function crc32(data::AbstractVector{UInt8})
    crc = 0xffffffff
    for b in data
        crc = (crc >> 8) ⊻ CRC32_TABLE[((crc ⊻ b) & 0xff) + 1]
    end
    return crc ⊻ 0xffffffff
end

be32(x::Integer) = UInt8[(UInt32(x) >> 24) & 0xff, (UInt32(x) >> 16) & 0xff,
                         (UInt32(x) >> 8) & 0xff, UInt32(x) & 0xff]

# One PNG chunk: length, type, payload, CRC over type+payload.
function png_chunk(kind::String, payload::Vector{UInt8})
    body = vcat(Vector{UInt8}(kind), payload)
    return vcat(be32(length(payload)), body, be32(crc32(body)))
end


"""
    set_png_dpi!(path, dpi)

Rewrite `path` with its pHYs chunk set to `dpi`, so the file reports the
resolution it was actually rendered at.  Idempotent: an existing pHYs chunk is
replaced rather than duplicated.

The pixels are not touched - this only records how large they are meant to be.
Returns the path.
"""
function set_png_dpi!(path::AbstractString, dpi::Real)
    dpi > 0 || throw(ArgumentError("dpi must be positive, got $dpi"))
    data = read(path)
    length(data) > 8 && data[1:8] == PNG_SIGNATURE ||
        throw(ArgumentError("$path is not a PNG"))

    # pixels per metre, which is the only unit PNG defines (unit = 1)
    ppm = round(UInt32, dpi / 0.0254)
    phys = png_chunk("pHYs", vcat(be32(ppm), be32(ppm), UInt8[0x01]))

    out = IOBuffer()
    write(out, PNG_SIGNATURE)
    inserted = false
    i = 9
    while i <= length(data)
        len = (UInt32(data[i]) << 24) | (UInt32(data[i+1]) << 16) |
              (UInt32(data[i+2]) << 8)  |  UInt32(data[i+3])
        kind = String(data[i+4:i+7])
        chunk = data[i:(i + 11 + len)]

        if kind == "pHYs"
            i += 12 + len          # drop any existing one; ours replaces it
            continue
        end
        # pHYs must precede IDAT, so insert at the first one we reach
        if !inserted && kind in ("IDAT", "IEND")
            write(out, phys)
            inserted = true
        end
        write(out, chunk)
        i += 12 + len
    end
    inserted || write(out, phys)   # degenerate file with neither IDAT nor IEND

    write(path, take!(out))
    return path
end


"""
    savefig_dpi(plt, path; dpi)

`Plots.savefig` followed by `set_png_dpi!` for PNG targets, so every figure this
project writes reports its real resolution.  Other formats (PDF, SVG) are vector
and pass straight through.
"""
function savefig_dpi(plt, path::AbstractString; dpi::Real)
    Plots.savefig(plt, path)
    endswith(lowercase(path), ".png") && set_png_dpi!(path, dpi)
    return path
end
