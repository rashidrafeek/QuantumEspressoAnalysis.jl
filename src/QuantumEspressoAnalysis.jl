module QuantumEspressoAnalysis

using DataFrames, DelimitedFiles, Printf
using Unitful, UnitfulAtomic

"""
    extract_final_energy(outfile)

Extracts the final energy (in eV) from a Quantum Espresso output file.

# Arguments
- `outfile::String`: Path to the Quantum Espresso output file.

# Returns
- `Float64`: Final energy in eV.
"""
function extract_final_energy(outfile)
    endat = last(readlines(`grep ! $outfile`))
    energy = uconvert(u"eV", parse(Float64, split(endat)[5])u"Ry")
    return energy
end

"""
    compute_enthalpy_from_output(filename; energy_unit=u"Ry")

Computes the enthalpy (energy + pressure * volume) from a Quantum Espresso output file.

# Arguments
- `filename::String`: Path to the Quantum Espresso output file.
- `energy_unit::Unitful.FreeUnits`: Unit for energy (default: Ry).

# Returns
- `Float64`: Enthalpy in the specified energy unit.
"""
function compute_enthalpy_from_output(filename; energy_unit=u"Ry")
    dat = split.(readlines(pipeline(
        `grep -E 'unit-cell volume|!|P=' $filename`, `tail -n3`
    )))

    volume = parse(Float64, dat[1][4])u"bohr^3"
    energy = parse(Float64, dat[2][5])u"Ry"
    pressure = parse(Float64, dat[3][6])u"kbar"

    enthalpy = uconvert(energy_unit, energy + pressure * volume)
    return enthalpy
end

"""
    extract_final_pressure(filename; pressure_unit=u"GPa")

Extracts the final pressure from a Quantum Espresso output file.

# Arguments
- `filename::String`: Path to the Quantum Espresso output file.
- `pressure_unit::Unitful.FreeUnits`: Unit for pressure (default: GPa).

# Returns
- `Float64`: Final pressure in the specified unit.
"""
function extract_final_pressure(filename; pressure_unit=u"GPa")
    dat = split(last(readlines(`grep -E 'P=' $filename`)))
    pressure = parse(Float64, dat[6])u"kbar"
    return uconvert(pressure_unit, pressure)
end

"""
    extract_final_volume(outfile)

Extracts the final unit-cell volume from a Quantum Espresso output file.

# Arguments
- `outfile::String`: Path to the Quantum Espresso output file.

# Returns
- `Float64`: Final unit-cell volume in bohr³.
"""
function extract_final_volume(outfile)
    voldat = split(last(readlines(`grep 'unit-cell volume' $outfile`)))
    volume = parse(Float64, voldat[4])u"bohr^3"
    return volume
end

"""
    extract_dispersion_correction_energy(filename; energy_unit=u"eV")

Extracts the dispersion correction energy from a Quantum Espresso output file.

# Arguments
- `filename::String`: Path to the Quantum Espresso output file.
- `energy_unit::Unitful.FreeUnits`: Unit for energy (default: eV).

# Returns
- `Float64`: Dispersion correction energy in the specified unit.
"""
function extract_dispersion_correction_energy(filename; energy_unit=u"eV")
    dat = split(last(readlines(`grep 'Dispersion Correction     =' $filename`)))
    energy = parse(Float64, dat[4])u"Ry"
    return uconvert(energy_unit, energy)
end

"""
    compute_enthalpies_for_multiple_files(directory, regex, file_format; natoms=40)

Computes enthalpies for multiple Quantum Espresso output files matching a pattern.

# Arguments
- `directory::String`: Directory containing the output files.
- `regex::Regex`: Regex to identify relevant files.
- `file_format::String`: Format string to construct file paths (e.g., "file_%d.out").
- `natoms::Int`: Number of atoms in the system (default: 40).

# Returns
- `DataFrame`: DataFrame with columns `:pressure`, `:enthalpy`, and `:enthalpy_per_atom`.
"""
function compute_enthalpies_for_multiple_files(directory, regex, file_format; natoms=40)
    file_paths = joinpath.(directory, filter(contains(regex), readdir(directory)))
    pressures = sort(parse.(Int, only.(match.(regex, file_paths))))
    formatted_files = Printf.format.(Ref(Printf.Format("$directory/$file_format")), pressures)

    if !all(isfile.(formatted_files))
        missing_file = findfirst(==(false), isfile.(formatted_files))
        error("File not found: $(formatted_files[missing_file])")
    end

    enthalpies = compute_enthalpy_from_output.(formatted_files)

    return DataFrame(
        :pressure => pressures .* u"GPa",
        :enthalpy => enthalpies,
        :enthalpy_per_atom => enthalpies ./ natoms
    )
end

"""
    parse_pdos_data(filename, l)

Parses a projected density of states (PDOS) file for a given angular momentum.

# Arguments
- `filename::String`: Path to the PDOS file.
- `l::Int`: Angular momentum quantum number (0: s, 1: p, 2: d, 3: f).

# Returns
- `DataFrame`: DataFrame with energy and PDOS columns.
"""
function parse_pdos_data(filename, l)
    mult = 2 * l + 1
    header_regex = raw"^# *E \(eV\) +ldosup\(E\) +ldosdw\(E\)" *
                   repeat(raw" +pdosup\(E\) +pdosdw\(E\)", mult) * raw" *$"
    header = readline(filename)
    !contains(header, Regex(header_regex)) && error("Unexpected header: $header")

    data = readdlm(filename; comments=true)
    size(data, 2) != 3 + mult * 2 && error("Unexpected number of columns")

    lname = l == 0 ? "s" : l == 1 ? "p" : l == 2 ? "d" : l == 3 ? "f" : error("Invalid angular momentum: $l")
    pdos_names = lname .* "_" .* string.(1:mult)
    pdos_spin_names = vec(stack(pdos_names) do name name .* ["_up", "_down"] end)

    return DataFrame(
        :energy => data[:, 1]u"eV",
        (Symbol.(pdos_spin_names) .=> eachcol(data[:, 4:end]))...
    )
end

"""
    combine_pdos_files(filenames...)

Combines PDOS data from multiple files into a single DataFrame.

# Arguments
- `filenames::Vararg{String}`: List of PDOS file paths.

# Returns
- `DataFrame`: Combined DataFrame with energy and PDOS data.
"""
function combine_pdos_files(filenames...)
    l_vals = replace.(only.(match.(r"_wfc#[0-9]\((.+)\)", filenames)), "s" => 0, "p" => 1, "d" => 2, "f" => 3)
    data_frames = map(parse_pdos_data, filenames, l_vals)

    if length(data_frames) > 1
        !allequal(getproperty.(data_frames, :energy)) && error("Energy columns differ across files")
        return innerjoin(data_frames...; on=:energy)
    else
        return only(data_frames)
    end
end

"""
    extract_magnetization_from_output(filename)

Extracts total and absolute magnetization from a Quantum Espresso SCF output file.

# Arguments
- `filename::String`: Path to the SCF output file.

# Returns
- `Tuple{Float64, Float64}`: Total and absolute magnetization values.
"""
function extract_magnetization_from_output(filename)
    tot_mag_line = only(readlines(pipeline(`grep "total magnetization" $filename`, `tail -n1`)))
    abs_mag_line = only(readlines(pipeline(`grep "absolute magnetization" $filename`, `tail -n1`)))

    total_magnetization = parse(Float64, split(tot_mag_line)[4])
    absolute_magnetization = parse(Float64, split(abs_mag_line)[4])

    return total_magnetization, absolute_magnetization
end

end # module
