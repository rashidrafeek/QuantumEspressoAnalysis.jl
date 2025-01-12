module QuantumEspressoAnalysis

using DataFrames, DelimitedFiles, Printf
using Unitful, UnitfulAtomic

"""
getfinalenergy(outfile)

Get the final energy from the QE output file, `outfile`.
"""
function getfinalenergy(outfile)
	endat = last(readlines(`grep ! $outfile`))
	energy = uconvert(u"eV", parse(Float64, split(endat)[5])u"Ry")
end

"""
    getenthalpy_QEout(fname)

Grep the final energy, pressure and volume from the given QE output file and 
return the enthalpy in Ry.
"""
function getenthalpy_QEout(fname; enunit=u"Ry")
    dat = split.(readlines(pipeline(
        `grep -E 'unit-cell volume|!|P=' $fname`, `tail -n3`
       )))

    volume = parse(Float64, dat[1][4])u"bohr^3"
    energy = parse(Float64, dat[2][5])u"Ry"
    pressure = parse(Float64, dat[3][6])u"kbar"

    enthalpy = uconvert(enunit, energy + pressure*volume)
    return enthalpy
end

function getfinalpressure(fname; prunit=u"GPa")
    dat = split(last(readlines(`grep -E 'P=' $fname`)))

    pressure = parse(Float64, dat[6])u"kbar"

    return uconvert(prunit, pressure)
end

"""
getfinalvolume(outfile)

Get the final volume from the QE output file, `outfile`.
"""
function getfinalvolume(outfile)
    voldat = split(last(readlines(`grep 'unit-cell volume' $outfile`)))
    volume = parse(Float64, voldat[4])u"bohr^3"

    return volume
end

function getfinaldispersionenergy(fname; enunit=u"eV")
    dat = split(last(readlines(`grep 'Dispersion Correction     =' $fname`)))

    en = parse(Float64, dat[4])u"Ry"

    return uconvert(enunit, en)
end

"""
    getallQEenthalpies(path, prregex, outformat; natom=40)

Get a `DataFrame` of enthalpies of a set of calculation at different pressures
given by the files which match `prregex` and QE output file format given by
`outformat`, a string which can be read by the `@sprintf` macro in `Printf`.
The output contains columns `:pressure`, `:enthalpy` and `:enthalpyperatom`
computed based on the given number of atoms, `natom`.
"""
function getallQEenthalpies(path, prregex, outformat; natom=40)
    prpaths = joinpath.(path, filter(contains(prregex), readdir(path)))
    pressures = sort(parse.(Int, only.(match.(prregex, prpaths))))
    fullformat = Printf.Format("$path/$outformat")
    outfiles = Printf.format.(Ref(fullformat), pressures)

    fileexists = isfile.(outfiles)

    if all(fileexists)
        enthalpies = getenthalpy_QEout.(outfiles)
    else
        notexists = findfirst(==(false), fileexists)
        error("File with path $(outfiles[notexists]) does not exist.")
    end

    enthalpydat = DataFrame(
        :pressure=>pressures.*u"GPa",
        :enthalpy=>enthalpies,
        :enthalpyperatom=>enthalpies./natom
    )

    return enthalpydat
end

"""
    get_pdoscolumns(filename, l)

Obtain a dataframe containing the energy and corresponing pdos
values from a PDOS output file with name, `filename` and containing data for
angular momentum, `l`.
"""
function get_pdoscolumns(filename, l)
    # Ensure header is in the expected format
    mult = 2*l + 1
    headerformat = raw"^# *E \(eV\) +ldosup\(E\) +ldosdw\(E\)"
    for m in 1:mult
        headerformat *= raw" +pdosup\(E\) +pdosdw\(E\)"
    end
    headerformat *= raw" *$"
    header = readline(filename)
    !contains(header, Regex(headerformat)) && error("Unexpected header: $header")

    dat = readdlm(filename; comments=true)
    size(dat, 2) != 3 + mult*2 && error("Unexpected number of columns")

    lname = l == 0 ? "s" :
            l == 1 ? "p" :
            l == 2 ? "d" :
            l == 3 ? "d" :
            l == 4 ? "f" : error("Unknown value for l: $l")
    pdos_colnames = lname .* "_" .* string.(1:mult)
    pdos_spin_colnames = vec(stack(pdos_colnames) do lname
        lname .* ["_up", "_down"]
    end)
     
    return DataFrame(
        :energy=>dat[:,1]u"eV",
        (Symbol.(pdos_spin_colnames) .=> eachcol(dat[:,4:end]))...,
    )
end

"""
    get_pdos_atomdataframe(filenames...)

Get a dataframe containing energy and pdos values from several files.
"""
function get_pdos_atomdataframe(filenames...)
    # Detect the angular momentums from filenames
    l_names = only.(match.(r"_wfc#[0-9]\((.+)\)", filenames))
    lvals = replace(l_names, "s"=>0, "p"=>1, "d"=>2, "f"=>3)
    
    alldata = map(filenames, lvals) do filename, l
        get_pdoscolumns(filename, l)
    end
    
    if length(alldata) > 1
        allenergies = getproperty.(alldata, :energy)
        !allequal(allenergies) && error("Energy column differ between files: $filenames")
    
        return innerjoin(alldata...; on=:energy)
    else
        return only(alldata)
    end
end

"""
    get_scfout_magnetization(filename)

Get a tuple contaning total and absolute magnetization values respectively from
the given QE SCF calculation output file.
"""
function get_scfout_magnetization(filename)
    dat_tot = only(readlines(pipeline(`grep "total magnetization" $filename`, `tail -n1`)))
    dat_abs = only(readlines(pipeline(`grep "absolute magnetization" $filename`, `tail -n1`)))

    tot_mag = parse(Float64, split(dat_tot)[4])
    abs_mag = parse(Float64, split(dat_abs)[4])

    return tot_mag, abs_mag
end

end # module
