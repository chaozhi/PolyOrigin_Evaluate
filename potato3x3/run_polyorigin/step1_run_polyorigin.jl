using PolyOrigin
cd(@__DIR__)
pwd()

dataid="potato_diallel"
datadir="..//step1_data"
pedfile = joinpath(datadir,string(dataid,"_ped.csv"))
genofile = joinpath(datadir,string(dataid,"_dose.csv"))
outstem=string(dataid,"_output")
@time polyancestry=polyOrigin(genofile,pedfile,
    isphysmap=true,
    recomrate=1.25,
    chrpairing_phase=22,
    chrpairing=44,
    refinemap=true,
    refineorder=false,
    outstem=outstem,
    isplot=true
)

22932/3600.
