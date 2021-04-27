using PolyOrigin
cd(@__DIR__)
pwd()

genofile = "TableS2_dose.csv"
pedfile =  "TableS3_ped.csv"
@time polyancestry=polyOrigin(genofile,pedfile;
    isphysmap=true,
    recomrate=1.25,
    chrpairing=44,
    refinemap=true,
    refineorder=false,
    outstem="potato_output",
    isplot=true
)
