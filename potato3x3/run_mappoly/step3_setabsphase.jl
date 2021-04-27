using PolyOrigin
cd(@__DIR__)
pwd()
Pkg.status("PolyOrigin")

@time for pop=1:3
    for chr=1:12
        println("pop=",pop,", chr=",chr)
        mappolyfile =string("potato_pop",pop,"_chr",chr,"_genoprob_mappoly_correctid.csv")
        # isfile(mappolyfile)
        pedfile =string("potato_pop",pop,"_ped.csv")
        pedfile = joinpath("..//step1_data",pedfile)
        # isfile(pedfile)
        refhapfile =string("refhap//potato_diallel_output_genoprob_pop",pop,"_chr",chr,".csv")
        # isfile(refhapfile)
        polyancestry = readPolyAncestry(mappolyfile,pedfile)
        PolyOrigin.setAbsPhase!(refhapfile,polyancestry)
        outfile=replace(mappolyfile,"_correctid.csv"=>"_correctid_absphase.csv")
        PolyOrigin.savegenoprob(outfile,polyancestry)
    end
end
