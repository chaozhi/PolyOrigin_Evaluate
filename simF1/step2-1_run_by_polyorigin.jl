using PolyOrigin
cd(@__DIR__)
pwd()

for byparent=[false,true]
    datadir=joinpath(pwd(),"F1data")
    workdir=joinpath(pwd(),"res_polyorigin")
    designid="F1"
    noffls=[10,15,20,30,50,100,200]
    sumfile=string("ressum_","F1_polyorigin_byparent",Int(byparent),".txt")
    logfile=string("F1_polyorigin_byparent",Int(byparent),".log")
    resio=open(joinpath(workdir,sumfile),"w+")
    id="ndoseerr,nphaseerr,nalleleerr,nparentgeno,assignerr,callerr,delfraction"
    msg = string("outid,tuse,designid,byparent,DR,noff,rep,",id)
    write(resio, msg,"\n")
    logio=open(joinpath(workdir,logfile),"w+")
    for dr = [0,0.5]
        for noff=noffls
            for rep=1:50
                dataid=string(designid,"_noff", noff,"_DR",dr==0 ? "0" : dr,"_rep",rep)
                outid = string(dataid,"_byparent",Int(byparent),"_output")
                genofile=joinpath(datadir,string(dataid,"_polyorigin_geno_snparray.csv"))
                pedfile =joinpath(datadir,string(dataid,"_polyorigin_pedigree.csv"))
                truefile = joinpath(datadir,string(dataid,"_polyorigin_truevalue_ancestral.csv"))
                tuse= @elapsed polyancestry=polyOrigin(genofile,pedfile,
                    chrpairing_phase=22,
                    chrpairing=44,
                    byparent=byparent,
                    workdir=workdir,
                    outstem=outid,
                    logfile=logio,
                    refhapfile=truefile
                )
                acc =calAccuracy!(truefile,polyancestry,workdir=workdir)
                acc2=join(collect(values(acc)),",")
                println("outid=",outid,", tuse=",round(tuse,digits=2), ",acc=", acc2)
                ls =[outid,tuse,designid,byparent,dr,noff,rep,acc2]
                write(resio, join(ls,","),"\n")
                flush(resio)
            end
        end
    end
    close(resio)
    close(logio)
end
