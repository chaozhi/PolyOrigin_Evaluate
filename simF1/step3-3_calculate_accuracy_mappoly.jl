using PolyOrigin
cd(@__DIR__)
pwd()
Pkg.status("PolyOrigin")

datadir=joinpath(pwd(),"F1data")
workdir=joinpath(pwd(),"res_mappoly")
designid="F1"
sumfile=string("accuracy_F1_mappoly.txt")
resio=open(joinpath(workdir,sumfile),"w+")
id="ndoseerr,nphaseerr,nalleleerr,nparentgeno,assignerr,callerr,delfraction"
msg = string("dataid,DR,noff,rep,",id)
write(resio, msg,"\n")
# [200,100,50,30,20,15,10]
for dr = [0,0.5]
    noffls=[200,100,50,30,20,15,10]
    for noff=noffls
        for rep=1:50
            dataid=string(designid,"_noff", noff,"_DR",dr==0 ? "0" : dr,"_rep",rep)
            truefile = joinpath(datadir,string(dataid,"_polyorigin_truevalue_ancestral.csv"))
            pedfile = joinpath(datadir,string(dataid,"_polyorigin_pedigree.csv"))
            probfile = joinpath(workdir,string(dataid,"_genoprob_mappoly.csv"))
            polyancestry = readPolyAncestry(probfile,pedfile,
                workdir=workdir,verbose=false)
            acc =calAccuracy!(truefile,polyancestry,workdir=workdir)
            acc2=join(collect(values(acc)),",")
            println("dataid=",dataid,",acc=", acc2)
            ls =[dataid,dr,noff,rep,acc2]
            write(resio, join(ls,","),"\n")
            flush(resio)
        end
    end
end
close(resio)
