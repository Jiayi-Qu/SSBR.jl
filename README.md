# SSBR

SSBR is a tool for single step Bayesian regression analyses.


####Quick-start

```Julia
using DataFrames
using QTL
using SSBR

PhenoFile  ="training.txt"
GenoFile   ="snp.txt"
PedFile    ="ped.txt"
Validation ="validation.txt"
fixed=FixedMatrix(ones(90,1),[0]) #change to read txt file later

input=InputParameters()
input.varGenotypic = 4.48
input.varResidual  = 6.72
input.probFixed    = 0.0
input.chainLength  =50000

ped,geno,hmats =SSBR.make_matrices_hybrid(PedFile,GenoFile,PhenoFile)

hmats.M.g=zeros(1,1) #move inside later
hmats.M.n=zeros(1,1)
gc()

###BayesC0
out =SSBR.ssBayesC0(hmats,geno,fixed,ped,input,outFreq=5000);
#out =SSBR.ssBayesC0_constantvariance(hmats,geno,fixed,ped,input,outFreq=100)

###check accuracy
df = readtable(Validation, eltypes =[UTF8String, Float64], separator = ' ',header=false,names=[:ID,:EBV]);
comp=join(out,df,on=:ID);
cor(comp[:EBV],comp[:EBV_1])
```

####More

* **homepage**: [QTL.rocks](http://QTL.rocks)
* **Installation**: at the Julia REPL, `Pkg.clone("https://github.com/QTL-rocks/SSBR.jl.git")`
* **Documentation**: [available here](https://github.com/QTL-rocks/SSBR.jl/wiki)
* **Authors**: [Hao Cheng](http://reworkhow.github.io),[Rohan Fernando](http://www.ans.iastate.edu/faculty/index.php?id=rohan)
