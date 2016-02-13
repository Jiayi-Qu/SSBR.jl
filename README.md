# SSBR

SSBR is a tool for single step Bayesian regression analyses.


####Quick-start

```Julia
using QTL
using SSBR

#data files form QTLDatasets package
pedfile    = QTLDatasets.dataset("test1","ped.txt")
genofile   = QTLDatasets.dataset("test1","genotype.txt")
phenofile  = QTLDatasets.dataset("test1","phenotype.txt")
Validation = QTLDatasets.dataset("test1","validation.txt")

#set up input parameters
input=InputParameters()
input.method       = "BayesC"
input.varGenotypic = 4.48
input.varResidual  = 6.72
input.probFixed    = 0.99
input.chainLength  = 50000

#run it
out=runSSBR(input,pedigree=pedfile,genotype=genofile,phenotype=phenofile);

#check accuracy
using DataFrames
df = readtable(Validation, eltypes =[UTF8String, Float64], separator = ' ',header=false,names=[:ID,:EBV]);
comp=join(out,df,on=:ID);
cor(comp[:EBV],comp[:EBV_1])
```

####More

* **homepage**: [QTL.rocks](http://QTL.rocks)
* **Installation**: at the Julia REPL, `Pkg.clone("https://github.com/QTL-rocks/SSBR.jl.git")`
* **Documentation**: [available here](https://github.com/QTL-rocks/SSBR.jl/wiki)
* **Authors**: [Hao Cheng](http://reworkhow.github.io),[Rohan Fernando](http://www.ans.iastate.edu/faculty/index.php?id=rohan)
