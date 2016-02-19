# SSBR

SSBR is a tool for single step Bayesian regression analyses.


####Quick-start

```Julia
using QTL
using QTLDatasets
using SSBR

#data files from QTLDatasets package
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

MCMCinfo(input)
#MCMC Information:
#seed                        314
#chainLength               50000
#method                   BayesC
#outFreq                    1000
#probFixed                 0.990
#varGenotypic              4.480
#varResidual               6.720
#estimateVariance           true
#estimatePi                false
#estimateScale             false
#dfEffectVar               4.000
#nuRes                     4.000
#nuGen                     4.000
#centering                 false


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
* **Installation**: at the Julia REPL, `Pkg.clone("https://github.com/QTL-rocks/SSBR.jl.git")`.  ( install [QTL.jl](https://github.com/QTL-rocks/QTL.jl) at first )
* **Documentation**: [available here](https://github.com/QTL-rocks/SSBR.jl/wiki)
* **Authors**: [Hao Cheng](http://reworkhow.github.io),[Rohan Fernando](http://www.ans.iastate.edu/faculty/index.php?id=rohan)
