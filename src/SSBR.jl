module SSBR

using QTL
using PedModule
using DataFrames

include("SSBR_types.jl")
include("getMatrices.jl")
include("sampleEpsilon.jl")
include("ssBayesC.jl")
include("ssBayesC0.jl")
include("ssBayesC0_constantvariance.jl")#for test
#include("ssMME.jl")

export calc_Ai
export make_MMats
export make_yVecs
export make_JVecs
export make_ZMats
export make_XWMats
#export ssGibbs
#export ssMME
#export PBLUP


end
