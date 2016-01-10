module SSBR

using QTL
using PedModule
using DataFrames

include("types.jl")
include("functions.jl")
#include("ssGibbs.jl")
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
