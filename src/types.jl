type Numbers
  ped::Int64      #individulals in pedigree
  pedn::Int64 #g1        #non-genotyped individuals in pedigree
  pedg::Int64 #g2        #genotyped individuals in pedigree
  y::Int64        #individuals with phenotypes
  yn::Int64       #non-genotyped individuals with phenotypes
  yg::Int64       #genotyped individuals with phenotypes
  markers::Int64  #number of markers
end

type ZMats
  Z::SparseMatrixCSC{Float64,Int64}    # Z= | Zn 0 |=|Z_n Z_g|
  Zn::SparseMatrixCSC{Float64,Int64}   #    | 0  Zg|
  Zg::SparseMatrixCSC{Float64,Int64}   #
  Z_n::SparseMatrixCSC{Float64,Int64}  # Z_n=|Zn|
  Z_g::SparseMatrixCSC{Float64,Int64}  #     |0 |
end

type AMats
  Ai::SparseMatrixCSC{Float64,Int64}
  Ai_nn::SparseMatrixCSC{Float64,Int64} #Ai11
  Ai_ng::SparseMatrixCSC{Float64,Int64} #Ai12
end

type YVecs
  y::Array{Float64,1}
  yn::Array{Float64,1} #y1
  yg::Array{Float64,1} #y2  #order of ids is same to order of y
  ids::Array{ASCIIString,1} #order of ids is nongeno then geno
end

type MMats
  M::Array{Float64,2}
  Mn::Array{Float64,2} #M1
  Mg::Array{Float64,2} #M2
end

type JVecs
  J::Array{Float64,2}
  Jn::Array{Float64,2}
  Jg::Array{Float64,2}
end

type XMats #sparse may be better
  X::Array{Float64,2}
  Xn::Array{Float64,2}
  Xg::Array{Float64,2}
end

type WMats
  W::Array{Float64,2}
  Wn::Array{Float64,2}
  Wg::Array{Float64,2}
end

type HybridMatrices
  zmats::ZMats
  amats::AMats
  yvecs::YVecs
  jvecs::JVecs
  xmats::XMats
  wmats::WMats
  mmats::MMats
  num::Numbers
end
