function ssGibbs(matrices::HybridMatrices,
                 geno::Genotypes,fixed::FixedMatrix,
                 ped::PedModule.Pedigree,
                 input::QTL.InputParameters;outFreq=5000)

    y   = matrices.y.full
    X   = matrices.X.full
    W   = matrices.W.full
    Zn  = matrices.Z.n
    Ai_nn = matrices.Ai.nn

    current   = Current(input,geno,fixed,y)
    current.β = zeros(length(current.β)+1) #add one for μ_g
    current.ϵ = zeros(matrices.num.pedn)

    output                 = QTL.Output(input,geno,fixed)
    output.meanFixdEffects = zeros(length(current.β)+1) #add one for μ_g
    output.meanEpsi        = zeros(matrices.num.pedn)


    wGibbs = GibbsMats(W)
    xGibbs = GibbsMats(X)

    if input.estimateVariance==false
      #construct lhs for sampleEpsilon!
      λ_ϵ             = input.varResidual/input.varGenotypic
      Z_n             = matrices.Z._n
      lhs_ϵ           = Z_n'Z_n+Ai_nn*λ_ϵ
      lhsCol_ϵ        = [lhs_ϵ[:,i] for i=1:size(lhs_ϵ,1)]
      lhsDi_ϵ         = 1.0./diag(lhs_ϵ)
      sd_ϵ            = sqrt(lhsDi_ϵ*current.varResidual)
      #construct lhs for sample marker effects
      λ               = current.varResidual/current.varEffect
      lhsDi           = [1.0/(wGibbs.xpx[i]+λ)::Float64 for i=1:size(wGibbs.xpx,1)]
      sd              = sqrt(lhsDi*current.varResidual)
    end

    println("running ",input.method," with a MCMC of length ",input.chainLength)

    for iter = 1:input.chainLength

      current.iter += 1
      # sample fixed effects
      sample_fixed!(xGibbs,current,output)
      # sample marker effects
      sample_random_ycorr!(wGibbs,current,output,lhsDi,sd)
      # sample epsilon
      sampleEpsi!(matrices,current,output,lhsCol_ϵ,lhsDi_ϵ,sd_ϵ)

      if (iter%outFreq ==0)
          println("This is iteration ",iter)
      end
    end

    estimatedMarkerEffects = output.meanMarkerEffects

    mu_g = output.meanFixdEffects[end]
    EBV = matrices.J.full*mu_g+matrices.M.full*estimatedMarkerEffects
    EBV[1:matrices.num.pedn,:] += output.meanEpsi

    IDs=PedModule.getIDs(ped);
    EBV=DataFrame(ID=IDs,EBV=vec(EBV))

    return EBV
end

export ssGibbs



