function ssBayesC0(matrices::HybridMatrices,
                 geno::Genotypes,fixed::FixedMatrix,
                 ped::PedModule.Pedigree,
                 input::QTL.InputParameters;outFreq=5000)

    y   = matrices.y.full
    X   = matrices.X.full
    W   = matrices.W.full
    Zn  = matrices.Z.n
    Ai_nn = matrices.Ai.nn

    current   = Current(input,geno,fixed,y)
    current.fixed_effects = zeros(length(current.fixed_effects)+1) #add one for μ_g
    current.imputation_residual = zeros(matrices.num.pedn)

    output                 = QTL.Output(input,geno,fixed)
    output.mean_fixed_effects = zeros(length(output.mean_fixed_effects)+1) #add one for μ_g
    output.mean_imputation_residual        = zeros(matrices.num.pedn)

    wGibbs = GibbsMats(W)
    xGibbs = GibbsMats(X)

    current.varEffect=input.varGenotypic/geno.sum2pq

    println("running ",input.method," with a MCMC of length ",input.chainLength)

    for iter = 1:input.chainLength

      current.iter += 1
      # sample fixed effects
      sample_fixed!(xGibbs,current,output)
      # sample marker effects
      sample_random_ycorr!(wGibbs,current,output)
      # sample epsilon
      sampleEpsi!(matrices,current,output)
      # sample residual vairance
      current.varResidual=sample_variance(current.yCorr, geno.nObs, input.nuRes, current.scaleRes)
      # sample marker vairance
      current.varEffect = sample_variance(current.α, geno.nMarkers, input.dfEffectVar, current.scaleVar)
      # sample marker vairance
      current.varGenotypic = sample_epsilon_variance(current.imputation_residual,
                                                     Ai_nn,
                                                     matrices.num.pedn,
                                                     input.dfEffectVar,
                                                     current.scaleVar)

      if (iter%outFreq ==0)
          println("This is iteration ",iter," with residual variance ",
                  current.varResidual," and marker variance ", current.varEffect)
      end
    end

    estimatedMarkerEffects = output.meanMarkerEffects

    mu_g = output.mean_fixed_effects[end]
    EBV = matrices.J.full*mu_g+matrices.M.full*estimatedMarkerEffects
    EBV[1:matrices.num.pedn,:] += output.mean_imputation_residual

    IDs=PedModule.getIDs(ped);
    EBV=DataFrame(ID=IDs,EBV=vec(EBV))

    return EBV
end

export ssGibbs



