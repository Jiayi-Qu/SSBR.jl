function ssGibbs(matrices::HybridMatrices,geno::Genotypes,inputParameters::QTL.InputParameters;outFreq=5000)

    y   = matrices.y.full
    X   = matrices.X.full
    W   = matrices.W.full
    Z11 = matrices.Z.n
    Ai11= matrices.Ai.nn

    sum2pq    = geno.sum2pq
    varEffects = inputParameters.varGenotypic/((1-inputParameters.probFixed)*sum2pq)


    mu   = mean(y)
    yCorr= y - mu
    β    = [mu, 0.0]
    α    = zeros(Float64,matrices.num.markers)
    ϵ    = zeros(Float64,matrices.num.pedn)



    meanBeta  = [0.0, 0.0]
    meanAlpha = zeros(Float64,matrices.num.mbarkers)
    meanEpsi  = zeros(Float64,matrices.num.pedn)

    wGibbs = GibbsMats(W)
    xGibbs = GibbsMats(X)

    zpz = diag(Z11'Z11)

    #construct lhs for sampleEpsilon!
    λ_epsilon       = vRes/vG
    Z_1             = all_Z.Z_1
    lhs_epsilon     = Z_1'Z_1+Ai11*λ_epsilon
    lhsCol_epsilon  = [lhs_epsilon[:,i] for i=1:size(lhs_epsilon,1)]
    lhsDi_epsilon   = 1.0./diag(lhs_epsilon)
    sd_epsilon      = sqrt(lhsDi_epsilon*vRes)

    if inputParameters.estimateVariance==false
      #construct lhs for sampleEpsilon!
      λ_epsilon       = vRes/vG
      Z_1             = all_Z.Z_1
      lhs_epsilon     = Z_1'Z_1+Ai11*λ_epsilon
      lhsCol_epsilon  = [lhs_epsilon[:,i] for i=1:size(lhs_epsilon,1)]
      lhsDi_epsilon   = 1.0./diag(lhs_epsilon)
      sd_epsilon      = sqrt(lhsDi_epsilon*vRes)

      λ               = vRes/vAlpha
      lhsDi           = [1.0/(wpw[i]+λ)::Float64 for i=1:size(wpw,1)]
      sd              = sqrt(lhsDi*vRes)
    end

    println("running ",input.method," with a MCMC of length ",input.chainlength)

    for iter = 1:nIter

      iIter = 1.0/iter
      # sample fixed effects
      sample_fixed!(xGibbs,current,out)
      # sample marker effects
      sample_random_ycorr!(wGibbs,current,out,lhsDi,sd)
      # sample epsilon
      sampleEpsi!(all_Z,lhsCol_epsilon,lhsDi_epsilon,sd_epsilon,yCorr,ϵ,meanEpsi,iIter)

      if (iter%outFreq ==0)
          println("This is iteration ",iter)
      end
    end

    mu_g = meanBeta[2]

    alpha_hat = meanAlpha
    epsi_hat  = meanEpsi
    aHat = all_J.J*mu_g+all_M.M*alpha_hat
    aHat[1:all_num.num_g1,:] += epsi_hat
    return aHat,meanAlpha,meanBeta,meanEpsi
end

function sampleEpsi!(mats::HybridMatrices,current::QTL.Current,out::QTL.Output)
    yCorr         = current.yCorr
    varRes        = current.varResidual
    λ             = current.varResidual/current.varEffects
    iIter         = 1/current.iter
    Z_n           = mats.Z._n
    Ai_nn         = mats.Ai.nn
    meanEpsi      = out.meanEpsi
    ϵ             = current.ϵ

    #add back {Zn'Zn}_{i,i} *ϵ
    yCorr[:] = yCorr[:] + Z_n*ϵ
    rhs = Z_n'yCorr

    lhs = Z_n'Z_n+Ai_nn*λ

    sample_random_rhs!(lhs,rhs,current,out) #use this general function for sample epsilon(Gianola Book)

    yCorr[:] = yCorr[:] - Z_n*ϵ
end

function sampleEpsi!(mats::HybridMatrices,current::QTL.Current,out::QTL.Output)
    yCorr         = current.yCorr
    varRes        = current.varResidual
    λ             = current.varResidual/current.varEffects
    iIter         = 1/current.iter
    Z_n           = mats.Z._n
    Ai_nn         = mats.Ai.nn
    meanEpsi      = out.meanEpsi
    ϵ             = current.ϵ

    #add back {Zn'Zn}_{i,i} *ϵ
    yCorr[:] = yCorr[:] + Z_n*ϵ
    rhs = Z_n'yCorr

    lhs = Z_n'Z_n+Ai_nn*λ

    sample_random_rhs!(lhs,rhs,current,out) #use this general function for sample epsilon(Gianola Book)

    yCorr[:] = yCorr[:] - Z_n*ϵ
end

function sampleEpsi!(mats::HybridMatrices,current::QTL.Current,out::QTL.Output,lhsDi,sd)
    Z_n  = mats.Z._n
    yCorr= current.yCorr

    yCorr[:] = yCorr[:] + Z_n*ϵ
    rhs = Z_n'yCorr

    sample_random_rhs!(lhsCol,rhs,lhsDi,sd,ϵ,meanEpsi,iIter)

    yCorr[:] = yCorr[:] - Z_1*ϵ
end

