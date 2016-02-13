function runSSBR(input;pedigree="pedfile",genotype="genofile",phenotype="phenofile")

      srand(input.seed)

      ped,geno,hmats =make_matrices_hybrid(pedigree,genotype,phenotype);

      hmats.M.g=zeros(1,1) #move inside later
      hmats.M.n=zeros(1,1)
      gc()

      fixed=QTL.FixedMatrix(ones(90,1),[0]); #modify later

      if input.method=="BayesC0" && input.estimateVariance==false
        out=ssBayesC0_constantvariance!(hmats,geno,fixed,ped,input,outFreq=5000)
      elseif input.method=="BayesC0"
        out=ssBayesC0(hmats,geno,fixed,ped,input,outFreq=5000)
      elseif input.method=="BayesB"
        out=ssBayesB(hmats,geno,fixed,ped,input,outFreq=5000)
      elseif input.method=="BayesC"
        out=ssBayesC(hmats,geno,fixed,ped,input,outFreq=5000)
      end
      return out
end


