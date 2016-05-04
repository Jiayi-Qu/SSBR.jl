function runSSBR(input;pedigree="pedfile",genotype="genofile",phenotype="phenofile",fixedfile="fixedeffectfile")

      srand(input.seed)

      ped,geno,hmats =make_matrices_hybrid(pedigree,genotype,phenotype,center=input.centering);

      hmats.M.g=zeros(1,1) #move inside later
      hmats.M.n=zeros(1,1)
      gc()

      fixed=QTL.make_fixed(fixedfile,ID_order=hmats.y.ids)
    
       Xn  = hcat(ones(hmats.num.yn),fixed.C[1:hmats.num.yn,:], hmats.Z.n*hmats.J.n)#intercept, fixed efeectsJ
       Xg  = hcat(ones(hmats.num.yg),fixed.C[(hmats.num.yn+1):end,:], hmats.Z.g*hmats.J.g)
       X   =[Xn;
             Xg]
    hmats.X = XMats(X,Xn,Xg) #Done, work

    
      if input.method=="BayesC0" && input.estimateVariance==false
        out=ssBayesC0_constantvariance(hmats,geno,fixed,ped,input,outFreq=input.outFreq)
      elseif input.method=="BayesC0"
        out=ssBayesC0(hmats,geno,fixed,ped,input,outFreq=input.outFreq)
      elseif input.method=="BayesB"
        out=ssBayesB(hmats,geno,fixed,ped,input,outFreq=input.outFreq)
      elseif input.method=="BayesC" && input.estimateVariance==false
        out=ssBayesC_constantvariance(hmats,geno,fixed,ped,input,outFreq=input.outFreq)
      elseif input.method=="BayesC"
        out=ssBayesC(hmats,geno,fixed,ped,input,outFreq=input.outFreq)
      end
      return out
end


