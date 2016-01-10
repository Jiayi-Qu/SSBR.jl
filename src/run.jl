

num=Numbers(0,0,0,0,0,0,0)

geno=genotypes("genotype.txt",id4row=true)


ped,A_Mats,numSSBayes = calc_Ai("example/ped.txt","example/genotype.ID")
df    = read_genotypes("example/genotype.txt",numSSBayes)
M_Mats = make_MMats(df,A_Mats,ped)
y_Vecs = make_yVecs("example/phenotype.txt",ped,numSSBayes)
J_Vecs = make_JVecs(numSSBayes,A_Mats)
Z_Mats = make_ZMats(ped,y_Vecs,numSSBayes)
X_Mats, W_Mats = make_XWMats(J_Vecs,Z_Mats,M_Mats,numSSBayes)

#Gibbs sampler
nIter  = 50000
vRes   = 1.0
vG     = 1.0
aHat,alphaHat,betaHat,epsiHat
=ssGibbs(M_Mats,y_Vecs,J_Vecs,Z_Mats,X_Mats,W_Mats,A_Mats,numSSBayes,vRes,vG,nIter);
