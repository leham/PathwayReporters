# Testing the pathway reporter method
using Distributions, Plots
import Main.PathFuncs; const PF=PathFuncs
"""
Time an example to estimate run time before the main loop.
Typically the t_max should be at least 20, and closer to
200 or 400 is more reliable (depending on other parameters).
To get useful reporters the Samplemax should be in the
hundreds or thousands.
Arguments:
lambda, mu, KN,KM,KP,deltap,a1,b1,a2,b2,a3,b3,t_max,samplesize
Each of KN has Beta(a1,b1) noise, scaled to the input mean value of KN
Each of lambda,mu,KM,KP have Beta(a2,b2) noise, scaled to the input mean value
deltap has Beta(a3,b3) noise, similarly scaled.
"""
@time PF.reporters(1,0,20,10,1,0.1,3,6,5,6,8,6,10,10)

""" Simulate scaled Beta Noise on parameters in a
burtsy model with maturation and translation.
Note that mRNA decay (delta) is fixed at 1.
For constitutve expression set lambda=1 and mu=0.
"""


"""
List of mean lambda values to be used:
"""
LList =[1.0]#2.0,4.0,8.0]

"""
List of mean mu values to be used:
"""
MList = [2.0]#4.0,8.0,16.0,32.0]


"""
List of hyper parameter values for Beta noise on transcription rate KN:
"""
KNBetaList = [
    [8.,6.],#eta^2 is 0.05
    [5.,6.],#eta^2 is 0.1
    [3.,6.],#eta^2 is 0.2
    [2.,6.],#eta^2 is 0.33
    [2.,12.]]#eta^2 is 0.4
"""
Peruse the Beta distribution shapes and squared coefficent of variation (noise strength):
"""
x=0:0.01:1
plot(x,betapdf, label="Beta PDF")
for ab in KNBetaLists
    println(PF.coeff_var_sqrd(ab[1], ab[2]))
end




"""
List of approximate intended Nascent mRNA averages.  These are estimates only.
"""
NList = [5]
"""
List of Maturation rates.
#Note that Mature mRNA copy number average
#will be the product of this and the nascent average
"""
MatList = [10]
"""
List of approximate intended protein averages.  These are estimates only.
"""
PList = [1000]

"""
List of protein decay rates.
"""
DeltaPList=[0.2]

"""
###################
 Simulation Results
###################
"""
Reporters= PF.ReporterResults(LList, MList,
    KNBetaList, NList, MatList, PList, DeltaPList, 500, 500,1.0)
"""
output the results as a csv file
"""
writedlm("Noise_Trials.csv", Reporters, ',')
