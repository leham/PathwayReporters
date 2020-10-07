# Various functions used in the module
using Distributions


"""
Function for plotting scaled Beta pdf
"""
betapdf(x)= pdf.(Beta(2,12),x)


""" Function to calculate the KN value
 for given model params and mean nascent mRNA"""

function KNval(lambda,mu,KM,Expected_nascent)
    return Expected_nascent * KM * ((lambda+mu)/lambda)
end

""" Function to calculate the KM value
for given model params and mean mature mRNA"""

function KMval(lambda,mu,Expected_mature)
    return Expected_mature * ((lambda+mu)/lambda)
end

""" Function to calculate the KP value
for given model params and mean protein"""

function KPval(lambda,mu,KN,delta,deltap,Expected_protein)
    return (Expected_protein * deltap * (lambda+mu)) / (KN * lambda)
end

""" Function to calculate the squared coefficent of variation
for a Beta distribution Beta(alpha, beta
"""

function coeff_var_sqrd(alpha, beta)
    b = alpha * (alpha + beta + 1)
    t = beta
    return t/b
end

""" Function for producing values sampled
 from an Expotential distribution"""

rexp(x)=rand(Exponential(x))

""" Mean nascent mRNA """
function nascentRNAmean(lambda,mu,KM,KN)
    return (KN/KM) * (lambda/(lambda+mu))
end

""" Mean protien"""
function Proteinmean(lambda,mu,K,Pr,deltap)
        return (K/delta) * (lambda/(lambda+mu)) * (Pr/deltap)
end

"""
Function to sample from a simulated telegraph system with
maturation and protein translation.
On rate: lambda
Off rate: mu
Transcription rate: Kn
Maturation rate: KM
Protein translation rate: KP
Protein degradation: deltap
mRNA decay is fixed at 1
NOTE: Maturation assumed high enough to ignore degradation in nascent phase.
Sample is taken once t_max time units is reached.
"""

function copynumbers(lambda,mu,KN,KM,KP,deltap,t_max)
    s = 1 # Set gene state on
    n = 0 # Initial mRNA copy number
    m = 0 # Initial mature copy number
    p = 0 # Initial protein copoy number
    t = 0 # Initial time
    delta = 1 # fixed mRNA decay rate

    while t <= t_max
        L = rexp(1/lambda) #turn on
        U = rexp(1/mu) # turn off
        mRNA = rexp(1/KN) # make mRNA
        MatmRNA = rexp(1/(n * KM)) #mature mRNA
        Md = rexp(1/(m*delta)) # decay mRNA
        makePr = rexp(1/(m * KP)) # translate a protein
        Pd = rexp(1/(p * deltap)) # decay protein
        if s == 0 && n == 0 && m == 0 && p == 0 #0000
            t += L
            s = 1
        elseif s == 0 && n == 0 && m == 0 && p >= 1 #0001
            mincase = min(Pd,L)
            t += mincase
            if mincase == L
                s = 1
            else
                p += -1
            end
            ##
        elseif s == 1 && n == 0 && m == 0 && p == 0 #1000
            mincase = min(U,mRNA)
            t += mincase
            if mincase== U
                s = 0
            else
                n += 1
            end
            ##
        elseif s == 1 && n == 0 && m == 0 & p >= 1 #1001
            mincase = min(Pd,U,mRNA)
            t += mincase
            if mincase == U
                s = 0
            elseif mincase == mRNA
                n += 1
            else
                p += -1
            end
            ##
        elseif s == 0 && n >= 1 && m == 0 && p == 0#0100
            mincase = min(L,MatmRNA)
            t += mincase
            if mincase == L
                s = 1
            else
                m += 1
                n += -1
            end
            ##
        elseif s == 0 && n >= 1 && m == 0 && p >= 1 #0101
            mincase = min(MatmRNA,L,Pd)
            t += mincase
            if mincase == L
                s = 1
            elseif mincase == Pd
                p += -1
            else
                m += 1
                n += -1
            end
            ##
        elseif s == 0 && n >= 1 && m == 0 && p == 0#0100
            mincase = min(L,MatmRNA)
            t += mincase
            if mincase == L
                s = 1
            else
                m += 1
                n += -1
            end
            ##
        elseif s == 0 && n >= 1 && m == 0 && p >= 1 #0101
            mincase = min(MatmRNA,L,Pd)
            t += mincase
            if mincase == L
                s = 1
            elseif mincase == Pd
                p += -1
            else
                m += 1
                n += -1
            end
            ##
        elseif s == 0 && n >= 1 && m == 0 && p == 0#0100
            mincase = min(L,MatmRNA)
            t += mincase
            if mincase == L
                s = 1
            else
                m += 1
                n += -1
            end
        elseif s == 0 && n >= 1 && m == 0 && p >= 1 #0101
            mincase = min(MatmRNA,L,Pd)
            t += mincase
            if mincase == L
                s = 1
            elseif mincase == Pd
                p += -1
            else
                m += 1
                n += -1
            end
            ##
        elseif s == 1 && n >= 1 && m == 0 && p == 0 #1100.
            mincase = min(U,MatmRNA,mRNA)
            t += mincase
            if mincase == U
                s = 0
            elseif mincase == MatmRNA
                m += 1
                n += -1
            else
                n += 1
            end
            ##
        elseif s == 1 && n >= 1 && m == 0 && p >= 1 #1101.
            mincase = min(U,MatmRNA,mRNA,Pd)
            t += mincase
            if mincase == U
                s = 0
            elseif mincase == MatmRNA
                m += 1
                n += -1
            elseif mincase == mRNA
                n += 1
            else
                p += -1
            end
            ##
        elseif s == 0 && n == 0 && m >= 1 && p == 0 #0010
            mincase = min(L,Md,makePr)
            t += mincase
            if mincase == L
                s = 1
            elseif mincase == makePr
                p += 1
            else
                m += -1
            end
            ##
        elseif s == 0 && n == 0 && m >= 1  && p >= 1#0011
            mincase = min(L,Md,Pd,makePr)
            t += mincase
            if mincase == L
                s = 1
            elseif mincase == Md
                m += -1
            elseif mincase == makePr
                p += 1
            else
                p += -1
            end
            ##
        elseif s == 1 && n == 0 && m >= 1 && p == 0 #1010
            mincase = min(U,mRNA,Md,makePr)
            t += mincase
            if mincase == mRNA
                n += 1
            elseif mincase == U
                s = 0
            elseif mincase == makePr
                p += 1
            else
                m += -1
            end
            ##
        elseif s == 1 && n == 0 && m >= 1 && p >= 1 #1011
            mincase = min(U,mRNA,Md,Pd,makePr)
            t += mincase
            if mincase == mRNA
                n += 1
            elseif mincase == U
                s = 0
            elseif mincase  == Pd
                p += -1
            elseif mincase == makePr
                p += 1
            else
                m += -1
            end
            ##
        elseif s == 0 && n >= 1 && m >= 1 && p == 0 #0110
            mincase = min(L,MatmRNA,Md,makePr)
            t += mincase
            if mincase == L
                s = 1
            elseif mincase == MatmRNA
                m += 1
                n += -1
            elseif mincase == makePr
                p += 1
            else
                m += -1
            end
            ##
        elseif s == 0 && n >= 1 && m >= 1 && p >= 1 #0111
            mincase = min(L,MatmRNA,Md,Pd,makePr)
            t += mincase
            if mincase == L
                s = 1
            elseif mincase == MatmRNA
                m += 1
                n += -1
            elseif mincase == Md
                m += -1
            elseif mincase == makePr
                p += 1
            else
                p += -1
            end
            ##
        elseif s == 1 && n >= 1 && m >= 1 && p == 0 #1110
            mincase = min(U,mRNA,MatmRNA,Md,makePr)
            t += mincase
            if mincase == U
                s = 0
            elseif mincase == mRNA
                n += 1
            elseif mincase == MatmRNA
                n += -1
                m += 1
            elseif mincase == makePr
                p += 1
            else
                m += -1
            end
            ##
        else s == 1 && n >= 1 && m >= 1 && p >= 1#1111
            mincase = min(U,mRNA,MatmRNA,Md,Pd,makePr)
            t += mincase
            if mincase == U
                s = 0
            elseif mincase == mRNA
                n += 1
            elseif mincase == MatmRNA
                n += -1
                m += 1
            elseif mincase == Md
                m += -1
            elseif mincase == makePr
                p += 1
            else
                p += -1
            end
        end
    end
    return (n,m,p)
end


"""
Function to sample reporters from noisy systems
-> Nascent-Mature Pathway Reporter (two examples are returned)
-> Mature-Protein Pathway Reporter (two examples are returned)
-> Nascent-Protein Pathway Reporter (two examples are returned)
-> Nascent-Nascent dual reporter (one is returned)
-> Mature-Mature dual reporter (one is returned)
-> Protein-Protein dual reporter (one is returned)
a1,b1 determine Beta parameter for KN
a2,b2 determine Beta parameters for KM, and KP
a3,b3 determine Beta parameters for deltap
"""


function reporters(lambda,mu,KN,KM,KP,deltap,a1,b1,a2,b2,a3,b3,t_max,samplemax)
    nasc1 = []
    nasc2 = []
    mRNA1 = []
    mRNA2 = []
    Protein1 = []
    Protein2 = []

    A1val = (a1+b1)/a1
    A2val = (a2+b2)/a2
    A3val = (a3+b3)/a3
    for i in 1:samplemax
        L = lambda# * A1val * rand(Beta(a1,b1))
        #Remove comment to place Beta(a1,b1) noise on lambda
        M = mu #* A1val * rand(Beta(a1,b1))
        #Remove comment to place Beta(a1,b1) noise on lambda
        Kn = KN * A1val * rand(Beta(a1,b1))
        #This samples from a scaled Beta(a,b),
        #with scaling to achieve given mean value KN
        Km = KM * A2val * rand(Beta(a2,b2))
        Kp = KP * A2val * rand(Beta(a2,b2))
        DELTAP = deltap * A3val * rand(Beta(a3,b3))
        GENE1 = copynumbers(L,M,Kn,Km,Kp,DELTAP,t_max)
        GENE2 = copynumbers(L,M,Kn,Km,Kp,DELTAP,t_max)#max(300,20*(1/lambda + 1/mu)))
        push!(nasc1, GENE1[1]) #Nascent mRNA from gene 1
        push!(nasc2, GENE2[1]) #Nascent mRNA from gene 2
        push!(mRNA1, GENE1[2]) #Mature mRNA from gene 1
        push!(mRNA2, GENE2[2]) #Mature mRNA from gene 1
        push!(Protein1, GENE1[3]) #Protein from gene 1
        push!(Protein2, GENE2[3]) #Protein from gene 2
    end
    exp_nasc1 = mean(nasc1)
    exp_nasc2 = mean(nasc2)
    exp_mRNA1 = mean(mRNA1)
    exp_mRNA2 = mean(mRNA2)
    exp_Prot1 = mean(Protein1)
    exp_Prot2 = mean(Protein2)
    N_M1 = exp_nasc1 * exp_mRNA1
    N_M2 = exp_nasc2 * exp_mRNA2
    M_P1 = exp_Prot1 * exp_mRNA1
    M_P2 = exp_Prot2 * exp_mRNA2
    N_P1 = exp_Prot1 * exp_nasc1
    N_P2 = exp_Prot2 * exp_nasc2
    N_N = exp_nasc1 * exp_nasc2
    M_M = exp_mRNA1 * exp_mRNA2
    P_P = exp_Prot1 * exp_Prot2
    return [exp_nasc1,exp_nasc2,exp_mRNA1,exp_mRNA2,
        exp_Prot1,exp_Prot2,cov(nasc1,mRNA1)/N_M1,
        cov(nasc2,mRNA2)/N_M1,cov(mRNA1,Protein1)/M_P1,
        cov(mRNA2,Protein2)/M_P1,cov(nasc1,Protein1)/N_P1,
        cov(nasc2,Protein2)/N_P2,cov(nasc1,nasc2)/N_N,
        cov(mRNA1,mRNA2)/M_M, cov(Protein1,Protein2)/P_P]
end


function ReporterResults(LList::Array, MList::Array,
    KNBetaList::Array{Array{Float64,1},1}, NList::Array, MatList::Array,
    PList::Array, DeltaPList::Array, tmax::Int64, samplemax::Int64, delta::Float64=1)
    ReporterList = []
    for L in LList
        for M in MList
            for Kn in KNBetaList
                for n in NList
                    for Km in MatList
                        for p in PList
                            for deltap in DeltaPList
                                    KN = KNval(L,M,Km,n)
                                    KM = Km
                                    KP = KPval(L,M,KN,delta,deltap,p)
                                    AList = [n,n*Km,p,Kn[1],Kn[2],coeff_var_sqrd(Kn[1],Kn[2])]
                                    BList = reporters(L,M,KN,KM,KP,deltap,Kn[1],Kn[2],5,6,8,6,tmax,samplemax)
                                    add_to_list = append!(AList,BList)
                                    push!(ReporterList,add_to_list)
                            end
                        end
                    end
                end
            end
        end
    end
    return ReporterList
end
