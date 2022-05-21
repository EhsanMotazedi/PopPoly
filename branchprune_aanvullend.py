# Written by Ehsan Motazedi, Wageningen UR, 06-09-2017.
# Last Updated: 07-12-2017

import copy
import collections
import itertools
import math
import numpy.random as nprnd
import sys
from branchprune import BlockException, makePermutation, reduce, GetLogProbH, GetProbReads
from genotypes import Genotype 
from haplotypes import Haplotypes, Hamming, getMEC, Gametogenesis
from logprob import diffVec, Hamming
from logprob2 import veclog as log, loge, GetlogProb
from math import sqrt, exp, isinf
from numpy import array, divide, vectorize, delete as npdel, argwhere as npwhere
from scipy import misc


genotype_err = 0.0001

def Branch_Founders(H, G, SReadsLST, ploidy_levels, rho, error, qscores=[None, None], G_offspring = None, Impute_Incompatible=True, Impute_Missing=True, redose=False):
        """ Branch the ordered pair of maternal and paternal haplotypes H = (Hm, Hf) to position s using the genotypes G = (Gm, Gf) at s, 
        and assign probability to each phasing extension through the semi-reads of all of the individuals in the family. For this purpose,
        an error rate, i.e. error, or the Q-scores of the reads are used to measure the likelihood of each read conditional on the haplotypes."""
        G = [_G for _G in G] # convert tuple to list
        for _id in range(0, 2):
            if G[_id].isMISSING(): # If the genotype is missing, extension will be skipped or the genotype has to be imputed!
                if not Impute_Missing: # without imputation, missing genotypes will simply be dropped from the final haplotypes
                    garbage = sys.stderr.write('WARNING: {1:s}\'s genotype is missing at s={0:d}, position {2:d}! Phasing extension will be escaped at s={0:d} for {1:s}!\n'.format(G[_id].GetS()+1, "Mother" if _id==0 else "Father", G[_id].GetPos()))
                else:
                    garbage = sys.stderr.write('WARNING: {1:s}\'s genotype is missing at s={0:d}, position {2:d}! It will be imputed anew!\n'.format(G[_id].GetS()+1, "Mother" if _id==0 else "Father", G[_id].GetPos()))
                G[_id] = Genotype(G[_id].GetS(), G[_id].GetPos(), *['-' for _homologue in H[_id].GetVS()])
            if all(r.isNULL() for r in SReadsLST[_id]):
                sys.stderr.write('WARNING: No semi-reads exist for {0:s} at SNP {1:d}, position {2:d}!\n'.format("Mother" if _id==0 else "Father", G[_id].GetS()+1, G[_id].GetPos()))
        G = tuple(_G for _G in G) # convert list back to tuple
        ProbH = H[0].GetRL() # H[0].GetRL() is the same as H[1].GetRL() as probability is assigned to a pair not to individuals.
        #for _id in range(0,2):
            #garbage = sys.stdout.write('Base Haplotype {2:s}, S= {1:d}, Pos= {3:d}:\n\t{0}\n'.format('\n\t'.join(('\t'.join(_x for _x in H[_id].GetGenotype(_pos))) for _pos in range(H[_id].GetStart(), H[_id].GetStop()+1)), G[_id].GetS(), "Mother" if _id==0 else "Father", G[_id].GetPos()))
        extend_H_branched = []
        extend_logprobs_branched = []
        uniques, priors, logrprobs, counts, Returned_Imputation = GetProbTot_Founders(H, G, G_offspring, SReadsLST, ploidy_levels, error, True, qscores, False, Impute_Incompatible, Impute_Missing, redose) # An imputation is returned if al of the parental extensions are incompatible with the offspring genotypes and Impute_Incompatible=True, if genotypes are missing and Impute_Missing = True, or if redose=True 
        if not uniques:
            return [[tuple(_H+Haplotypes(_G.GetS(), _G.GetS(), 0, 0, None, None, *['-' for _homologue in _H.GetVS()]) for _H, _G in zip(H, G))], Returned_Imputation] # skip extension if no extension has been possible
        logrprobs_adj = [_x+_y for (_x, _y) in zip(logrprobs, log(priors))] # adjust P[SR(s)|Hp, H, eps] by its prior P[Hp|H, eps] 
        _norm = max(logrprobs_adj)
        logrprobs_adj = [_x - _norm for _x in logrprobs_adj] # subtract the max log(prob) from the set of logprobs to prevent numerical underflow 
        _norm = loge(sum(exp(_x) for _x in logrprobs_adj)) 
        if isinf(_norm):
            logHpprobs = [-loge(len(logrprobs_adj)) for _x in _x in logrprobs_adj]  
        else:
            logHpprobs = [_x - _norm for _x in logrprobs_adj] # obtain p[Hp|SR(s), H, eps] by P[SR(s)|Hp, H, eps]
        myrho = loge(rho) # change rho to log scale 
        Candid_Offspring_Extensions = []
        for _n, Hcandid in enumerate(uniques): # remove duplicate extensions that occur due to presence of similar homologues in H
            #garbage = sys.stdout.write('\tCandidate Extension:\n\t    {0}\n'.format('\t'.join(str(_x) for _x in Hcandid.GetGenotype(Hcandid.GetStop()))))
            #garbage = sys.stdout.write("\t    prob={:7.19f}, logprob= {:7.19f}\n".format(2**logHpprobs[_n], logHpprobs[_n]))
            if logHpprobs[_n]>=myrho: # cut the extensions with an adjusted reads-probability lower than the threshold
                extend_H_branched.append(tuple(_Hcandid.GetCopy() for _Hcandid in Hcandid))
                extend_logprobs_branched.append(logHpprobs[_n])
                #garbage = sys.stdout.write("\t    Candidate Accepted!\n")
            else:
                #garbage = sys.stdout.write("\t    Candidate Rejected by rho!\n")
                pass 
        if not extend_H_branched:
            garbage = sys.stderr.write('WARNING: No founder extension survived the threshold at SNP {0:d}, position {1:d}!\n'.format(G[0].GetS()+1, G[0].GetPos()))
            _maxindex = logHpprobs.index(max(logHpprobs))
            extend_H_branched.append(uniques[_maxindex])
            for _n in range(0, len(extend_H_branched[-1])):
                extend_H_branched[-1][_n].SetRL(logHpprobs[_maxindex])
        for _H, _prob in zip(extend_H_branched, extend_logprobs_branched):  # Update the stored RL value of H during branching to\ 
            for _n in range(0, len(_H)):
                _H[_n].SetRL(ProbH+_prob) # Update the RL of Hp for each founder
        return [extend_H_branched, Returned_Imputation]

def Check_Balance_Genome(G_H, sample_name):
    """Checks if a genotype (single marker or haplotype genotype) contains an even number of chromosomes and hence
    can go through balanced meiosis."""
    if isinstance(G_H, Haplotypes):
        if len(G_H.GetVS()) % 2!=0:
            sys.stderr.write('WARNING: balanced meiosis is not possible for the {0:s}!\n'.format(sample_name)) 
            return False
        return True
    else:
        try:
            if len(G_H.GetGenes()) % 2!=0:
                sys.stderr.write('WARNING: balanced meiosis is not possible for the {0:s}!\n'.format(sample_name)) 
                return False
            return True
        except AttributeError as e:
            if "object has no attribute 'GetGenes'" in ("'str' object has no attribute 'GetGenes'",)[0]:
                sys.stderr.write("ERROR: invalid input to the function! Only accepts Haplotypes or Genotype objects for the first argument!\n"+
                                "Original error message:\n"+"\n".join(e.args)+"\n")

def Check_Genotype_Compatibility(H, Glst, Number_of_Tolerable_Errors=0, give_number_of_incompatibles=False):
    """Checks if a given phasing is compatible with the called unphased genotypes over the same region. Maximum number of tolerated
    differences between the called unphased genotypes and the genotypes deduced from the phasing is equal to the parameter 
    Number_of_Tolerable_Errors (default 0), and the function returns true if the number of inconsistencies is less than this value."""
    differences = 0
    for _G in Glst:
        if _G.isMISSING(): # Ignore missing genotypes while checking compatiblity
            continue
        _s = _G.GetS()
        if _s<H.GetStart() or _s>H.GetStop() or set(_G.GetGenes()).issubset(set(['.','-'])) or set(H.GetGenotype(_s)).issubset(set(['-','.'])) or Genotype(_s, _G.GetPos(), *H.GetGenotype(_s))==_G:
            pass
        else:
            differences += 1
    if give_number_of_incompatibles:
        return differences
    return differences<=Number_of_Tolerable_Errors

def GetProbReads_Founders(ReadsLST, H, eps = 0.0005, pplog = False, QualsLST = [None, None, None], getcounts=False, min_read_length = 2):
        """ Probability of a set of reads gathered from all of the family members, i.e. P[Rm,Rf,Rc1,...,Rcn|Hm, Hf, eps] = 
        Mult(Mult(P[r|Hm, Hf, eps] for r in R) for R in (Rm,Rf,Rc1,...,Rcn)), assuming independence & using GetlogProb(r, Vset, eps) (Berger et al. 2014, p. 4). 
        P[r|Hm, Hf, eps] = P[r|Hm, eps] if r in Rm,  P[r|Hm, Hf, eps] = P[r|Hf, eps] if r in Rf & P[r|Hm, Hf, eps] = 1/2*P[r|Hm, eps]+1/2*P[r|Hf, eps] if r 
        in Rci, i=1,...n. If getcounts if True, also calculate the number of reads assigned to each homologue."""
        probLST = [] # list to store the probability of R in (Rm, Rf, Rc1, ...,Rcn)
        if getcounts:
            countLST = [] # list to store the number of reads mapped to each homologue in [h1m,...hkmm, h1f,...mkff] for R in (Rm, Rf, Rc1, ...,Rcn)
        for _id, Reads in enumerate(ReadsLST):#[0:2]):
            try:
                if (_id == 0) or (_id == 1): # if the read set belongs to one of the founders
                    if QualsLST[_id]:
                        if getcounts:
                            probs, counts = list(zip(*[GetlogProb(_Read, H[_id], eps, _Qual, getcounts, min_read_length) for _Read, _Qual in zip(Reads, QualsLST[_id])]))
                        else:
                            probs = [GetlogProb(_Read, H[_id], eps, _Qual, False, min_read_length) for _Read, _Qual in zip(Reads, QualsLST[_id])]
                    else:
                        if getcounts:
                            probs, counts = list(zip(*[GetlogProb(_Read, H[_id], eps, None, getcounts, min_read_length) for _Read in Reads]))
                        else:
                            probs = [GetlogProb(_Read, H[_id], eps, None, False, min_read_length) for _Read in Reads]
                    if getcounts:
                        counts = [_count+[0 for _k in range(0, len(H[1].GetVS()))] for _count in counts] if (_id == 0) else [[0 for _k in range(0, len(H[0].GetVS()))]+_count for _count in counts]
                else: # if the read set belongs to an offspring
                    if QualsLST[_id]:
                        if getcounts:
                            probs_m, counts_m = list(zip(*[GetlogProb(_Read, H[0], eps, _Qual, getcounts, min_read_length) for _Read, _Qual in zip(Reads, QualsLST[_id])]))
                            probs_f, counts_f = list(zip(*[GetlogProb(_Read, H[1], eps, _Qual, getcounts, min_read_length) for _Read, _Qual in zip(Reads, QualsLST[_id])]))
                        else:
                            probs_m = [GetlogProb(_Read, H[0], eps, _Qual, False, min_read_length) for _Read, _Qual in zip(Reads, QualsLST[_id])]
                            probs_f = [GetlogProb(_Read, H[1], eps, _Qual, False, min_read_length) for _Read, _Qual in zip(Reads, QualsLST[_id])]
                    else:
                        if getcounts:
                            probs_m, counts_m = list(zip(*[GetlogProb(_Read, H[0], eps, None, getcounts, min_read_length) for _Read in Reads]))
                            probs_f, counts_f = list(zip(*[GetlogProb(_Read, H[1], eps, None, getcounts, min_read_length) for _Read in Reads]))
                        else:
                            probs_m = [GetlogProb(_Read, H[0], eps, None, False, min_read_length) for _Read in Reads]
                            probs_f = [GetlogProb(_Read, H[1], eps, None, False, min_read_length) for _Read in Reads]
                    if getcounts:
                        counts = [[1./2*_count for _count in _countsm+_countsf] for _countsm, _countsf in zip(counts_m, counts_f)]
                    #probs = [loge(1./2*(exp(_logprobm)+exp(_logprobf))) for _logprobm, _logprobf in zip(probs_m, probs_f)]
                    probs = [max(_logprobm,_logprobf) for _logprobm, _logprobf in zip(probs_m, probs_f)]
                    #probs = [loge(exp(max(_logprobm,_logprobf))/(exp(_logprobm)+exp(_logprobf))) for _logprobm, _logprobf in zip(probs_m, probs_f)]
                if getcounts:
                    countLST.append([sum(_counts[_n] for _counts in counts) for _n in range(0, len(H[0].GetVS())+len(H[1].GetVS()))])
                if pplog:
                    probLST.append(sum(probs))
                else:
                    probLST.append(exp(sum(probs)))
            except IndexError as e: # Error that occurs at the event that the Reads set is empty
                if "index 0 is out of bounds for axis 0 with size 0" in e.args[0]:
                    if getcounts:
                        countLST.append([0 for _h in range(0, len(H[0].GetVS())+len(H[1].GetVS()))])
                    if pplog:
                        probLST.append(0)
                    else:
                        probLST.append(1)
                else:
                        raise
        if getcounts:
            if pplog:
                return sum(probLST), reduce(lambda x, y: [_x+_y for _x, _y in zip(x,y)], countLST)
            else:
                return reduce(lambda x, y: x*y, probLST), reduce(lambda x,y: [_x+_y for _x, _y in zip(x,y)], countLST)
        else:
            if pplog:
                return sum(probLST)
            else:
                return reduce(lambda x, y: x*y, probLST)

def GetProbTot_Founders(H, G, G_offspring, ReadsLST, ploidy_levels, error_rate, plog = False, QscoresLST = None, usecounts=False, Impute_Incompatible=True, Impute_Missing=True, redose=False):
        """Determine the set of distinct extensions for a pair of founder base haplotypes H=(Hm, Hf) to position s using their genotypes at s G=(Gm, Gf), 
        calculate their prior weights and report the read-probabilities conditional on each extension. In case usecounts is True, also report the number 
        of observed reads compatible with each homologue, so that the upstream functions may set a threshold on the minimum number of reads compatible with 
        each homologue."""
        global genotype_err
        Returned_Imputation = None
        done, attempt = False, 1
        Imputation = [Genotype(_G.GetS(), _G.GetPos(), *_G.GetGenes()) for _G in G] # Parental genotypes are needed at s to extend H
        Imputation = Imputation + [[Genotype(_Gc.GetS(), _Gc.GetPos(), *_Gc.GetGenes()) for _Gc in _Glst] for _Glst in G_offspring] # Offspring genotypes are needed from SNP 1 to s to check compatiblity with H. Nevertheless, only the SNPs at s-1 and s are currently used to check compatiblity. 
        Impute_Missing = Impute_Missing and (any([_Gp.isMISSING() for _Gp in G]) or any([_Glst[-1].isMISSING() for _Glst in G_offspring])) # if no genotype is missing, Impute_Missing can be set to False
        if redose or Impute_Missing: # if redose, reassign all of the genotypes. Otherwise, if Impute_Missing is True and missing genotypes exist at s (for some of the parents or the offspring), try to impute them!
            if not redose:
                garbage = sys.stderr.write('WARNING: Missing offspring genotypes will be imputed at s={0:d}, position {1:d}!\n'.format(G[0].GetS()+1, G[0].GetPos())) 
            try:
                if not redose:
                    _Imputations = ImputeGenotype(G[0].GetS(), ReadsLST, ploidy_levels, error_rate, QscoresLST, FixedGenos = Imputation[0:2]+[_Glst[-1] for _Glst in Imputation[2:]])
                else:
                    _Imputations = ImputeGenotype(G[0].GetS(), ReadsLST, ploidy_levels, error_rate, QscoresLST, FixedGenos = None)
            except BlockException as e:
                sys.stderr.write('WARNING: '+''.join(e.args)+" Failed to {2:s} at SNP {0:d}, position {1:d}!\n".format(G[0].GetS()+1, G[0].GetPos(), "estimate dosages" if redose else "impute missing genotypes"))
                if redose:
                    _Imputations = [('-',)*ploidy_levels[0], ('-',)*ploidy_levels[1]]+[tuple('-' for _x in range(0, sum(ploidy_levels)//2)) for _Glst in G_offspring]
                else:
                    _Imputations = [_G.GetGenes() for _G in G]+[tuple('-' for _x in range(0, sum(ploidy_levels)//2)) if _Glst[-1].isMISSING() else _Glst[-1].GetGenes() for _Glst in G_offspring]
            except:
                raise
            Imputation = [Genotype(G[0].GetS(), G[0].GetPos(), *_alleles) for _alleles in _Imputations[0:2]]+[[Genotype(_Gc.GetS(), _Gc.GetPos(), *_Gc.GetGenes()) for _Gc in _Glst[:-1]]+[Genotype(_Glst[-1].GetS(), _Glst[-1].GetPos(), *_alleles)] for _Glst, _alleles in zip(G_offspring, _Imputations[2:])]
            Returned_Imputation = [Genotype(_Gp.GetS(), _Gp.GetPos(), *_alleles) for _Gp, _alleles in zip(G, _Imputations[0:2])]+[Genotype(_Glst[-1].GetS(), _Glst[-1].GetPos(), *_alleles) for _Glst, _alleles in zip(G_offspring, _Imputations[2:])] # This will be passed to the caller function to show imputationhas occured
        while not done and attempt<=2:
            perm = [makePermutation(_G) for _G in Imputation[0:2]] # (distinct) permutations of maternal and paternal genoptypes
            probs = []   # the probability of Semi Reads at s conditional on (Hp, H, eps)
            weights = [] # the prior pobability of Hp conditional on (H, eps)  
            Uniques = [] # distinct Hp's
            Counts = []  # Number of reads compatible with each homologue        
            for P1 in perm[0]: # evaluate all of the possible extensions for the base haplotypes of each parent
                for P2 in perm[1]:
                    Hp = (H[0] + P1, H[1] + P2)
                    if Hp not in Uniques:
                        Uniques.append(Hp)
                        if usecounts:
                            _prob, _counts = GetProbReads_Founders(ReadsLST, Hp, error_rate, plog, QscoresLST, True)
                            probs.append(_prob)
                            Counts.append(_counts)
                        else:
                            probs.append(GetProbReads_Founders(ReadsLST, Hp, error_rate, plog, QscoresLST))
                            Counts.append(None)
                        ploidies = [len(Hp[0].GetVS()), len(Hp[1].GetVS())]
                        Candid_Offspring_Extensions_Hp = []
                        for _megagamete in Gametogenesis(Haplotypes(Hp[0].GetStop()-1, Hp[0].GetStop(), 0, 0, None, None, *[tuple(_h[len(_h)-2:len(_h)]) for _h in Hp[0].GetVS()])): # obtain and store all of the possible offspring haplotypes from the parents assuming no recombination
                            for _microgamete in Gametogenesis(Haplotypes(Hp[1].GetStop()-1, Hp[1].GetStop(), 0, 0, None, None, *[tuple(_h[len(_h)-2:len(_h)]) for _h in Hp[1].GetVS()])): 
                                Candid_Offspring_Extensions_Hp.append(Haplotypes(Hp[0].GetStop()-1, Hp[0].GetStop(), 0, 0, None, None, *(_megagamete+_microgamete)))
                        _prior = 0 # the prior weight of a candidate founder extension
                        #_prior = 1 # the prior weight of a candidate founder extension
                        total_offspring_phasings_possible = len(Candid_Offspring_Extensions_Hp)
                        for _id in range(2, len(Imputation)):
                            #_min_number_of_errors = 2 # Mininum number of genotype incompatibilities for each offspring at s-1 and s, assumign a candidate parental extension Hp (min = 0 and max = 2, naturally!)
                            _number_of_compatible_phasings = 0
                            for _Hc in Candid_Offspring_Extensions_Hp:
                                #_error_Hc = Check_Genotype_Compatibility(_Hc, Imputation[_id], 0, give_number_of_incompatibles=True)
                                #if _error_Hc < _min_number_of_errors:
                                 #        _min_number_of_errors = _error_Hc
                                if Check_Genotype_Compatibility(_Hc, Imputation[_id], 0):
                                    _number_of_compatible_phasings +=1
                                else:
                                    pass
                            _prior+=(_number_of_compatible_phasings*_number_of_compatible_phasings)
                            #print(_prior)
                            #_prior*=misc.comb(2,_min_number_of_errors)*genotype_err**_min_number_of_errors*(1-genotype_err)**(2-_min_number_of_errors)
                        #weights.append(_prior/(1e-60+len(Candid_Offspring_Extensions_Hp)*(len(Imputation)-2))) # P(Hm,Hf|ReadLST)=P(ReadLST|Hm,Hf)P(Hm,Hf)=P(ReadsLST|Hm,Hf)P(Offspring Genotypes|Hm, Hf)
                        #if weights[-1]>1e-10:
                        #    weights[-1]=1
                        weights.append(float(_prior)/(total_offspring_phasings_possible*total_offspring_phasings_possible))
                        for _Hp in Hp:
                            _npVset = []
                            for _v in _Hp.GetVS():
                                _npv = array(_v)
                                _npVset.append(npdel(_npv, npwhere(_npv=='-')).tolist())
                            try:
                                weights[-1]*=(2**GetLogProbH(Haplotypes(1, 2, 1, loge(len(set(itertools.permutations(tuple((_v[-2],_v[-1]) for _v in _npVset))))),None, None, *tuple((_v[-2],_v[-1]) for _v in _npVset))))
                            except IndexError:
                                pass
                        #weights[-1]*=4**sum(1 for x in Hp[0].GetVS() if str(x[-1])==str(x[-2]) and str(x[-1])=='0')
                        #weights[-1]*=4**sum(1 for x in Hp[1].GetVS() if str(x[-1])==str(x[-2]) and str(x[-1])=='0')
                        #weights.append(_prior) # P(Hm,Hf|ReadLST)=P(ReadLST|Hm,Hf)P(Hm,Hf)=P(ReadsLST|Hm,Hf)P(Offspring Genotypes|Hm, Hf)
                        #if _prior>1e-10:
                        #        weights.append(1) # uninformative prior
                        #else:
                        #        weights.append(0) # incompatible extension
                    else:
                        pass
            #all_incompatible_prob = (1+1e-10)*(genotype_err**2)**(len(Imputation)-2)
            if all(_x<1e-60 for _x in weights):# if no offspring extension derived from the parental extensions is compatible with the offspring genotypes, estimate the offspring genotype at s anew. This condition is NOT expected to occur with "redose" set to True. Check if all weights are zero taking numerical uncertainty into account.
            #if all(_x<all_incompatible_prob for _x in weights):# if no offspring extension derived from the parental extensions is compatible with the offspring genotypes, estimate the offspring genotype at s anew. This condition is NOT expected to occur with "redose" set to True. Check if all weights are zero taking numerical uncertainty into account.
                if not Impute_Incompatible:
                    sys.stderr.write("WARNING: Parental genotypes were incompatible with the offspring genotypes! Extension will be skipped at SNP {0:d}, position {1:d}!\n".format(Imputation[0].GetS()+1, Imputation[0].GetPos()))
                    Uniques, weights = [], []
                    done = True
                else:
                    attempt+=1
                    if attempt<=2:
                        sys.stderr.write("WARNING: Parental genotypes were incompatible with the offspring genotypes! All of the genotypes will be imputed anew at SNP {0:d}, position {1:d}!\n".format(Imputation[0].GetS()+1, Imputation[0].GetPos()))
                        _Imputations = Imputation[0:2]+[_Glst[-1] for _Glst in Imputation[2:]]
                        try:
                            _Imputations = [Genotype(_G.GetS(), _G.GetPos(), *_alleles) for _G, _alleles in zip(_Imputations, ImputeGenotype(G[0].GetS(), ReadsLST, ploidy_levels, error_rate, QscoresLST, None))]
                        except BlockException as e:
                            sys.stderr.write('WARNING: '+''.join(e.args)+" Extension will be skipped at SNP {0:d}, position {1:d}!\n".format(Imputation[0].GetS()+1, Imputation[0].GetPos()))
                            Uniques, weights = [], []
                            done = True
                        except:
                            raise
                        else:
                            Returned_Imputation = [_G for _G in _Imputations]
                            Imputation = _Imputations[0:2] + [[Genotype(_Gc.GetS(), _Gc.GetPos(), *_Gc.GetGenes()) for _Gc in _Glst[:-1]]+[_Impute] for _Glst, _Impute in zip(G_offspring,_Imputations[2:])]
                    #weights = [1./len(weights) for _w in weights] # Uninformative prior
                    #done = True
            else:
                done = True
        if not done:
            sys.stderr.write("WARNING: Parental genotypes were still incompatible with the offspring after imputation! Extension will be therefore skipped at SNP {0:d}, position {1:d}!\n".format(Imputation[0].GetS()+1, Imputation[0].GetPos()))
            Uniques, weights = [], []
        _norm = float(sum(weights))
        weights = [_w/_norm for _w in weights]
        return Uniques, weights, probs, Counts, Returned_Imputation

def ImputeGenotype(Position, ReadLST, PloidyLevels, ErrorRate = 0.005, Quals = None, FixedGenos=None):
    """ Impute the genotypes at a specific position for a sample of reads (ReadLST) with ploidy levels given in PloidyLevels, assuming 
    Hardy-Weinberg equilibrium to obtain prior frequencies of each genotype from the allele frequencies obtained from the reads. These 
    priors are updated using the Bayes's formula to produce posterior probabilities for each parental and offspring genotype. To allow
    partial imputation, a set of genotypes (FixedGenos) can be specified (of length equal to the number of samples), so that a sample
    is not imputed if its corresponding FixedGeno is not None and not missing."""
    assert FixedGenos is None or len(FixedGenos)==len(ReadLST), "Fixed genotypes must be given as a list of genotypes equal in length to the number of samples!"
    assert FixedGenos is None or all([_G is None or isinstance(_G, Genotype) for _G in FixedGenos]), "Fixed genotypes list can only include genotype objects!"
    if FixedGenos is not None:
        if all([_G is not None and not _G.isMISSING() for _G in FixedGenos]):
            return [_G.GetGenes() for _G in FixedGenos]
    else:
        FixedGenos = [None for _Read in ReadLST]
    all_alleles = []
    for _Reads in ReadLST:
        for _r in _Reads:
            if Position in _r.GetPos():
                all_alleles.append(_r.GetDict()[Position]) # extract all of the alleles from the reads at the specified position
    n = len([_x for _x in all_alleles if (_x==1 or _x=='1' or _x==0 or _x=='0')]) # The total population coverage at the specified position
    if n == 0:
        raise BlockException("No genotype could be imputed as the population coverage was zero!")
    if not Quals or all(not _qual for _qual in Quals): # binomial model Z ~ Binomial(n, p) 
        Z = len([_x for _x in all_alleles if (_x==1 or _x=='1')])
        p_hat = float(Z)/n
        p_true_hat = (p_hat - ErrorRate)/(1 - 2 * ErrorRate) # p_hat = Z/n = (1-ErrorRate)P_true + ErrorRate(1-p_true), var_p_hat_true = ((1./(1 - 2 * ErrorRate))**2) * p_hat * (1 - p_hat) / n
    else:
        all_quals = []
        for _Reads, _Quals in zip(ReadLST, Quals):
            if _Quals is None:
                all_quals.extend([Error_Rate] * len(_Reads))
            else:
                for _r, _qual in zip(_Reads, _Quals):
                    if Position in _r.GetPos():
                        try:
                            all_quals.append(10**(-float(_qual[Position])/10)) 
                        except ValueError:
                            all_quals.append(ErrorRate)
        Z = 0 # In case quality scores, yielding error rates, are provided, we assume Prob(Allele(i)==1), i=1,...,n  has a normal distribution
        for _n in range(0, len(all_alleles)):
            if all_alleles[_n] == 1 or all_alleles[_n] =='1':
                Z += 1-all_quals[_n]
            elif all_alleles[_n] == 0 or all_alleles[_n] =='0':
                Z += all_quals[_n]
        p_true_hat = float(Z)/n
    total_genos = [] # The list of ML imputed genotypes 
    _Genotypes_Possible = dict()
    _Genotypes_Priors = dict()
    for _ploidy in set(PloidyLevels):
        _Genotypes_Possible[_ploidy] = [_g for _g in itertools.product(*[[0,1]]*_ploidy)] # all possible bi-allelic SNP genotypes 
        _Genotypes_Priors[_ploidy] = [1]*len(_Genotypes_Possible[_ploidy])
    for _ploidy, _Genos in _Genotypes_Possible.iteritems():
        for _n, _G in enumerate(_Genos):
            _number_of_ones = 0
            for _allele in _G:
                if _allele==1:
                    _number_of_ones += 1
                    _Genotypes_Priors[_ploidy][_n]*=p_true_hat
                else:
                    _Genotypes_Priors[_ploidy][_n]*=(1-p_true_hat)
            _Genotypes_Priors[_ploidy][_n] *= misc.comb(_ploidy, _number_of_ones)
    _Geno_Parents = []
    _Prior_Parents = []
    for _id in range(0, 2): # obtain first parental posteriors P(Gm,Gf|R)=P(R|Gm,Gf)P(Gm,Gf)/...
        if FixedGenos[_id] is not None and not FixedGenos[_id].isMISSING(): # in case a fixed genotype is given for a parent, just use that!
            _Geno_Parents.append((FixedGenos[_id].GetGenes(),))
            _Prior_Parents.append((1,))
        else:
            _Geno_Parents.append(tuple(_G for _G in _Genotypes_Possible[PloidyLevels[_id]]))
            _Prior_Parents.append(tuple(_P for _P in _Genotypes_Priors[PloidyLevels[_id]]))
    _Geno_Parents = [_GmGf for _GmGf in itertools.product(_Geno_Parents[0], _Geno_Parents[1])]
    _Prior_Parents =  [_PmPf for _PmPf in itertools.product(_Prior_Parents[0], _Prior_Parents[1])]
    Population_GenotypeLST = []
    for _n, _GmGf in enumerate(_Geno_Parents):
        _Child_GenotypeLST = [] # obtain the list of possible offspring genotypes for each choice of the parents
        for _megagamete in Gametogenesis(Haplotypes(Position, Position, loge(_Prior_Parents[_n][0]), 0, None, None, *_GmGf[0])):  # obtain and store all of the possible offspring haplotypes from the parents assuming no recombination
            for _microgamete in Gametogenesis(Haplotypes(Position, Position, loge(_Prior_Parents[_n][1]), 0 , None, None, *_GmGf[1])):
                _Child_GenotypeLST.append(Haplotypes(Position, Position, 0, 0, None, None, *(_megagamete+_microgamete)))
        num_all_zygotes = len(_Child_GenotypeLST)
        _Child_GenotypeLST = collections.Counter(_Child_GenotypeLST)
        for _Hc in _Child_GenotypeLST.keys():
            _Hc.SetRL(loge(float(_Child_GenotypeLST[_Hc])/num_all_zygotes))
        Population_GenotypeLST.append((Haplotypes(Position, Position, loge(_Prior_Parents[_n][0]), 0 , None, None, *_GmGf[0]), Haplotypes(Position, Position, loge(_Prior_Parents[_n][1]), 0 , None, None, *_GmGf[1]))+tuple(_Hc for _Hc in _Child_GenotypeLST.keys()))
    if len(Population_GenotypeLST) > 1: # Calculate the posterior of the parental genotype pairs only if more than one estimated genotype pair is acceptable
        Posterior_Parents = []
        for _Impute_Num in range(0, len(Population_GenotypeLST)): # Obtain the posterior probability of every population imputation. First get the posterior of each parent
            _LogLik_Reads_Parents = GetProbReads_Founders(ReadLST, (Population_GenotypeLST[_Impute_Num][0], Population_GenotypeLST[_Impute_Num][1]), eps = ErrorRate, pplog = True, QualsLST=Quals, getcounts=False, min_read_length=1) # Get the likelihood of all the reads conditional on the proposed imputation of the parents
            Posterior_Parents.append(Population_GenotypeLST[_Impute_Num][0].GetRL()+Population_GenotypeLST[_Impute_Num][1].GetRL()+_LogLik_Reads_Parents) # P(GmGf|R) = P(R|GmGf)P(GmGf)/... = P(R|GmGf)P(Gm)P(Gf)/...
        _norm_parents = 0
        for _Impute_Num in range(0, len(Population_GenotypeLST)):
            _norm_parents += exp(Posterior_Parents[_Impute_Num])
        _log_norm_parents = loge(_norm_parents)
        Posterior_Parents = [_Post - _log_norm_parents for _Post in Posterior_Parents]
        Population_GenotypeLST, Posterior_Parents  = [_Hp_Post for _Hp_Post in zip(*sorted(zip(Population_GenotypeLST, Posterior_Parents), key = lambda x: x[1], reverse = True))] # Choose the most probable estimate for parental genotypes
        for _imputation_number in range(0, len(Population_GenotypeLST)):
            for _id in range(0, 2):
                Population_GenotypeLST[_imputation_number][_id].SetRL(Posterior_Parents[_imputation_number])
        #Imputations = [Population_GenotypeLST[0][0], Population_GenotypeLST[0][1]] # Consider the most likely parental imputations as the new parental genotypes
        Imputations = [[Population_GenotypeLST[_imputation_number][0], Population_GenotypeLST[_imputation_number][1]] + [None for _cid in range(2, len(ReadLST))] for _imputation_number in range(0, len(Population_GenotypeLST))] # Consider all parental imputations as the new parental genotypes
    else:
        Imputations = [[Population_GenotypeLST[0][0], Population_GenotypeLST[0][1]]+ [None for _cid in range(2, len(ReadLST))]]
    parental_imputation_number = 0
    Max_a_posteriori = 0
    while parental_imputation_number < len(Imputations):
        for _cid in range(2, len(ReadLST)): # assume offspring independence conditional on the parents, obtain the posterior of each offspring
            if FixedGenos[_cid] is None or FixedGenos[_cid].isMISSING(): # do not impute an offspring genotype if it is given in FixedGenos
                _LogLik_Offspring = [GetProbReads(ReadLST[_cid], _Child_Geno, eps = ErrorRate, pplog = True, Quals = Quals[_cid], getcounts=False, min_read_length=1) for _Child_Geno in Population_GenotypeLST[parental_imputation_number][2:]] # Calculate the probability of reads conditional on each candidate offspring haplotype
                Posterior_Offspring = [_Hc.GetRL()+_LogLik for _Hc, _LogLik in zip(Population_GenotypeLST[parental_imputation_number][2:], _LogLik_Offspring)] # P(Gci|Ri, Gm, Gf) = P(Ri|Gci, Gm, Gf)P(Gci|Gm, Gf)/... = P(Ri|Gci)P(Gci|Gm, Gf)/...
                lognorm_Offspring = loge(sum(map(exp, Posterior_Offspring)))
                Posterior_Offspring = [_x - lognorm_Offspring for _x in Posterior_Offspring]
                Imputations[parental_imputation_number][_cid] = Population_GenotypeLST[parental_imputation_number][2+Posterior_Offspring.index(max(Posterior_Offspring))] # Choose the child genotype with maximum posterior probability
            else:
                Imputations[parental_imputation_number][_cid] = Haplotypes(Position, Position, 0, 0 , None, None, *FixedGenos[_cid].GetGenes())
        if sum(_H.GetRL() for _H in Imputations[parental_imputation_number][1:]) > sum(_H.GetRL() for _H in Imputations[Max_a_posteriori][1:]):
            Max_a_posteriori = parental_imputation_number
        parental_imputation_number+=1
    return [_H.GetVS() for _H in Imputations[Max_a_posteriori]]

def Prune_Founders(Hlst, kappa, error_rate, alpha = None, subreadsLST=[[], []], qsLST = [(), ()], Het=False):
        """ Prune a set of founder haplotypes using their log Relative Likelihoods and the pruning rate kappa (Berger et al. 2014 p. 6)."""
        try:
                #RLs = [GetProbReads_Founders(subreadsLST, _H, error_rate, True, qsLST) for _H in Hlst] #Calculate RL[H|R(SubReads), error_Rate] directly using formula (3) in Berger et al. 2014 p. 5.
                RLs = [_H[0].GetRL() for _H in Hlst]
                _norm = max(RLs)
                RLs = [_x - _norm for _x in RLs]
                _norm = loge(sum(exp(_x) for _x in RLs))
                if isinf(_norm):
                    RLs = [-loge(len(RLs)) for _x in RLs]
                else:    
                    RLs = [_x - _norm for _x in RLs]
                maxprob = max(RLs)
        except ValueError as e:
                raise BlockException(''.join(e.args)+'\n'+"Pruning could not be done! "+'\n')
        pruned = []
        k = loge(kappa)
        for _num, _H in enumerate(Hlst):
                #garbage = sys.stdout.write('****Pruning Candidate Mother {1}, RL={2:7.19f}:\n\t{0}\n'.format('\n\t'.join(('\t'.join(_x for _x in _H[0].GetGenotype(_pos))) for _pos in range(_H[0].GetStart(), _H[0].GetStop()+1)), _num, RLs[_num]))
                #garbage = sys.stdout.write('    Pruning Candidate Father {1}, RL={2:7.19f}:\n\t{0}\n-------------------\n'.format('\n\t'.join(('\t'.join(_x for _x in _H[1].GetGenotype(_pos))) for _pos in range(_H[1].GetStart(), _H[1].GetStop()+1)), _num, RLs[_num]))
                if RLs[_num] >= k + maxprob:
                        pruned.append(_H)
                        #garbage = sys.stdout.write("\t    Candidate Accepted!\n")
                        for _id in range(0, len(pruned[-1])):
                            pruned[-1][_id].SetRL(RLs[_num])
                else:
                        #garbage = sys.stdout.write("\t    Candidate Rejected!\n")
                        pass 
        non_empty_reads = [[_r for _r in subreads if not _r.isNULL()] for subreads in subreadsLST]
        MEC_Scores = [(getMEC(non_empty_reads[0], _H_pruned[0], Het)+getMEC(non_empty_reads[1], _H_pruned[1], Het)) for _H_pruned in pruned] # The MEC score is calculated for each candidate founder haplotype pair using founder reads
        HAPandMEC = sorted(zip(pruned, MEC_Scores), key = lambda x: (x[1], -1*round(x[0][0].GetRL(), 4)), reverse=False)
        if (alpha is not None) and len(pruned)>(1-alpha)*len(Hlst):
                PrunedMEC = list(zip(*HAPandMEC[0:int((1-alpha)*len(Hlst))+2])) # aggressive pruning
                return PrunedMEC[0], PrunedMEC[1]
        PrunedMEC = list(zip(*HAPandMEC))
        return PrunedMEC[0], PrunedMEC[1]
