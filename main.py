#!/usr/bin/env python
# Haplotype estimation method for polyploid F1 populations
# Written by Ehsan Motazedi, Wageningen UR, 08-09-2017.
# Last Updated: 31-07-2018

import bamprocess
import copy
import functools
import getopt
import itertools
import os
import random as rnd
import re
import sys
import textwrap
import traceback
from bamprocess import InputError
from branchprune import BlockException, makePermutation, SetGenos, GetProbReads
from branchprune_aanvullend import Branch_Founders, Prune_Founders, Check_Genotype_Compatibility, ImputeGenotype, GetProbReads_Founders
from genotypes import Genotype, getGenotypesPop, dropHomozygous
from haplotypes import Haplotypes, getMEC, Gametogenesis
from logprob2 import loge
from math import exp, isinf
from numpy import random as nprnd
from reads import Read, Extract_Components, getReads, SubReads, SemiRead

if __name__=='__main__':
	try:
                top = False
                _max_len = len('Read/Fragment file'+'[0,1]') # used in the formatting of the help message
                hlpmsg = (['']+textwrap.wrap("Haplotype estimation and haplotype based dosage calling tool from short sequence reads for polyploid F1-populations, based on a Bayesian framework. The algorithm uses the short sequence reads and the inheritance pattern to obtain local founder haplotypes with the maximum a posteriori probability (MAP). The offspring haplotypes will be inferred from the founder haplotypes and the reads of each offspring, using the MEC criterion (Motazedi et al. 2017).", width = 95)+
                ['Written by Ehsan Motazedi, Wageningen UR, 10-09-2017.']+
                ['\nPositional Arguments:']+
                [' '*4+'Read/Fragment file'+' '*4+' STR '+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('The name of the multi-sample BAM file containing the input aligned reads of the members of the population. The first two sample names in the BAM file are considered the names of the maternal and parental samples, respectively, and the rest samples are considered as progeny.', width=64))]+
                [' '*4+'VCF file          '+' '*4+' STR '+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('The name of the multi-sample VCF file specifying the variants in the population, corresponding to the input BAM file.', width=64))]+
                [' '*4+'Output            '+' '*4+' STR '+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('The path name of the output directory to write the results, as Output/PopPolySolution_{samplename} and Output/MEC_{samplename}', width=64))]+
                ['\nOptional Arguments:']+
                [' '*4+'-a, --aggressive  '+' '*4+'[0,1]'+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('With this option, aggressive pruning will be performed with the given rate, '+u"\u03B1"+', after each extension of the founder haplotypes, in case the original pruning scheme is not able to rule-out at least (100'+u"\u00D7"+
u"\u03B1"+')% of the candidate haplotypes, so that at most 1 + 100'+u"\u00D7"+'(1-'+u"\u03B1"+')% of the founder haplotypes remain after pruning (default is the original scheme using maximum RL).', width=64))]+
                [' '*4+'-e, --error       '+' '*4+'[0,1]'+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('The fixed base calling error rate in the read-fragments, between 0 and 1. If not specified, the software will first try to use the bam quality scores for each position within each fragment. (default = 0.015).', width=64))]+
                [' '*4+'-r, --rho         '+' '*4+'[0,1]'+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('The branching threshold of the algorithm for the founder haplotypes, '+u'\u03C1'+', as explained in Motazedi et al. (2017) (default = 0.3).', width=64))]+
                [' '*4+'-k, --kappa       '+' '*4+'[0,1]'+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('The pruning rate of the algorithm, '+ u'\u03BA'+', using the relative likelihood of the founder haplotypes as explained in Motazedi et al. (2017). (default = 0.7).', width=64))]+
                [' '*4+'-w, --warmup      '+' '*4+' INT '+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('The genomic length in base-pairs from the beginning of each haplotype block, along which no candidate extension is pruned or excluded during branching. This initial warm-up phase increases the precision and efficiency of the pruning and branching steps at the cost of memeory/speed (default = 0).', width=64))]+
                [' '*4+'-v, --verror      '+' '*4+'[0,1]'+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap("The tolerable inconsistency rate of the called (or imputed/reassigned) genotypes with the final offspring haplotypes inferred from the founder haplotypes. A proportion between 0 and 1 (default = 0.4).", width=64))]+ 
                [' '*4+'--P1              '+' '*4+' STR '+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap("The name of the first parent (the mother) in the population (default first sample in the BAM header). Can be the same as the name of the second parent in case of self-fertilization. If P2 is set but not P1, P1 will be assumed the same as P2", width=64))]+
                [' '*4+'--P2              '+' '*4+' STR '+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap("The name of the second parent (the father) in the population (default second sample in the BAM header). Can be the same as the name of the first parent in case of self-fertilization. If P1 is set but not P2, P2 will be assumed the same as P1.", width=64))]+
                [' '*4+'--mmq             '+' '*4+'[0,1]'+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('Minimum read mapping quality to consider a read for phasing (default 20).', width=64))]+
                [' '*4+'--mbq             '+' '*4+'[0,1]'+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('Minimum phred score to consider a base for haplotype fragment (default 13).', width=64))]+
                [' '*4+'--maxIS           '+' '*4+'[0,1]'+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('Maximum insert-size for a paired-end read to be considered as a single fragment for phasing (default 3000).', width=64))]+
                [' '*4+'--redose          '+' '*4+'     '+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('Setting this option forces the software to ignore the dosages detected in the input VCF file and to reassign all of the dosages using maximum a posteriori (MAP) estimation. See the explanation for option --impute for details.', width=64))]+
                [' '*4+'--filter          '+' '*4+'     '+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('With this option, the phasing is only reported for the heterozygous SNPs for each individual, i.e. the homozygous SNPs of each individual are filtered out of its haplotype estimate (default is to report all SNPs).', width=64))]+
                [' '*4+'--skip            '+' '*4+'     '+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('With this option, SNPs will be ommited from the final phasing, instead of being estimated anew by default (see explanation for --redose), if no offspring extension that is compatible with the offspring genotypes can be derived from any of the candidate parental extensions.', width=64))]+
                [' '*4+'--impute          '+' '*4+'     '+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('With this option, the genotypes missing in the VCF file will be imputed in a Bayesian way, using the allele frequencies among all of the population reads to obtain the founder priors (assuming Hardy-Weinberg equilibrium), as well as the inheritance pattern and the likelihood of the reads to obtain genotypes with the maximum a posteriori probability (MAP).', width=64))]+
                [' '*4+'-t, --top         '+' '*4+'     '+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('With this option, only the most likely haplotype is reported from the final set of haplotypes that have survived the prunings (default is to report all).', width=64))]+
                ['\n'])
                optlist, args = getopt.gnu_getopt(sys.argv[1:], 'he:k:r:a:w:tv:', ["help", "error=", "kappa=", "rho=", "aggressive=", "warmup=", "top", "verror=", "P1=", "P2=", "mmq=", "mbq=", "maxIS=", "redose", "filter", "skip", "impute"])
                valid_options = set([("-a", "--aggressive"),("-h", "--help"), ("-w", "--warmup"), ("--error", '-e'), ("--kappa", '-k'), ("--rho", '-r'), ("-t","--top"), ("--mmq", ), ("--mbq", ), ("--maxIS", ), ("-v", "--verror"), ("--P1", ), ("--P2", ), ("--redose", ), ("--filter", ), ("--skip", ), ("--impute", )])	
                valid_set = set(_y for _x in valid_options for _y in _x)
                if set(_v[0] for _v in optlist).issubset(valid_set):
                        if set(["-h", "--help"]) & set(_v[0] for _v in optlist):
                                try:
                                        garbage = sys.stdout.write('{}\n'.format('\n'.join(hlpmsg)))
                                except UnicodeEncodeError as e:
                                        import codecs
                                        UTF8Writer = codecs.getwriter('utf8')
                                        sys.stdout = UTF8Writer(sys.stdout)
                                        for _x in hlpmsg:
                                                garbage = sys.stdout.write(_x+'\n')
                                finally:
                                        sys.exit(0)
                        valid_options.remove(("-h", "--help"))
                else:
                        raise InputError('Unrecognized optional arguments!')
                if len(args)<1:
                        raise InputError('No BAM and VCF file given!')
                if len(args)<2:
                        raise InputError('No VCF file given!')
                if len(args)<3:
                        raise InputError('Output file name not specified!')
                if len(args)>3:
                        raise InputError('Unrecognized positional arguments found! Check the input to HapPoly! Detected positional arguments:\n'+'\t'+'\t'.join(args))
                if set(["-t", "--top"]) & set(_v[0] for _v in optlist):
                        top = True
                valid_options.remove(("-t", "--top"))
                if set(["--redose"]) & set(_v[0] for _v in optlist):
                        GenoConstraint = False
                        sys.stderr.write("WARNING: all of the dosages will be estimated anew! The dosages given in the input VCF file will be ignored!\n")
                else:
                        GenoConstraint = True
                valid_options.remove(("--redose", ))
                if set(["--filter"]) & set(_v[0] for _v in optlist):
                        filter_on = True
                else:
                        filter_on = False
                valid_options.remove(("--filter", ))
                if set(["--impute"]) & set(_v[0] for _v in optlist):
                        if GenoConstraint:
                            sys.stderr.write("WARNING: Missing genotypes in the population will be imputed!\n")
                        Impute_Missing = True
                else:
                        Impute_Missing = False 
                valid_options.remove(("--impute", ))
                if set(["--skip"]) & set(_v[0] for _v in optlist):
                        Impute_Incompatible = False
                else:
                        Impute_Incompatible = True
                valid_options.remove(("--skip", ))
                valid_options = sorted(valid_options, key = lambda x: x[0].lstrip('-').lower())
                opts = [None, '1.5e-02', 0.7, 3000, 13, 20, None, None, 3e-01, 0.4, 0]	# default values for aggressive pruning rate, base calling error, pruning rate, maxIS, minimum base quality, minimum mapping quality, parent1, parent2, branching threshold, variant inconsistency rate and warm-up length.
                for _n, _pair in enumerate(valid_options):
                        _lst = list(_v[1] for _v in optlist if _v[0] in _pair)
                        if len(_lst)>1:
                                raise ValueError('Duplicate values for {0}/{1}!'.format(*_pair))
                        try:
				if _n in set([0, 1, 2, 8, 9]):
					opts[_n] = float(_lst[0])
				elif _n in set([3, 4, 5, 10]):
					opts[_n] = int(_lst[0])
				else:
					opts[_n] = _lst[0]
                        except IndexError: # if no value has been specified for a parameter 
                                pass
                        except ValueError as e:
                                if "could not convert string to float" in e.args[0] or "invalid literal for int" in e.args[0] or "invalid literal for float" in e.args[0]:
                                        e.args=("Invalid value '"+_lst[0]+"' passed to "+','.join(_pair)+"! Use -h, --help for help!",)+e.args[1:]
                                        raise
                                raise
                alpha, error_rate, k, maxIS, mbq, mmq, P1, P2, rho, v_error, warmup = opts
                if alpha and (alpha>1 or alpha<0):
                        raise ValueError('The aggressive pruning rate, '+u"\u03B1"+', must be between 0 and 1!')
                if isinstance(error_rate, str): # happens if no fixed error rate is specified
                        error_fixed = False
                        error_rate = float(error_rate)
                elif error_rate>1 or error_rate<0:
                        raise ValueError('The variant calling error must be between 0 and 1!')
                else:
                        error_fixed = True
                if k>1 or k<0:
                        raise ValueError('The pruning rate, '+u'\u03BA'+', must be between 0 and 1!')
                if rho>1 or rho<0:
                        raise ValueError('The branching threshold, '+u'\u03C1'+', must be between 0 and 1!')
                if maxIS < 0:
                        raise ValueError('The maximum insert-size, maxIS, cannot be negative!')
                if v_error>1 or v_error<0:
                        raise ValueError('The tolerable inconsistency rate of the offspring variants with their haplotypes must be between 0 and 1!')
                if warmup<0:
                        raise ValueError('The warm-up length cannot be negative!')
		if P1 is None and P2 is None:
			sys.stderr.write("WARNING: No sample name has been specified for the parents! The first two samples in the BAM header are automatically counted as population parents!\n")
		elif P1 is not None and P2 is None:
			sys.stderr.write("WARNING: Only P1 has been specified! P2 is assumed to be the same as P1 (self-fertilization).\n")
			P2 = P1
		elif P1 is None and P2 is not None:
			sys.stderr.write("WARNING: Only P2 has been specified! P1 is assumed to be the same as P2 (self-fertilization).\n")
			P1 = P2
		elif P1==P2:
			sys.stderr.write("WARNING: P1 and P2 are the same! This corresponds to self-pollination!\n")
                if not all(isinstance(_x, str) for _x in args):
                        raise InputError("The input/output file names must be valid strings!")
                try:
                        Frag_lsts, Qual_lsts, SampleNames = bamprocess.get_frags(args[0], args[1], maxIS, mbq, mmq)  # Extract a separate SNP-fragment/quality score list for each member of the population
                        if error_fixed: # do not use the quality scores if a fixed base calling error is given by the user
                                Qual_lsts = [None for _i in range(0, len(Frag_lsts))]
                except IOError:
                        raise InputError('The input BAM file was not found!')
                except:
                        raise
		if P1 is not None: # rearrange the order of sample data in case parents have been specified
			try:
				P1indx = SampleNames.index(P1)
			except ValueError:
				raise ValueError('The given name of the first parent (P1) does not exist in the samples!')
			try:
				P2indx = SampleNames.index(P2)
			except ValueError:
				raise ValueError('The given name of the second parent (P2) does not exist in the samples!')
			Frag_lsts = [Frag_lsts[P1indx], Frag_lsts[P2indx]] + [_x for _indx, _x in enumerate(Frag_lsts) if _indx not in (P1indx, P2indx)]
			Qual_lsts = [Qual_lsts[P1indx], Qual_lsts[P2indx]] + [_x for _indx, _x in enumerate(Qual_lsts) if _indx not in (P1indx, P2indx)]
                try:
                        GenosORIGINAL, contig_names = getGenotypesPop(args[1], SampleNames, True) # extract the genotypes of each population from the VCF file
			if P1 is not None:
				SampleNames = [SampleNames[P1indx], SampleNames[P2indx]] + [_x for _indx, _x in enumerate(SampleNames) if _indx not in (P1indx, P2indx)]
				GenosORIGINAL = [GenosORIGINAL[P1indx], GenosORIGINAL[P2indx]] + [_x for _indx, _x in enumerate(GenosORIGINAL) if _indx not in (P1indx, P2indx)]
                        if len(set(contig_names))>1:
                            raise BlockException("The VCF file contains more than one contig! Please run TriPoly separately for each contig!")
                        if set(contig_names).issubset(set(['.',''])):
                            sys.stderr.write("WARNING: No contig name is specified in the VCF file! All of the variants are assumed located on the same contig!\n")
                        ploidies_m_f = []                               
                        for _GENOS in GenosORIGINAL[0:2]: # determine the ploidy of the parents
                                for _geno in _GENOS:
                                        if set(['-']).intersection(set(_geno.GetGenes())) or len(_geno.GetGenes())<2:
                                                continue
                                        ploidies_m_f.append(len(_geno.GetGenes()))
                                        break
                        if len(ploidies_m_f)!=2:
                                raise ValueError("Could not detect ploidy levels of the parents from the VCF file or ploidy is less than 2 for some! Check the VCF file!")
                        if ploidies_m_f[0] % 2!=0:
                            sys.stderr.write('WARNING: balanced meiosis is not possible for the mother!\n') 
                            m_b = False
                        else:
                            m_b = True
                        if ploidies_m_f[1] % 2!=0:
                            sys.stderr.write('WARNING: balanced meiosis is not possible for the father!\n') 
                            f_b = False
                        else:
                            f_b = True
                        ploidies_c = []                               
                        exclude_samples = re.split('\s+', os.getenv('exclude', '').strip()) # obtain the list of the samples to be excluded from the analysis
			def get_index(lst, element):
				try:
					index_value = lst.index(element)
				except ValueError:
					garbage = (not element) or sys.stderr.write("WARNING: excluded sample '{0:s}' does not exist in the bam and vcf files!\n".format(element))
					index_value = -1
				return index_value
			exclude_indices = sorted([get_index(SampleNames, _sam) for _sam in exclude_samples])
			exclude_indices = [_indx for _indx in exclude_indices if _indx>-1]
			if any([_indx<2 for _indx in exclude_indices]):
				raise InputError('One or more parents among the specified to be excluded samples! PopPoly cannot work without parents!')
			GenosORIGINAL, Qual_lsts, SampleNames, Frag_lsts = [list(_data) for _data in zip(*[(GenosORIGINAL[_indx], Qual_lsts[_indx], SampleNames[_indx], Frag_lsts[_indx]) for _indx in range(0, len(GenosORIGINAL)) if _indx not in exclude_indices])] # exclude the desired samples
			if len(SampleNames)<3:
				raise InputError('All of the offspring have been excluded! PopPoly needs at least one!')
			if P1 is not None:
				if P1==P2:
					SampleNames[0]+="_Megasporocyte"
					SampleNames[1]+="_Microsporocyte"
                        for _GENOS in GenosORIGINAL[2:]: # determine the ploidy of the offspring
                            got_it = False
                            for _geno in _GENOS:
                                if set(['-']).intersection(set(_geno.GetGenes())) or len(_geno.GetGenes())<2:
                                    continue
                                ploidies_c.append(len(_geno.GetGenes()))
                                got_it = True
                                break
                            if not got_it: # if there is no genotype information, determine the child ploidy from its parents
                                ploidies_c.append(sum(ploidies_m_f)//2)
                        if m_b and f_b:
                            gamete_ploidies = (None, None)
                        elif m_b:
                            gamete_ploidies = (None, _ploidy_c-ploidies_m_f[0]//2) # (None, ploidies_m_f[1]//2+_ploidy_c-(ploidies_m_f[0]//2+ploidies_m_f[1]//2))
                        elif f_b:
                            gamete_ploidies = (_ploidy_c-ploidies_m_f[1]//2, None) # (ploidies_m_f[0]//2+_ploidy_c-(ploidies_m_f[0]//2+ploidies_m_f[1]//2), None))
                        else:
                            if rnd.random > 0.5:
                                gameploidies = (ploidies_m_f[0]//2, _ploidy_c-ploidies_m_f[0]//2)
                            else:
                                gamete_ploidies = (_ploidy_c-ploidies_m_f[1]//2, ploidies_m_f[1]//2)
                        ploidy_levels =  ploidies_m_f + ploidies_c
                except IOError:
                        raise InputError('The VCF file was not found!')
                except:
                        raise
		if not Impute_Missing:  # eliminate SNPs with missing genotypes in the population, if Impute is not set, from the fragments
			SNP_to_delete = []
			for _n_snp in range(0, len(GenosORIGINAL[0])):
				if not all(not _Genos[_n_snp].isMISSING() for _Genos in GenosORIGINAL[0:2]):
					SNP_to_delete.append(GenosORIGINAL[0][_n_snp].GetS())
					sys.stderr.write('WARNING: SNP {0:d}, position {1:d} is missing in one or both parents! It will be therefore skipped from phasing!\n'.format(SNP_to_delete[-1]+1, GenosORIGINAL[0][_n_snp].GetPos()))
			for _x in range(0, len(Frag_lsts)):
				_tmpLST = [_r for _r in Frag_lsts[_x]]
				for _y in range(0, len(_tmpLST)):
					_tmp =_tmpLST[_y].GetDict()
					for n_snp in SNP_to_delete:
						garbage = _tmp.pop(n_snp, None)
					_tmpLST[_y] = Read(_tmp)
				Frag_lsts[_x] = tuple(_r for _r in _tmpLST)
		if [] in GenosORIGINAL:
			sys.exit('Program terminates as there is no genotypes to phase!')
                Imputed_Genotypes = dict() # Imputed_Genotype[Poistion] = Genotype, used to adjust original genotype list later 
                Number_of_Founders = 2
		GenosORIGINAL_backup = copy.deepcopy(GenosORIGINAL) # make a backup of the list of the genotypes to be used in reporting of the output, as the latter gets modified later on
		BLOCKS, BLOCK_SUBREADS, MEC_Scores, first_block = [], [], [], True
		Frags =[] # the set of all fragments in the population
		for _frags in Frag_lsts: # aggregate all of the SNP-fragments to build the SNP-fragment matrix from. As pedigree info is used, a SNP is considered connected\
			Frags+=_frags # if it contained in at least one informative fragment, which can originate from either parent or the progeny.
		MECwritten = False
		for comp in Extract_Components(Frags):	 # phase each connected-component separately for the founders of the population
			sys.stderr.write("{0}aplotype block started at SNP {1:d}...\n".format("H" if first_block else "Next h", min(comp)+1))
			first_block = False
			Genos = []
			for _GenosORIGINAL in GenosORIGINAL:
                            Genos.append(sorted((_geno for _geno in _GenosORIGINAL if _geno.GetS() in comp), key = lambda x: x.GetS())) # extract the genotypes that belong to each component and\
			GenosORIGINAL = [list(set(_GenosORIGINAL)-set(_Genos)) for _GenosORIGINAL, _Genos in zip(GenosORIGINAL, Genos)]     # eliminate those genotypes from the total list
                        if (not GenoConstraint) or (Impute_Missing and any([Genos[_id][0].isMISSING() for _id in range(0, len(Genos))])): # Impute the first SNP in the block if it is missing in some individuals and Impute_Missing is set to true.
                            try:
                                if not GenoConstraint:
                                    ImpAlleles = ImputeGenotype(Genos[0][0].GetS(), tuple(map(lambda x: SemiRead(x, Genos[0][0].GetS()), _Frags) for _Frags in Frag_lsts), ploidy_levels, Qual_lsts)
                                else:
                                    ImpAlleles = ImputeGenotype(Genos[0][0].GetS(), tuple(map(lambda x: SemiRead(x, Genos[0][0].GetS()), _Frags) for _Frags in Frag_lsts), ploidy_levels, Qual_lsts, [Genos[_id][0] for _id in range(0, len(Genos))])
                            except BlockException as e:
                                sys.stderr.write('WARNING: '+''.join(e.args)+" Failed to {2:s} at SNP {0:d}, position {1:d}!\n".format(Genos[0][0].GetS()+1, Genos[0][0].GetPos(), "estimate dosages" if not GenoConstraint else "impute missing genotypes"))
                                if not GenoConstraint:
                                    ImpAlleles = [tuple('-' for _x in range(0, max(len(_Geno.GetGenes()) for _Geno in Genos[_id]))) for _id in range(0, len(Genos))]
                                else:
                                    ImpAlleles = [Genos[_id][0].GetGenes() for _id in range(0, len(Genos))]
                            except:
                                raise
                            Current_Position = Genos[0][0].GetPos()
                            for _n in range(0, len(GenosORIGINAL_backup[0])): # If a genotype has been imputed at some position, replace the called genotype with the imputation
                                if GenosORIGINAL_backup[0][_n].GetPos() == Current_Position:
                                    for _id in range(0, len(GenosORIGINAL_backup)):
                                        GenosORIGINAL_backup[_id][_n] = Genotype(GenosORIGINAL_backup[_id][_n].GetS(), Current_Position, *ImpAlleles[_id])
                                        Genos[_id][0] = Genotype(GenosORIGINAL_backup[_id][_n].GetS(), Current_Position, *ImpAlleles[_id])
                        H_pruned_founders = [tuple(Haplotypes(_Genos[0].GetS(), _Genos[0].GetS(), 0, makePermutation(_Genos[0], lognum=True), None, None, *_Genos[0].GetGenes()) for _Genos in Genos[0:Number_of_Founders])] # Start the haplotypes with the first SNP of the block
			for _id in range(0, len(H_pruned_founders[0])):
                            if '.' in H_pruned_founders[0][_id].GetVS():
                                H_pruned_founders[0][_id].ChangeVS(*['-' for _x in range(0, max(len(_Geno.GetGenes()) for _Geno in Genos[_id]))])
			MEC = [0]
			_s, block_start = (0, 0)
			while _s < (len(Genos[0])-1):
                            _s+=1
                            SemiFrags = []
                            SemiFrags = tuple(map(lambda x: SemiRead(x, Genos[0][_s].GetS()), _Frags) for _Frags in Frag_lsts) # Extract the semi-reads for position _s from the total list of the fragments
                            try:
                                Current_Position = Genos[0][_s].GetPos()
			        Hap_Impute=[[], None] # list to store candidate haplotype extensions at s and eventually genotype imputations at s!
			        for _H_base_founders in H_pruned_founders:
                                    _Hap_Impute = Branch_Founders(_H_base_founders, tuple(_Genos[_s] for _Genos in Genos[0:Number_of_Founders]), SemiFrags, ploidy_levels, rho if Current_Position-Genos[0][block_start].GetPos()>=warmup else 0, error_rate, Qual_lsts, [[_Genos[_s-1], _Genos[_s]] for _Genos in Genos[Number_of_Founders:]], Impute_Incompatible, Impute_Missing, not GenoConstraint)
				    Hap_Impute[0].extend(_Hap_Impute[0])	
				    Hap_Impute[1]=_Hap_Impute[1]
                                    if _Hap_Impute[1] is not None:
                                        for _n in range(0, len(GenosORIGINAL_backup[0])): # If a genotype has been imputed at some position, replace the called genotype with the imputation
                                            if GenosORIGINAL_backup[0][_n].GetPos() == Current_Position:
                                                for _id in range(0, len(GenosORIGINAL_backup)):
                                                    GenosORIGINAL_backup[_id][_n] = Genotype(GenosORIGINAL_backup[_id][_n].GetS(), Current_Position, *Hap_Impute[1][_id].GetGenes()) # If some genotypes have been imputed, adjust the initial list of the genotypes accordingly
                                                    Genos[_id][_s] = Genotype(GenosORIGINAL_backup[_id][_n].GetS(), Current_Position, *Hap_Impute[1][_id].GetGenes())
                                subfrags_current = list(SubReads(_Frags, [Genos[0][_x].GetS() for _x in range(block_start, _s+1)]) for _Frags in Frag_lsts)
                                H_pruned_founders, MEC = Prune_Founders(list(_H_founders_ext for _H_founders_ext in Hap_Impute[0]), k if Genos[0][_s].GetPos()-Genos[0][block_start].GetPos()>=warmup else 0, error_rate, alpha if Genos[0][_s].GetPos()-Genos[0][block_start].GetPos()>=warmup else None, subfrags_current, Qual_lsts, filter_on)
                            except BlockException as e:
                                garbage = sys.stderr.write("WARNING: New Block started at SNP number {:d} due to extension error:\n{}\n".format(Genos[0][_s].GetS()+1, e))
                                BLOCKS.append(H_pruned_founders)
                                BLOCK_SUBREADS.append([[_r for _r in _Reads] for _Reads in subfrags_current])
                                MEC_Scores.append(MEC)
                                block_start = _s
                                if (not GenoConstraint) or (Impute_Missing and any([Genos[_id][block_start].isMISSING() for _id in range(0, len(Genos))])): # Impute the first SNP in the block if it is missing in some individuals and Impute_Missing is set to true.
                                    try:
                                        if GenoConstraint:
                                            ImpAlleles = ImputeGenotype(Genos[0][block_start].GetS(), tuple(map(lambda x: SemiRead(x, Genos[0][block_start].GetS()), _Frags) for _Frags in Frag_lsts), ploidy_levels, Qual_lsts)
                                        else:
                                            ImpAlleles = ImputeGenotype(Genos[0][block_start].GetS(), tuple(map(lambda x: SemiRead(x, Genos[0][block_start].GetS()), _Frags) for _Frags in Frag_lsts), ploidy_levels, Qual_lsts, [Genos[_id][block_start] for _id in range(0, len(Genos))])
                                    except BlockException as e:
                                        sys.stderr.write('WARNING: '+''.join(e.args)+" Failed to {2:s} at SNP {0:d}, position {1:d}!\n".format(Genos[0][block_start].GetS()+1, Genos[0][block_start].GetPos(), "estimate dosages" if not GenoConstraint else "impute missing genotypes"))
                                        if not GenoConstraint:
                                            ImpAlleles = [tuple('-' for _x in range(0, max(len(_Geno.GetGenes()) for _Geno in Genos[_id]))) for _id in range(0, len(Genos))]
                                        else:
                                            ImpAlleles = [Genos[_id][block_start].GetGenes() for _id in range(0, len(Genos))]
                                    except:
                                        raise
                                    Current_Position = Genos[0][block_start].GetPos()
                                    for _n in range(0, len(GenosORIGINAL_backup[0])): # If a genotype has been imputed at some position, replace the called genotype with the imputation
                                        if GenosORIGINAL_backup[0][_n].GetPos() == Current_Position:
                                            for _id in range(0, len(GenosORIGINAL_backup)):
                                                GenosORIGINAL_backup[_id][_n] = Genotype(GenosORIGINAL_backup[_id][_n].GetS(), Current_Position, *ImpAlleles[_id])
                                                Genos[_id][block_start] = Genotype(GenosORIGINAL_backup[_id][_n].GetS(), Current_Position, *ImpAlleles[_id])
                                H_pruned_founders, MEC = [tuple(Haplotypes(_Genos[block_start].GetS(), _Genos[block_start].GetS(), 0, makePermutation(_Genos[block_start], lognum=True), None, None, *_Genos[block_start].GetGenes()) for _Genos in Genos[0:Number_of_Founders])], [0]
			BLOCKS.append(H_pruned_founders)
			BLOCK_SUBREADS.append([[_r for _r in _Reads] for _Reads in subfrags_current])
			MEC_Scores.append(MEC)
                Candid_Offspring = []
                for _block_number, _founder_hapblocks in enumerate(BLOCKS): # Having phased the founders using all of the reads in the population, we shall proceed with phasing each offspring using the estimated founder haplotypes
                    Candid_Offspring.append([]) # Candid_Offspring[_n] => Candid offsprings for the _n'th block
                    for _founder_hap_solutions in _founder_hapblocks:
			_parental_prob = _founder_hap_solutions[0].GetRL()
			#_parental_prob = GetProbReads_Founders(BLOCK_SUBREADS[_block_number], _founder_hap_solutions, error_rate, True, Qual_lsts) #sum(_Hp.GetRL() for _Hp in _founder_hap_solutions) 
			#for _id in range(0, len(_founder_hap_solutions)):
                        #    _founder_hap_solutions[_id].SetRL(_parental_prob) # Set the probability of each founder to P(Reads_parent|H_parent, error_rate)
			Candid_Offspring[-1].append([]) # Candid_Offspring[_n][_m] => Candid offsprings for the _m'th solution of the _n'th block
                        for _megagamete in Gametogenesis(_founder_hap_solutions[0], gamete_ploidies[0]):  # obtain and store all of the possible offspring haplotypes from the parents assuming no recombination
                            for _microgamete in Gametogenesis(_founder_hap_solutions[1], gamete_ploidies[1]):
                                Candid_Offspring[-1][-1].append(Haplotypes(_founder_hap_solutions[0].GetStart(), _founder_hap_solutions[0].GetStop(), _parental_prob, 0, None, None, *(_megagamete+_microgamete)))
                BLOCK_OFFSPRING = []
                for _block_number in range(0, len(Candid_Offspring)): # For each offspring, throw-away the candidate haplotypes that are not compatible with its genotypes
                    hap_accepted = []
                    for _solution_number in range(0, len(Candid_Offspring[_block_number])):
                        hap_accepted.append([]) # add offspring solutions per each founder solution for the block number _block_number
                        for _id in range(Number_of_Founders, len(Genos)):
                            hap_accepted[-1].append([]) # add the list of candidate haplotypes for each offspring
                            for _H in Candid_Offspring[_block_number][_solution_number]:
                                if (Check_Genotype_Compatibility(_H, Genos[_id], Number_of_Tolerable_Errors=int(v_error*(_H.GetStop()-_H.GetStart()+1)))): # Each haplotype could differ with at most "Number_of_Tolerable_Genotype_Errors" genotypes. Implicitely, we assume that founder genotypes are correct so that any discrepency between them and the offspring genotypes should be due to error in the offspring 
                                    hap_accepted[-1][-1].append(_H.GetCopy())
                            if hap_accepted[-1][-1]==[]:
                                garbage = sys.stderr.write("WARNING: No deduced haplotype was consistent with the called genotypes for offspring number {0:d} given the tolerable inconsistency rate {1:f}! All of the possible haplotypes will be considered to obtain MEC for this offspring regardless of its called genotypes!\n".format(_id-1, v_error))
                                for _H in Candid_Offspring[_block_number][_solution_number]:
                                    hap_accepted[-1][-1].append(_H.GetCopy())
			    for _H in hap_accepted[-1][-1]:
			        #_H.SetRL(_H.GetRL()+GetProbReads(BLOCK_SUBREADS[_block_number][_id], _H, error_rate, True, Qual_lsts[_id])) # Set P(H_ci|Reads_ci, Hm, Hf, error_rate) to P(Reads_ci|H_ci, error_rate)P(H_ci|Hm, Hf) #(_H.GetRL()) 
			        _H.SetRL(GetProbReads(BLOCK_SUBREADS[_block_number][_id], _H, error_rate, True, Qual_lsts[_id])) # Set P(H_ci|Reads_ci, Hm, Hf, error_rate) to P(Reads_ci|H_ci, error_rate)P(H_ci|Hm, Hf) #(_H.GetRL()) 
                            #hap_accepted[-1][-1] = sorted(hap_accepted[-1][-1], key = lambda x: (getMEC([_r for _r in BLOCK_SUBREADS[_block_number][_id] if not _r.isNULL()], x), -x.GetRL()))
                            hap_accepted[-1][-1] = [sorted(hap_accepted[-1][-1], key = lambda x: (getMEC([_r for _r in BLOCK_SUBREADS[_block_number][_id] if not _r.isNULL()], x), -x.GetRL()))[0]]
                    BLOCK_OFFSPRING.append(hap_accepted)# add solution for the offspring, block number _block_number [block][solution][off1lst, off2slst,...]
                NEW_BLOCKS = []  # append the offspring haplotypes to the already obtained founder haplotypes to get complete population haplotypes and store it in NEW_BLOCKS
                NEW_MEC_SCORES = []
                max_number_of_accepted_H_offspring = 0
                for _block_number in range(0, len(BLOCK_OFFSPRING)):
                    NEW_BLOCKS.append([])
                    NEW_MEC_SCORES.append([]) # obtain the MEC scores for each candidate founder and offspring haplotype using its own reads and store it in NEW_MEC_SCORES
                    for _solution_number in range(0, len(BLOCK_OFFSPRING[_block_number])):
                        max_number_of_accepted_H_offspring = max(len(_Hlst) for _Hlst in BLOCK_OFFSPRING[_block_number][_solution_number])
                        for _offspring_id in range(0, len(BLOCK_OFFSPRING[_block_number][_solution_number])):
                            _defficiency =  max_number_of_accepted_H_offspring - len(BLOCK_OFFSPRING[_block_number][_solution_number][_offspring_id])
                            if _defficiency > 0:
                                BLOCK_OFFSPRING[_block_number][_solution_number][_offspring_id]+= [BLOCK_OFFSPRING[_block_number][_solution_number][_offspring_id][0]]*_defficiency # make the lengths of all offspring solution lists equal by repeating their best solutions if necessary
                        for _candidate in range(0, max_number_of_accepted_H_offspring): # Consider all of the combinations possible for the offspring haplotypes
                            _all_offspring = [_Hlst[_candidate] for _Hlst in BLOCK_OFFSPRING[_block_number][_solution_number]]
                            NEW_BLOCKS[_block_number].append(tuple(_individual_haplotype.GetCopy() for _individual_haplotype in ([_founder_haplotype for _founder_haplotype in BLOCKS[_block_number][_solution_number]] + [_offspring_haplotype for _offspring_haplotype in _all_offspring])))
                            NEW_MEC_SCORES[_block_number].append(tuple(getMEC([_r for _r in _subreads if not _r.isNULL()], _individual_haplotype, filter_on) for _subreads, _individual_haplotype in zip(BLOCK_SUBREADS[_block_number], NEW_BLOCKS[_block_number][-1])))
		    for _all_Candidate_H_for_an_individual in map(list, zip(*NEW_BLOCKS[_block_number])): # Normalize the RL of each individual's solutions for each block
                       _norm = loge(sum(exp(_Candidate_H.GetRL()) for _Candidate_H in _all_Candidate_H_for_an_individual))
		       if isinf(_norm): # Happens when all solutions for a family member have zero likelihood (or very very small)
		    	    for _Candidate_H in _all_Candidate_H_for_an_individual:
		    	        _Candidate_H.SetRL(loge(1./len(_all_Candidate_H_for_an_individual)))
                       else:
		    	    for _Candidate_H in _all_Candidate_H_for_an_individual:
		    	        _Candidate_H.SetRL(_Candidate_H.GetRL()-_norm)
                BLOCKS = NEW_BLOCKS
                MEC_Scores = NEW_MEC_SCORES
                del NEW_BLOCKS, NEW_MEC_SCORES
		if filter_on and GenoConstraint:
			Genos = [dropHomozygous(_Genolst) for _Genolst in GenosORIGINAL_backup] # Restore the original genomes and throw away homozygous variants for each population member
		else:
			Genos = [_Genolst for _Genolst in GenosORIGINAL_backup] # Restore the original genomes
		OriginalS = [[] for _id in range(0, len(Genos))] # Store the original SNP numbers to be written to the final output
		for _id in range(0, len(Genos)):
			for _Geno in Genos[_id]:
				for _Org in GenosORIGINAL_backup[_id]:
					if _Org.GetPos()==_Geno.GetPos():
						OriginalS[_id]+=[_Org.GetS()]
						break
		TOP_SOLUTION_BLOCKS = []
		TOP_MEC_BLOCKS = []
		for _block_number in range(0, len(BLOCKS)):
		    TOP_SOLUTION_BLOCKS.append([])
		    TOP_MEC_BLOCKS.append([])
		    if top:
		        _parental_probs = [_solutions[0].GetRL() for _solutions in BLOCKS[_block_number]]
		        _probs_max = max(_parental_probs)
		        for _sloution_number in range(0, len(BLOCKS[_block_number])):
			    if BLOCKS[_block_number][_sloution_number][0].GetRL() == _probs_max:
			        TOP_SOLUTION_BLOCKS[-1].append(BLOCKS[_block_number][_sloution_number])
			        TOP_MEC_BLOCKS[-1].append(MEC_Scores[_block_number][_sloution_number])
		    else:
		        for _sloution_number in range(0, len(BLOCKS[_block_number])):
			    TOP_SOLUTION_BLOCKS[-1].append(BLOCKS[_block_number][_sloution_number])
			    TOP_MEC_BLOCKS[-1].append(MEC_Scores[_block_number][_sloution_number])
		    TOP_SOLUTION_BLOCKS[-1], TOP_MEC_BLOCKS[-1] = [_sol_MEC for _sol_MEC in zip(*sorted(zip(TOP_SOLUTION_BLOCKS[-1], TOP_MEC_BLOCKS[-1]), 
									key = lambda x: (-x[0][0].GetRL(), sum(x[1]))))]
		    if top:
		        _tmp = []
		        _tmp_MEC = []
		        _top_MEC = sum(TOP_MEC_BLOCKS[-1][0])
		        for _sol_num in range(0, len(TOP_SOLUTION_BLOCKS[-1])):
		            if sum(TOP_MEC_BLOCKS[-1][_sol_num]) == _top_MEC:
			        _tmp.append(TOP_SOLUTION_BLOCKS[-1][_sol_num])
			        _tmp_MEC.append(TOP_MEC_BLOCKS[-1][_sol_num])
		        _tmp = [_sol_and_MEC for _sol_and_MEC in zip(_tmp, _tmp_MEC)]
		        garbage = nprnd.shuffle(_tmp)
		        TOP_SOLUTION_BLOCKS[-1]=[_tmp[0][0]]
		        TOP_MEC_BLOCKS[-1]=[_tmp[0][1]]
#	            _id = -1 # choose the best estimate for each individual regardless of its family
#		    for _all_Candidate_H_for_an_individual in map(list, zip(*BLOCKS[_block_number])):
#	 	        _id+=1
#			_all_MECscores_for_an_individual = [_MEC[_id] for _MEC in MEC_Scores[_block_number]]
#			TOP_SOLUTION_MEC = sorted(zip(_all_Candidate_H_for_an_individual, _all_MECscores_for_an_individual), key=lambda x: (x[1], -x[0].GetRL()))
#			_tmp = []
#			_tmp_MEC = []
#			for _sol_num in range(0, len(TOP_SOLUTION_MEC)):
#			    if TOP_SOLUTION_MEC[_sol_num][0] in _tmp:
#			        _tmp[_tmp.index(TOP_SOLUTION_MEC[_sol_num][0])].SetRL(loge(exp(_tmp[_tmp.index(TOP_SOLUTION_MEC[_sol_num][0])].GetRL())+exp(TOP_SOLUTION_MEC[_sol_num][0].GetRL())))
#		            else:	
#				_tmp.append(TOP_SOLUTION_MEC[_sol_num][0])
#				_tmp_MEC.append(TOP_SOLUTION_MEC[_sol_num][1])
#			TOP_SOLUTION_BLOCKS[-1].append([_SOLUTION for _SOLUTION in _tmp])
#			TOP_MEC_BLOCKS[-1].append([_MEC for _MEC in _tmp_MEC])
#		        TOP_SOLUTION_BLOCKS[-1][_id], TOP_MEC_BLOCKS[-1][_id] = [_x for _x in zip(*sorted([(_sol, _mec) for  _sol, _mec in zip(TOP_SOLUTION_BLOCKS[-1][_id], TOP_MEC_BLOCKS[-1][_id])], key=lambda x: (-x[0].GetRL(), x[1])))] # sort the solutions based on their RL and MEC
#		        if top: # if only one solution is to be reported, select the solution with the least MEC that has the highest posterior probability for each member
#			    _max_prob =  TOP_SOLUTION_BLOCKS[-1][_id][0].GetRL()  # if multiple solution have the higehst...
#		            _min_MEC =  TOP_MEC_BLOCKS[-1][_id][0]	          # probability and the smallest MEC, choose one of them...
#			    final_solutions, final_MECs = [], []                  # at random
#			    for _solution_num in range(0, len(TOP_SOLUTION_BLOCKS[-1][_id])):
#			        if TOP_SOLUTION_BLOCKS[-1][_id][_solution_num].GetRL() == _max_prob and TOP_MEC_BLOCKS[-1][_id][_solution_num]==_min_MEC:
#			            final_solutions.append(TOP_SOLUTION_BLOCKS[-1][_id][_solution_num])
#		    		    final_MECs.append(TOP_MEC_BLOCKS[-1][_id][_solution_num])
#				else:
#				    break
#			    _tmp = [_sol_and_MEC for _sol_and_MEC in zip(final_solutions, final_MECs)]
#			    garbage = nprnd.shuffle(_tmp)
#			    TOP_SOLUTION_BLOCKS[-1][_id] = [_tmp[0][0]]
#		            TOP_MEC_BLOCKS[-1][_id] = [_tmp[0][1]]
		BLOCKS = TOP_SOLUTION_BLOCKS 
		MEC_Scores = TOP_MEC_BLOCKS
		_dir = True
		try:
			os.mkdir(args[2])
		except OSError as e:
			from errno import EEXIST as DIR_EXISTS # IF the dir already exists, just write the files in it! Otherwise, raise the error.
			if e.errno != DIR_EXISTS:
				_dir = False
				e.args = e.args[:-1]+(e.args[-1]+" '"+args[2]+"'",)
				e.args+=('Could not write the results to or make the specified output directory!',)
				raise OSError(*e.args)
		except:
			_dir = False
			raise
		finally:
			garbage = _dir and os.chdir(args[2]) # Only try to change the directory if the directory of destination exists
		BLOCKS_NEW = []  # Transpose BLOCKS for easier processing
		MEC_NEW = []
		for _block_num in range(0, len(BLOCKS)):
			BLOCKS_NEW.append([[] for _sample in SampleNames])
			MEC_NEW.append([[] for _sample in SampleNames])
			for _solution_num in range(0, len(BLOCKS[_block_num])):
				for _id in range(0, len(BLOCKS[_block_num][_solution_num])):
					BLOCKS_NEW[-1][_id].append(BLOCKS[_block_num][_solution_num][_id])
					MEC_NEW[-1][_id].append(MEC_Scores[_block_num][_solution_num][_id])
		BLOCKS = BLOCKS_NEW
		MEC_Scores = MEC_NEW
		for _id in range(0, len(SampleNames)):
			blcknum = 0
			BLOCKS_ID = [_BLOCK[_id] for _BLOCK in BLOCKS] 
			line_start = 1
			with open('PopPolySolution_'+SampleNames[_id], 'w') as outhap, open('MEC_'+SampleNames[_id], 'w') as outMEC:
				for _blck, _H_pruned in enumerate(BLOCKS_ID):
					solution_number = -1
					blckrng = []
					n_snps_blck = []
					all_Invalid_Genos = []
					dropped_solutions = 0
					for _HH_pruned in _H_pruned:
						solution_number += 1
						n_snps_blck.append(0)
						_rng_start = [_n for _n, _s in enumerate(OriginalS[_id]) if _s >_HH_pruned.GetStart()-1]
						_rng_end = [_n for _n, _s in enumerate(OriginalS[_id]) if _s >_HH_pruned.GetStop()]
						if len(_rng_start)>0:
							_rng_start = _rng_start[0]
						else:
							_rng_start = -1
						if len(_rng_end)>0:
							_rng_end = _rng_end[0]
						else:
							_rng_end = len(OriginalS[_id])
						if _rng_start >= 0:
							rng = [_num for _num in range(_rng_start, _rng_end)]
						else:
							rng = []
						blckrng.append(rng)
						if filter_on: # filter out homozygous SNPs
							Invalid_Genos = [_s for _s in range(_HH_pruned.GetStart(), _HH_pruned.GetStop()+1) if '-' in _HH_pruned.GetGenotype(_s) or len(set(_HH_pruned.GetGenotype(_s)))<2] # throw away unphased as well as homozygous SNPs
							Invalid_Genos = [_s for _s in Invalid_Genos if _s in [OriginalS[_id][_n] for _n in rng]]
						else:
							Invalid_Genos = []
						all_Invalid_Genos.append(Invalid_Genos)
						n_snps_blck[-1] = len(rng)-len(Invalid_Genos)
						if len(rng)-len(Invalid_Genos)<2:
							dropped_solutions +=1
							continue
                                                if (solution_number == 0):
                                                    garbage = outhap.write('-'*19+' BLOCK '+str(blcknum+1)+' '+'-'*19+'\n')
                                                garbage = outhap.write("{3}\t{0}\t{1}\t{2}\n".format(Genos[_id][rng[0]].GetPos(),
                                                                line_start, line_start+n_snps_blck[-1]-1, 
                                                                "Top Solution" if top else "Solution "+str(solution_number-dropped_solutions+1)))
                                                for _n in rng:
                                                    if OriginalS[_id][_n] not in Invalid_Genos:
                                                        _new_line = '{2}\t{0:d}\t{1}\n'.format(Genos[_id][_n].GetPos(), '\t'.join(str(_x[OriginalS[_id][_n]-_HH_pruned.GetStart()]) for _x in _HH_pruned.GetVS()), contig_names[0])
                                                        garbage = outhap.write(_new_line)
					if all((len(_x)-len(_Invalid_Genos))<2 for _x, _Invalid_Genos in zip(blckrng, all_Invalid_Genos)):
						line_start += max(_n_snps_blck for _n_snps_blck in n_snps_blck)
						continue
					blcknum+=1
					start_snp_pos = min(Genos[_id][blckrng[_x][0]].GetPos() for _x in range(0, len(blckrng)))
					stop_snp_pos = max(Genos[_id][blckrng[_x][-1]].GetPos() for _x in range(0, len(blckrng)))
					garbage = outMEC.write('-'*19+' BLOCK '+str(blcknum)+' '+'-'*19+'\n')
					garbage = outMEC.write('{0:d} SNPs over a genomic length of {1:d} nucleotides.\n'.format(max(_n_snps_blck for _n_snps_blck in n_snps_blck), stop_snp_pos-start_snp_pos+1))
					garbage = outMEC.write('Start SNP number, coordinate: {0:d}, {1:d}(bp)\tStop SNP number, coordinate:{2:d}, {3:d}(bp)\n'.format(line_start, start_snp_pos, line_start+max(_n_snps_blck for _n_snps_blck in n_snps_blck)-1, stop_snp_pos))
					_solution_num = 0
					line_start += max(_n_snps_blck for _n_snps_blck in n_snps_blck)
					for _k, _rng in enumerate(blckrng):
						if len(_rng)-len(all_Invalid_Genos[_k])<2:
							continue
						_solution_num+=1
                                                garbage = outMEC.write("MEC for {0}:\t{1:3d}\tRL of {0} (accounting for all population reads and haplotypes):\t{2:5.3f}\n".format(("the top solution" if top else "solution ")+str(_solution_num), MEC_Scores[_blck][_id][_k], exp(_H_pruned[_k].GetRL())))
	except getopt.GetoptError as e:
		garbage = sys.stderr.write("GetoptError: {0}\n".format(str(e)))
		sys.exit(2)
	except OSError as e:
		garbage = sys.stderr.write("OSError: {0}\n".format(str(e)))
		sys.exit(4)
	except SystemExit as e:
		sys.exit(e.args[0])
	except (TypeError, ValueError) as e:
		try:
			garbage = sys.stderr.write(repr(e).split('(')[0]+': '+e.args[0]+'\n')
		except UnicodeEncodeError:
			import unicodedata
			e.args = tuple(unicodedata.normalize('NFKD', _arg).encode('ascii','ignore') for _arg in e.args) # drop the unicode characters
			garbage = sys.stderr.write(repr(e).split('(')[0]+': '+' '.join(_x.strip() for _x in e.args[0].split(',') if _x!=' ')+'\n')	# eliminate the empty place of the dropped characters
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_tb(exc_traceback, limit=1, file=sys.stderr)
		traceback.print_exception(exc_type, exc_value, exc_traceback, limit = 20, file=sys.stderr)
		sys.exit(2)
	except InputError as e:
		garbage = sys.stderr.write(str(e))
		#exc_type, exc_value, exc_traceback = sys.exc_info()
		#traceback.print_tb(exc_traceback, limit=1, file=sys.stderr)
		#traceback.print_exception(exc_type, exc_value, exc_traceback, limit = 20, file=sys.stderr)
		sys.exit(2)
	except:
		garbage = sys.stderr.write('Unexpected error:\n')
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_tb(exc_traceback, limit=1, file=sys.stderr)
		traceback.print_exception(exc_type, exc_value, exc_traceback, limit = 4, file=sys.stderr)
		sys.exit(3)
