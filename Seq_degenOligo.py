import math, struct, os, re, sys, csv, itertools, time
from operator import itemgetter
#from episweep import *
#from energies import *
#from epitopes import *
from Seq_MSA import *
#from AnalysisTools import *
from stats import *
from modules import *
from Covariance import *

class Seq_degenOligo:
    
    
    def __init__(self, wt, sequences_filtered, choices, conservation, pairs, coupling):
        print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        print 'Degenerate oligonucleotides choice statistics:'
        print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        
        self.wt_seq = wt
        self.AA_possible_choices = choices
        self.msa_conservation_score = conservation
        self.possible_pairs = pairs
        self.msa_coupling_score = coupling
        self.filtered_seq_list = sequences_filtered
        
        self.thres_filtered_seq_choices = []
        self.all_tubes_dict = {}
        self.non_redundant_tubes = {}
        self.tube_group_ratios = {}
        self.nonjunk_tubes = {}
        self.aa_pos_deimm = {}
        self.CDM_tubes = {}
        self.selected_tubes = {}
        
        self.onebody_tubescores = {}
        self.onebody_tube_scoreslist = {}
        self.twobody_tubescores = {}
        self.possible_tubechoices = []
        self.possible_tubepairs = []
        self.AAs = "ACDEFGHIKLMNPQRSTVWY"
        self.CYO_flag = False
        
        ################## DO NOT EDIT THIS BLOCK ########################################################
        # Degenerate nucleotides:--------------->>>> A,B,C,D,G,H,K,M,N,R,S,T,V,W,Y
        #                                            | | | | | | | | | | | | | | |
        self.dN = [[1,0,0,1,0,1,0,1,1,1,0,0,1,1,0], [0,1,1,0,0,1,0,1,1,0,1,0,1,0,1], [0,1,0,1,1,0,1,0,1,1,1,0,1,0,0], [0,1,0,1,0,1,1,0,1,0,0,1,0,1,1]]
        self.nucleotide = ['A','B','C','D','G','H','K','M','N','R','S','T','V','W','Y']
        #
        # codons[base1][base2][base3] = AminoAcid in order, where 0 = A and 21 = Z (stop codon).
        self.codon = [[[8,11,8,11],[16,16,16,16],[14,15,14,15],[7,7,10,7]],[[13,6,13,6],[12,12,12,12],[14,14,14,14],[9,9,9,9]],[[3,2,3,2],[0,0,0,0],[5,5,5,5],[17,17,17,17]],[[21,19,21,19],[15,15,15,15],[21,1,18,1],[9,4,9,4]]]
        #
        self.AA = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X','Z']
        #
        ##################################################################################################




    ############################################################################# THRESHOLD FILTER ###
    #
    # THESE FILTERED AA CHOICES CAN IMPLICITLY RESTRICT TUBE SIZE & CHOICES
    # Options for filtering AA seq choices using the criteria of:
    #   1. Conservation score of AA at each position.
    #       Enter 0 for considering all in the list (recommended).
    #       The conservation thres is converted to log and then compared for filtering.
    #       Example: Entering 0.10 conservation threshold would be -log(0.10) = 2.303
    #   2. Fixed number of mutations allowed at each position.
    #       If mutations = 0, then only wildtype is considered.
    #       Or enter mutations=20 to include all possible choices.
    #   3. If filter = False, then nothing is processed. Output choices = Input choices.
    #
    def thresholdFilter(self, cons_thresh=0.00001, mutations=2, filter=True):
        print 'Applying threshold filter on sequences...'
        
        if filter==True:
            if cons_thresh == 0: cons_thresh = 0.00001
            cons_thresh = (-1.0) * math.log(cons_thresh)
            
            ################################
            # Conservation threshold allowed mutation filtering
            temp_filtered_choices = []
            for i in xrange(len(self.AA_possible_choices)):
                if self.msa_conservation_score[tuple(self.AA_possible_choices[i])] <= cons_thresh or self.AA_possible_choices[i][1] == self.wt_seq[self.AA_possible_choices[i][0]]:
                    temp_filtered_choices.append(self.AA_possible_choices[i])

            ################################
            # Filtering based on number of mutations allowed
            #  Iterating each position of the sequence:
            for i in xrange(len(self.wt_seq)):
                # Non-wildtype residues are assessed and added in decreasing order of their
                #  conservation score; i.e., a choice with a higher conservation score
                #  is given a higher preference.
                temp_dict = {}
                for j in xrange(len(temp_filtered_choices)):
                    if temp_filtered_choices[j][0] == i and temp_filtered_choices[j][1] != self.wt_seq[i]:
                        temp_dict[(temp_filtered_choices[j][0], temp_filtered_choices[j][1])] = self.msa_conservation_score[tuple(temp_filtered_choices[j])]
                s = list(sorted(temp_dict.items(), key=itemgetter(1)))
                temp_dict.clear()
        
                # Keep adding until number of mutations threshold is reached.
                for m in xrange(len(s)):
                    if len(temp_dict) >= mutations:
                        break
                    temp_dict[s[m][0]] = s[m][1]
                
                # Also add Wildtype to the above list (always included, no matter what!)
                for j in xrange(len(temp_filtered_choices)):
                    if temp_filtered_choices[j][0] == i and temp_filtered_choices[j][1] == self.wt_seq[i]:
                        temp_dict[(temp_filtered_choices[j][0], temp_filtered_choices[j][1])] = self.msa_conservation_score[tuple(temp_filtered_choices[j])]

                # This is the final sequence choices list
                for keys in temp_dict.iterkeys():
                    self.thres_filtered_seq_choices.append(list(keys))
            
            self.thres_filtered_seq_choices = map(list, set(map(tuple, self.thres_filtered_seq_choices)))
            self.thres_filtered_seq_choices.sort()
    
        elif filter==False:
            self.thres_filtered_seq_choices = self.AA_possible_choices
        print self.thres_filtered_seq_choices
        print 'Total AA thres_filtered_seq_choices: %d' % (len(self.thres_filtered_seq_choices))
            #print self.msa_conservation_score[tuple(self.AA_possible_choices[i])]




    ############# DO NOT EDIT THE MODULES BELOW ######################################################
    #
    # THESE TWO MODULES ARE USED FOR CONSTRUCTING TUBES USING DIFFERENT METHODS
    #
    def generateTubes(self, method = '', maxTubeSize=50):
        if maxTubeSize < 2:
            print "maxTubeSize should be 2 or greater"
            sys.exit()
        
        if method == 'Standard': self.generate_all_tubes(restrict_tube_size=maxTubeSize)
        elif method == 'CYO': self.chooseYourOwnAA(restrict_tube_size=maxTubeSize)
        else:
            print "Enter method='Standard' or method='CYO' as argument. Atleast one is required."
            print "You can also restrict tube size, e.g., maxTubeSize=2"
            sys.exit()
    #
    #
    ################# This module is used to generate tubes using the STANDARD method ################
    # DO NOT EDIT THIS MODULE
    def generate_all_tubes(self, restrict_tube_size=50):
        print "Generating tubes using Standard Method"
        for d1 in xrange(15):
            for d2 in xrange(15):
                for d3 in xrange(15):
                    self.degenerateOligo(d1, d2, d3, restrict_tube_size)
        #print len(self.all_tubes_dict)
    def degenerateOligo(self, dn1, dn2, dn3, restrict_tube_size):
        for b1 in xrange(4):
            if self.dN[b1][dn1] == 1:
                for b2 in xrange(4):
                    if self.dN[b2][dn2] == 1:
                        for b3 in xrange(4):
                            if self.dN[b3][dn3] == 1:
                                nucleotide_key = ''.join([self.nucleotide[dn1], self.nucleotide[dn2], self.nucleotide[dn3]])
                                if len(self.AA[self.codon[b1][b2][b3]]) <= restrict_tube_size:
                                    if self.all_tubes_dict.has_key(nucleotide_key):
                                        self.all_tubes_dict[nucleotide_key].append(self.AA[self.codon[b1][b2][b3]])
                                    else:
                                        self.all_tubes_dict[nucleotide_key] = [self.AA[self.codon[b1][b2][b3]]]
    #
    #
    ################### This module is used to generate tubes using the CYO method ##################
    # DO NOT EDIT THIS MODULE
    def chooseYourOwnAA(self, restrict_tube_size=50):
        print "Generating tubes using Choose Your Own AA Method"
        temp_dict = {}
        for f in self.thres_filtered_seq_choices:
            if temp_dict.has_key(f[0]): temp_dict[f[0]].append(f[1])
            else: temp_dict[f[0]] = [f[1]]
        for key, values in temp_dict.iteritems():
            self.nonjunk_tubes[key] = []
            possible_tubes = [c for i in range(len(values)) for c in itertools.combinations(values, i+1)]
            for t in xrange(len(possible_tubes)):
                if len(possible_tubes[t]) == 1 and possible_tubes[t][0] == self.wt_seq[key]:
                    self.nonjunk_tubes[key].append("".join(possible_tubes[t][0]))
                    self.all_tubes_dict["".join(possible_tubes[t][0])] = list(possible_tubes[t][0])
                    self.tube_group_ratios[key, "".join(possible_tubes[t])] = [0.0]
                elif self.wt_seq[key] in list(possible_tubes[t]) and len(list(possible_tubes[t])) <= restrict_tube_size:
                    self.nonjunk_tubes[key].append("".join(possible_tubes[t]))
                    self.all_tubes_dict["".join(possible_tubes[t])] = list(possible_tubes[t])
                    self.tube_group_ratios[key, "".join(possible_tubes[t])] = [0.0]
        self.non_redundant_tubes = self.all_tubes_dict
        ################################
        # This flag is used to prevent user
        # from accessing junkFilter and removeRedundant
        # modules when using this CYO module
        self.CYO_flag = True
        ################################
        # Counting selected non-junk tubes
        total_nonjunk_tubes = 0
        for keys, values in self.nonjunk_tubes.iteritems():
            total_nonjunk_tubes += len(values)
        print 'Total nonjunk_tubes: %d' % (total_nonjunk_tubes)
        print self.nonjunk_tubes
    #
    #
    ############# DO NOT EDIT THE ABOVE MODULES ######################################################




    ############################################################################## AVOID AA FILTER ###
    #
    # Tubes are filtered based upon redundancy:
    # 1. USER PROVIDED CODON or AA: Tubes coding for stop codon are removed. Or just add more residues
    #     in the list for discarding them. Argument: avoid_aa
    # 2. PROPORTION-REDUNDANCY: Tubes coding for same types of amino acids with same
    #     proportions are redundant. As such, only one of the many is kept in the list
    #     based upon their size. Tube with larger size is removed.
    # 3. TUBE-SIZE: Tubes can also be restricted based upon the number amino acids
    #     coded for. Default is a very high number for considering tubes of all sizes.
    #
    def removeRedundant(self, avoid_aa=[], restrict_tube_size=50):

        if self.CYO_flag == True:
            print "removeRedundant() cannot be used with Choose Your Own AA"
            print "Turn off this module"
            sys.exit()
        
        print 'Removing redundant tubes...'
        ################################
        # avoid AA list 
        for key_1,value_1 in self.all_tubes_dict.iteritems():
            if len(list(set(value_1) & set(avoid_aa))) > 0 and len(value_1) > 1:
                pass
            elif len(value_1) > restrict_tube_size:
                pass
            else:
                self.non_redundant_tubes[key_1] = value_1
        #print 'self.non_redundant_tubes', len(self.non_redundant_tubes)
        
        ################################
        # Proportion redundancy check
        temp_dict = {}
        temp_dict['XXX'] = ['B']
        tube_dict = temp_dict.copy()
        for key_1,value_1 in self.non_redundant_tubes.iteritems():
            aa_coded_list_1 = value_1
            redundancy_flag = 0
            for key_2,value_2 in temp_dict.iteritems():
                aa_coded_list_2 = value_2
                if set(aa_coded_list_1) == set(aa_coded_list_2):
                    # they are redundant
                    redundancy_flag = 1
                    size_1 = len(aa_coded_list_1)
                    size_2 = len(aa_coded_list_2)
                    if size_1 > size_2:
                        tube_dict[key_2] = value_2
                        break
                    else:
                        del tube_dict[key_2]
                        tube_dict[key_1] = value_1
                        break
            if redundancy_flag == 0:
                tube_dict[key_1] = value_1
            del temp_dict
            temp_dict = tube_dict.copy()
            
        self.non_redundant_tubes = tube_dict.copy()
        del temp_dict, tube_dict
        del self.non_redundant_tubes['XXX']
        #print '\nnon_redundant_tubes length',
        #print len(self.non_redundant_tubes)
        



    ################################################################################### JUNK FILTER ###
    #
    # Tubes are filtered based upon the following criteria:
    #   1. If tubes with junk mutations should be selected or not.
    #       0 = junk not allowed
    #           OR
    #       1 = junk allowed with user provided junk percent
    #
    def junkFilter(self, junk, junk_percent_allowed = 0.34):
    
        if self.CYO_flag == True:
            print "junkFilter() cannot be used with Choose Your Own AA"
            print "Turn off this module"
            sys.exit()
        
        print 'Applying junk filter...'
        ################################
        # Dictionary keys are positions and values are tubes list
        # Junk filter is applied here
        for i in xrange(len(self.wt_seq)):
            choice_list = []
            for j in xrange(len(self.thres_filtered_seq_choices)):
                if self.thres_filtered_seq_choices[j][0] == i:
                    choice_list.append(self.thres_filtered_seq_choices[j][1])
            for keys,values in self.non_redundant_tubes.iteritems():
                degen_list = values
                match = 0.0
                for choice in choice_list:
                    match = match + degen_list.count(choice)
                if junk == 1 and self.wt_seq[i] in self.non_redundant_tubes[keys]:
                    if (len(degen_list) - match)/(len(degen_list)) <= junk_percent_allowed:
                        self.tube_group_ratios[i, keys] = [(len(degen_list) - match)/(len(degen_list))]
                        if self.nonjunk_tubes.has_key(i):
                            self.nonjunk_tubes[i].append(keys)
                        else:
                            self.nonjunk_tubes[i] = [keys]
                if junk == 0 and self.wt_seq[i] in self.non_redundant_tubes[keys]:
                    if (len(degen_list) - match)/(len(degen_list)) == 0.0:
                        self.tube_group_ratios[i, keys] = [(len(degen_list) - match)/(len(degen_list))]
                        if self.nonjunk_tubes.has_key(i):
                            self.nonjunk_tubes[i].append(keys)
                        else:
                            self.nonjunk_tubes[i] = [keys]
            del choice_list
            
        ################################
        # Counting selected non-junk tubes
        total_nonjunk_tubes = 0
        for keys, values in self.nonjunk_tubes.iteritems():
            total_nonjunk_tubes += len(values)
        print 'Total nonjunk_tubes: %d' % (total_nonjunk_tubes)
        print self.nonjunk_tubes




    ################################################################################################ RENAME TUBES ###
    #
    # This is the final step for preprocessing tubes.
    # Renaming tubes so that they can be used further for calculating potentials.
    #
    def renameTubes(self):
        ################################
        # Tubes coding for wildtype residue are replaced by single letter tube code
        # The choices are stored in variable self.possible_tubechoices
        total_selected_tubes = 0
        for keys, values in self.nonjunk_tubes.iteritems():
            total_selected_tubes += len(values)
            for i in xrange(0, len(values)):
                #print self.non_redundant_tubes[values[i]], self.wt_seq[keys]
                if self.non_redundant_tubes[values[i]][0] == self.wt_seq[keys] and len(self.non_redundant_tubes[values[i]]) == 1:
                    self.possible_tubechoices.append([keys, self.non_redundant_tubes[values[i]][0]])
                else:
                    self.possible_tubechoices.append([keys, values[i]])
        print 'Total singleton TUBE choices: %d' % (total_selected_tubes)
        #print 'Number of possible_tubechoices: %d' % (len(self.possible_tubechoices))

        ################################
        # Final list of selected tubes stored as dictionary
        # Keys = Residue positions, Values = Tube choices
        for i in xrange(len(self.possible_tubechoices)):
            if self.selected_tubes.has_key(self.possible_tubechoices[i][0]):
                self.selected_tubes[self.possible_tubechoices[i][0]].append(self.possible_tubechoices[i][1])
            else:
                self.selected_tubes[self.possible_tubechoices[i][0]] = [self.possible_tubechoices[i][1]]
        #print self.selected_tubes, len(self.selected_tubes)
        #for key, value in self.selected_tubes.iteritems():
        #    print key, value
        #    for v in value:
        #        if len(v)>1: print v, self.non_redundant_tubes[v]
        #        else: print v, v




    ################################################################################## ONEBODY SCORES ###
    #
    # Calculating onebody tubescores
    # MSA conservation score is required for calculation - see MSA.py
    #
    def onebody_tube(self, normalization = True):
        print 'Implementing onebody tube...'
        for keys,values in self.selected_tubes.iteritems():
            for i in xrange(len(values)):
                if values[i] == self.wt_seq[keys] and len(values[i]) == 1:
                    aa_coded = [self.wt_seq[keys]]
                else:
                    aa_coded = self.non_redundant_tubes[values[i]]
                aa_count = {}
                for m in set(aa_coded): aa_count[m] = float(aa_coded.count(m))
                tube_score = 0.0
                tube_score_list = []
                for aa,count in aa_count.iteritems():
                    if normalization == False: norm = 1
                    else: norm = len(aa_coded)
                    
                    if (keys, aa) in self.msa_conservation_score:
                        cons = self.msa_conservation_score[(keys, aa)]/norm
                    else:
                        cons = (1*self.penalty)/norm

                    tube_score += cons
                    for nt in xrange(int(count)): tube_score_list.append(cons*norm)
                self.onebody_tubescores[keys, values[i]] = tube_score
                self.onebody_tube_scoreslist[keys, values[i]] = tube_score_list
                #var = Covariance()
                #print aa_count, keys, values[i], tube_score_list, aa_coded, tube_score, var.variance(tube_score_list)
                #print aa_count, keys, values[i], tube_score_list, aa_coded, tube_score
                del aa_count
        
        #print self.onebody_tubescores
        print 'Length of onebody_tubescores: %d' % (len(self.onebody_tubescores))




    ################################################################################## TWOBODY SCORES ###
    #
    # Calculating twobody tubescores
    # MSA conservation and coupling scores are required for calculation - see MSA.py
    # Also required is the possible_pairs list from MSA.py obtained from Chisquare test
    # The scores can be normalized based upon their size.
    #
    def twobody_tube(self, normalization = True):
        print 'Implementing twobody tube...'
        count_zero_cov = 0
        rank1_list = []; rank2_list = []; rank3_list = []; rank4_list = []
        for keys_i,values_i in self.selected_tubes.iteritems():
            for i in xrange(len(values_i)):
                if values_i[i] == self.wt_seq[keys_i] and len(values_i[i]) == 1:
                    aa_coded_i = [self.wt_seq[keys_i]]
                else:
                    aa_coded_i = self.non_redundant_tubes[values_i[i]]
                aa_count_i = {}
                for m in set(aa_coded_i): aa_count_i[m] = float(aa_coded_i.count(m))
                
                for keys_j,values_j in self.selected_tubes.iteritems():
                    # if they are in coupled list from MSA
                    if (keys_i < keys_j) and ([keys_i, keys_j] in self.possible_pairs):
                        for j in xrange(len(values_j)):
                            if values_j[j] == self.wt_seq[keys_j] and len(values_j[j]) == 1:
                                aa_coded_j = [self.wt_seq[keys_j]]
                            else:
                                aa_coded_j = self.non_redundant_tubes[values_j[j]]
                            aa_count_j = {}
                            for m in set(aa_coded_j): aa_count_j[m] = float(aa_coded_j.count(m))
    
                            tube_score = 0.0
                            for aa_i, count_i in aa_count_i.iteritems():
                                for aa_j, count_j in aa_count_j.iteritems():
                                
                                    if normalization == False: norm_a = 1; norm_b = 1
                                    else: norm_a = len(aa_coded_i); norm_b = len(aa_coded_j)

                                    if self.msa_coupling_score.has_key((keys_i, aa_i, keys_j, aa_j)):
                                        cov = (self.msa_coupling_score[(keys_i, aa_i, keys_j, aa_j)])/(norm_a*norm_b)
                                        rank1_list.append(cov)
                                    elif (keys_i, aa_i) not in self.msa_conservation_score and (keys_j, aa_j) not in self.msa_conservation_score:
                                        cov = (2*self.penalty)/(norm_a*norm_b)
                                        rank4_list.append(cov)
                                    elif (keys_i, aa_i) in self.msa_conservation_score and (keys_j, aa_j) not in self.msa_conservation_score:
                                        cov = (2*self.penalty)/(norm_a*norm_b)
                                        rank3_list.append(cov)
                                    elif (keys_i, aa_i) not in self.msa_conservation_score and (keys_j, aa_j) in self.msa_conservation_score:
                                        cov = (2*self.penalty)/(norm_a*norm_b)
                                        rank3_list.append(cov)
                                    elif (keys_i, aa_i) in self.msa_conservation_score and (keys_j, aa_j) in self.msa_conservation_score:
                                        cov = ((2*self.penalty) - self.msa_conservation_score[(keys_i, aa_i)] - self.msa_conservation_score[(keys_j, aa_j)])/(norm_a*norm_b)
                                        rank2_list.append(cov)

                                    tube_score += cov
                        
                            self.twobody_tubescores[keys_i, values_i[i], keys_j, values_j[j]] = tube_score
                            del aa_count_j
                del aa_count_i
        
        print 'Length of twobody_tubescores: %d' % (len(self.twobody_tubescores))
        #print 'count_zero_cov', count_zero_cov

        
