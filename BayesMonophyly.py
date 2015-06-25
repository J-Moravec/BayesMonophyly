from __future__ import division
import sys, re
import argparse as arg

"""This script will perform bayesian monophyly test on output from MrBayes or BEAST."""

def ParseArgs():
 parser = arg.ArgumentParser(prog="BayesMonophyly",description="Perform bayesian monophyly test on output from MrBayes or BEAST. " + \
                                                               "Test is simple comparison of number of trees that have monophyly " + \
                                                               "versus number of trees without monophyly, corrected by prior chance " + \
                                                               "of obtaining ratio of monophyletic/nonmonophyletic trees at random " + \
                                                               "(i.e., from data that contain no information regarding monophyly).")
 parser.add_argument("-s","--species",required=True, nargs="*", help="Species that should be monophyletic")
 parser.add_argument("-i","--input",required=True,nargs="+", help="One or more MrBayes or BEAST input files." + \
                                                                  "More files should be used only from the same analysis" + \
                                                                  "(i.e., where it actually make sense, such as standard" + \
                                                                  " two runs from MrBayes analysis).")
 parser.add_argument("-b","--burnin",required=False, default=0.2,type=int)
 args=parser.parse_args()
 return(args)

def ParseTreeFile(treefile,species):
 TEXT=[]
 try:
  FILE=open(treefile,"r")

 except IOError:
  print "ERROR: Couldn't open file, does file exists?" 
  sys.exit()
 else:
  with FILE:
   TEXT=FILE.readlines()


 #Various checks:
 #first line must be #NEXUS
 if(TEXT[0].strip("\n\t ").lower()!="#nexus"):
  print "ERROR: is file NEXUS?"
  sys.exit()
 #find begin trees
 begin_block_start=0
 for num,line in enumerate(TEXT):
  if line.strip("\n\t ").lower()=="begin trees;":
   begin_block_start=num
   break
 if begin_block_start==0:
  print "ERROR: Begin trees block not found!"
  sys.exit()
 #check if following one is translate:
 if TEXT[begin_block_start+1].strip("\n\t ").lower()!="translate":
  print "ERROR: Misformed Begin trees block, \"translate\" not found."
  sys.exit()
 #translate block, numbers from 1 to ntaxa
 #but because taxa block is not required 
 #number of taxa is not known and must be estimated from translate
 translated_taxa=dict()
 begin_block_end=0
 for num,line in enumerate(TEXT[begin_block_start+2:]):
  pair=line.strip("\n\t, ").split()
  if len(pair)!=2:
   begin_block_end=num+begin_block_start+2
   break
  else:
   translated_taxa[int(pair[0])]=pair[1]
 #check if begin_block_end has changed:
 if begin_block_end==0:
  print "ERROR: end of translation block not found."
  sys.exit()

 #test if all species are in translated_taxa
 taxa_in_tree=translated_taxa.values()
 for specie in species:
  if not specie in taxa_in_tree:
   print "ERROR: Species \"{0}\" is not in tree taxa!".format(specie)
   sys.exit()

 #now, every tree should start with "tree", so find a first tree, if not next:
 trees_start=0
 for num,line in enumerate(TEXT[begin_block_end+1:]):
  if line.strip("\n\t ")[0:4].lower()=="tree":
   trees_start=num+begin_block_end+1
   break
 #test if trees_start changed:
 if trees_start==0:
  print "ERROR: no tree was found!"
  sys.exit()

 trees=[]
 #read all trees and put them into list
 for line in TEXT[trees_start:]:
  #get tree
  line=line.strip("\n\t ;")
  if line.lower()=="end":
   #end of tree block
   break
  tree=line.split(" = ")[1].strip()
  if(tree[0:5]=="[&U] "): #remove "[&U] ", if present
   tree=tree[5:]
  #delete [&something=number] tags from BEAST
  #TODO better matching is required
  #in my file, I am currently matching only [&rate=number]
  if "[" in tree:
   tree=re.sub("\[&rate=[0-9]*\.?[0-9]*([eE]-[0-9]+)?\]","",tree)
  #get cladogram
  tree=re.sub(":[0-9]+\.?[0-9]*([eE]-[0-9]+)?","",tree)
  trees.append(tree)
 return({"species":translated_taxa,"trees":trees})

def CheckSpeciesEquivalency(list_of_dicts):
 if len(list_of_dicts)==0:
  print "INTERNAL ERROR: In CheckSpeciesEquivalency: For some reason, no dictionary of translated_taxa was passed."
  sys.exit()
 elif len(list_of_dicts)==1:
  #with only one dict, no checking for equivalence is necessary
  pass
 else:
  #test every other dict against first dict:
  template_dict=list_of_dicts[0]
  template_dict_items=list_of_dicts[0].items().sort()
  for num,matched_dict in enumerate(list_of_dicts[1:]):
   if len(template_dict)!=len(matched_dict):
    print "ERROR: Number of taxa in file 1 and file {0} is different!".format(num+2)
    sys.exit()
   
   matched_dict_items=matched_dict.items().sort()
   for template,matched in zip(template_dict_items,matched_dict_items):
    if template!=matched:
     print "ERROR: Number or taxa name differs in file 1 {1} and file {0} {2}! They are probably not equivalent!".format(num+2,str(template),str(matched))
     sys.exit()

def ete2solution(trees,translated_species):
 monophyletic_counter=0
 total_trees=len(trees)
 import ete2
 for tree in trees:
  try:
   ete2_tree=ete2.Tree(tree+";")
  except ete2.parser.newick.NewickError:
   print "INTERNAL ERROR: Problem with turning text into tree with ete2!"
   print tree
   sys.exit()

  try:
   if ete2_tree.check_monophyly(values=translated_species,target_attr="name")[0]:
    monophyletic_counter+1
  except ValueError:
   print "INTERNAL ERROR: Species are not in tree. Error in translating?"
   print translated_species
   print tree
   sys.exit()
 return(monophyletic_counter,total_trees)

def TranslateSpecies(translate_dict,species):
 reversed_dict={value:key for key,value in translate_dict.iteritems()}
 translated_species=[str(reversed_dict[item]) for item in species]
 return(translated_species)

def IterativeGrouping(tree,species):
 """Iteratively takes taxons around first species in species.
 I think that solution would be with re.sub"""
 pass #TODO

def NumberOfUnrootedTrees(n):
 import math
 return math.factorial(2*n-5)/(2^(n-3)*math.factorial(n-3))

def BayesFactor(num_monophyletic,num_total,num_taxa,num_species):
 posterior=num_monophyletic/num_total
 prior=NumberOfUnrootedTrees(num_taxa-num_species+1)/NumberOfUnrootedTrees(num_taxa)
 bayes_factor=(posterior/(1-posterior))/(prior/(1-prior))
 print "Prior probability: {0}\nPosterior probability: {1}\nBayes Factor: {2}".format(prior,posterior,bayes_factor)

if __name__ == "__main__":
 args=ParseArgs()
 #analysis makes sense only for 2 and more species:
 if len(args.species)<2:
  print "ERROR: Must specify at least 2 species."
  sys.exit()
 #also, the species must be different:
 if len(set(args.species))<len(args.species):
  print "ERROR: Please, make sure that species are unique."
  sys.exit()

 results=ParseTreeFile(args.input[0],args.species) #TODO run for every input file
 #TODO I really need to get results quickly.
 #burnin is 0.2 ergo 20%
 burnin=0.2
 trees=results["trees"]
 burn_trees=trees[int(len(trees)*burnin):]
 translated_species=TranslateSpecies(results["species"],args.species)
 (num_monophyletic,num_total)=ete2solution(burn_trees,translated_species)
 print num_monophyletic,num_total
 BayesFactor(num_monophyletic,num_total,len(results["species"]),len(translated_species))






