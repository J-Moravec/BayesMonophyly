"""This script will perform bayesian monophyly test on output from MrBayes or BEAST."""
from __future__ import division
import sys
import re
import math
import argparse as arg


def parse_args():
    parser = arg.ArgumentParser(
        prog="BayesMonophyly",
        description=(
            "Perform bayesian monophyly test on output from MrBayes"
            " or BEAST. Test is simple comparison of number of"
            " trees that have monophyly versus number of trees"
            " without monophyly, corrected by prior chance of"
            " obtaining ratio of monophyletic/nonmonophyletic"
            " trees at random (i.e., from data that contain no"
            " information regarding monophyly)."
            )
        )
    parser.add_argument(
        "-s", "--species", required=True, nargs="*",
        help="Species that should be monophyletic"
        )
    parser.add_argument(
        "-i", "--input", required=True, nargs="+",
        help=(
            "One or more MrBayes or BEAST input files."
            "More files should be used only from the same analysis"
            "(i.e., where it actually make sense, such as standard"
            " two runs from MrBayes analysis)."
            )
        )
    parser.add_argument(
        "-b", "--burnin", required=False, default=0.2,
        type=int, help="Number of trees ignored as burnin phase"
        )
    parser.add_argument(
        "-r", "--rooted", required=False, default=False,
        action="store_true",
        help=(
            "Will treat input trees as rooted"
            "(e.g. BEAST always assume some molecular clock)."
            )
        )
    args = parser.parse_args()
    return(args)


class ParsingError(Exception):
    """Unexpected string while parsing file."""
    pass


def parse_tree_file(treefile):
    """Parse posterior tree sample file from BEAST or MrBayes.

    Parse posterior tree sample file generated by MrBayes or BEAST software.
    These files are of NEXUS format and contain several block. Of these blocks,
    only the Translate and Tree blocks are of interest. Translate block contain
    list of species in trees and their translation, as names are translated
    into numbers. Tree block contains posterior sample of trees.

    This parser will search this file and returns Translate block and trees.
    There are several checks employed to ensure, that parsing is correct.

    Parameters
    ----------
    treefile : string
        path to file that is to be parsed

    Returns
    -------
    translated_taxa : dictionary
        original names of species in file and their numeric translation
    """
    tree_file_text = []
    try:
        tree_file = open(treefile,"r")
    except IOError:
        raise ParsingError("Couldn't open file, does file exists?")
    else:
        with tree_file:
            tree_file_text = tree_file.readlines()
    #Various checks:
    #first line must be #NEXUS
    if(tree_file_text[0].strip("\n\t ").lower() != "#nexus"):
        raise ParsingError("No NEXUS header. Is file NEXUS?")

    #find begin trees
    begin_block_start = 0
    for num,line in enumerate(tree_file_text):
        if line.strip("\n\t ").lower() == "begin trees;":
            begin_block_start=num
            break
    if begin_block_start == 0:
        raise ParsingError("Begin trees block not found!")

    #check if following one is translate:
    if tree_file_text[begin_block_start+1].strip("\n\t ").lower() != "translate":
        raise ParsingError(
                "ERROR: Misformed Begin trees block,"
                " \"translate\" not found."
                )

    #translate block, numbers from 1 to ntaxa
    #but because taxa block is not required 
    #number of taxa is not known and must be estimated from translate
    translated_taxa = dict()
    begin_block_end = 0
    for num,line in enumerate(tree_file_text[begin_block_start+2 : ]):
        pair = line.strip("\n\t, ").split()
        if len(pair) != 2:
            begin_block_end = num + begin_block_start + 2
            break
        else:
            translated_taxa[int(pair[0]) ] = pair[1]
    #check if begin_block_end has changed:
    if begin_block_end == 0:
        raise ParsingError("ERROR: end of translation block not found.")

    #now, every tree should start with "tree", so find a first tree, if not next:
    trees_start = 0
    for num,line in enumerate(tree_file_text[begin_block_end + 1 : ]):
        if line.strip("\n\t ")[0:4].lower() == "tree":
            trees_start = num + begin_block_end + 1
            break
    #test if trees_start changed:
    if trees_start == 0:
        raise ParsingError("ERROR: no tree was found!")
    trees = []
    #read all trees and put them into list
    for line in tree_file_text[trees_start:]:
        #get tree
        line=line.strip("\n\t ;")
        if line.lower() == "end":
            #end of tree block
            break
        tree=line.split(" = ")[1].strip()
        if(tree[0:5] == "[&U] "): #remove "[&U] ", if present
            tree = tree[5:]
        #delete [&something=number] tags from BEAST
        #TODO better matching is required
        #in my file, I am currently matching only [&rate=number]
        if "[" in tree:
            tree = re.sub("\[&rate=[0-9]*\.?[0-9]*([eE]-[0-9]+)?\]", "", tree)
        #get cladogram
        tree = re.sub(":[0-9]+\.?[0-9]*([eE]-[0-9]+)?", "", tree)
        trees.append(tree)
    return(translated_taxa, trees)


def check_species_in_taxa(species,translated_taxa):
    """Check if specified species are in dictionary.

    Simple check if species required for monophyly are in species contained
    contained in tree.

    Parameters
    ----------
    species : list of strings
        list of species for monophyly
    translated_taxa : dictionary
        dictionary of species and their number from taxa block of nexus file

    Returns
    -------
    """
    #test if all species are in translated_taxa
    taxa_in_tree = translated_taxa.values()
    for specie in species:
        if not specie in taxa_in_tree:
            raise RuntimeError("Species \"{0}\" is not in tree taxa!".format(specie))


def check_species_equivalency(list_of_dicts):
    """Check if taxa in multiple tree files are exactly the same.

    Parameters
    ----------
    list_of_dicts : list of dictionaries
        list of dictionaries of several translated_taxa from taxa block of
        multiple nexus files

    Returns
    -------
    """
    if len(list_of_dicts) == 0:
        raise RuntimeError(" For some reason, no dictionary of translated_taxa" +
                           " was passed.")
    elif len(list_of_dicts) == 1:
        #with only one dict, no checking for equivalence is necessary
        pass
    else:
        template_dict = list_of_dicts[0]
        template_dict_items=list_of_dicts[0].items()
        template_dict_items.sort()
        for num,matched_dict in enumerate(list_of_dicts[1:]):
            if len(template_dict) != len(matched_dict):
                raise RuntimeError(("Number of taxa in file 1 and " +
                                   "file {0} is different!").format(num + 2))
            matched_dict_items = matched_dict.items()
            matched_dict_items.sort()
            for template,matched in zip(template_dict_items, matched_dict_items):
                if template != matched:
                    raise RuntimeError(("Number or taxa name differs in" + 
                                        " file 1 {1} and file {0} {2}!" +
                                        " They are probably not equivalent!")
                                        .format(num + 2, str(template), str(matched))
                                        )


def ete2solution(trees, translated_species):
    """Return number of monophyletic trees for input species.

    Converts trees from string to ete2.Tree and check if specific species
    are monophyletic.

    Parameters
    ----------
    trees : list of strings
        trees, cladograms in text form
    translated_species : list of strings
        list of species for monophyly, translated into numeric form

    Returns
    -------
    monophyletic_counter : int
        number of monophyletic trees
    """
    monophyletic_counter = 0
    import ete2
    for tree in trees:
        try:
            ete2_tree=ete2.Tree(tree + ";")
        except ete2.parser.newick.NewickError:
            print tree
            raise RuntimeError("Problem with turning text into tree with ete2!")
        try:
            if ete2_tree.check_monophyly(values=translated_species,
                                         target_attr="name")[0]:
                monophyletic_counter += 1
        except ValueError:
            print translated_species
            print tree
            raise RuntimeError("Species are not in tree. Error in translating?")
    return(monophyletic_counter)

def translate_species(translate_dict, species):
    """Translate input species to numbers as they appear in tree file."""
    reversed_dict = {value : key for key,value in translate_dict.iteritems()}
    translated_species = [str(reversed_dict[item]) for item in species]
    return(translated_species)

def non_ete2solution(tree, species):
    """ I think that solution would be with re.sub"""
    pass #TODO

def n_unrooted_trees(n):
    """Returns number of unrooted trees for n taxa."""
    return math.factorial(2*n-5) / (2**(n-3) * math.factorial(n-3))


def n_rooted_trees(n):
    """Returns number of rooted trees for n taxa."""
    return math.factorial(2*n-3) / (2**(n-2) * math.factorial(n-2))


def bayes_factor(prior, posterior):
    """Compute standard Bayes factor from prior and posterior."""
    if prior in [0,1]:
        #raise some error
        pass
    elif posterior == 1:
        #raise warning?
        pass
    bayes_factor = (posterior/(1-posterior)) / (prior/(1-prior))
    return(bayes_factor)


def compute_prior(num_taxa, num_species, rooted):
    """Compute prior probability for trees, either rooted or unrooted."""
    if rooted:
        prior = (n_rooted_trees(num_taxa-num_species+1)* \
                n_rooted_trees(num_species)) / n_rooted_trees(num_taxa)
    else:
        prior = (n_unrooted_trees(num_taxa-num_species+1)* \
                n_rooted_trees(num_species)) / n_unrooted_trees(num_taxa)

    #Can't actually happen, this equation does not work for some special cases:
    #TODO
    #Question is, how to treat those values?
    #
    if prior == 1:
        raise ValueError("Prior is one.")
    if prior == 0:
        raise ValueError("Prior is zero.")
    return(prior)


def compute_posterior(num_monophyletic, num_total):
    """Compute posterior probability of monophyletic trees."""
    posterior = num_monophyletic / num_total
    #posterior 0 can happen and is legitimate
    #posterior 1 however destroys bayes factor (division by zero)
    if posterior == 1:
        raise ValueError("Posterior is one")
    return(posterior)

if __name__ == "__main__":
    args=parse_args()
    #analysis makes sense only for 2 and more species:
    if len(args.species) < 2:
        print "ERROR: Must specify at least 2 species."
        sys.exit()
    #also, the species must be different:
    if len(set(args.species)) < len(args.species):
        print "ERROR: Please, make sure that species are unique."
        sys.exit()

    #read all input files
    all_translated_taxa = []
    all_trees = []
    for input_file in args.input:
        (translated_taxa,trees) = parse_tree_file(input_file)
        check_species_in_taxa(args.species,translated_taxa)
        all_translated_taxa.append(translated_taxa)
        all_trees.append(trees)
    
    check_species_equivalency(all_translated_taxa)
    #if equivalent, every file has same species, can use the first one
    #apply burnin
    all_trees_burned = [trees[int(len(trees)*args.burnin):] for trees in all_trees]
    all_trees_burned = [inner for outer in all_trees_burned for inner in outer]
    num_total = len(all_trees_burned)
    translated_species = translate_species(all_translated_taxa[0], args.species)
    num_monophyletic = ete2solution(all_trees_burned, translated_species)
    prior = compute_prior(len(all_translated_taxa[0]),
                          len(translated_species), args.rooted)
    posterior = compute_posterior(num_monophyletic, num_total)
    bayes_factor=bayes_factor(prior, posterior)
    #output:
    all_trees_burned_num = len(all_trees_burned)
    all_trees_num = sum([len(item) for item in all_trees])
    expected_monophyletic = int(round(prior*all_trees_burned_num))
    output=("Total trees read: {0}\n"
            "Trees after burnin: {1}\n"
            "Monophyletic trees found: {5}\n"
            "Monophyletic trees expected: {3}\n"
            "(in the case of noninformative data)\n\n"
            "Prior: {2:.4f}\n"
            "Posterior: {4:.4f}\n"
            "Bayes factor: {6:.4f}\n"
            ).format(
                    all_trees_num,
                    all_trees_burned_num,
                    float(prior),
                    expected_monophyletic,
                    float(posterior),
                    num_monophyletic,
                    bayes_factor
                    )
    print output
    if bayes_factor==0:
        print("Probability of this by chance alone given prior: {0:.2e}"
              .format((1-prior)**all_trees_burned_num)
             )



