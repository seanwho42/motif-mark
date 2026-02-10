#!/usr/bin/env python

# Author: <YOU> <optional@email.address>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "1"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNAbases = set('ATGCNatcgn')
RNAbases = set('AUGCNaucgn')

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter) - 33

def qual_score(phred_score: str) -> float:
    """
    Calculates the average quality score for a given phred_score string
    :param phred_score: str
    :return: float
    """
    score = 0
    for i, char in enumerate(phred_score):
        score += convert_phred(char)
        seq_len = i+1
    return score/seq_len

def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

def gc_content(DNA):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(DNA), "String contains invalid characters - are you sure you used a DNA sequence?"
    
    DNA = DNA.upper()
    return (DNA.count("G")+DNA.count("C"))/len(DNA)

def calc_median(lst): # PS4
    '''Given a sorted list, returns the median value of the list'''
    if len(lst) % 2 == 0:
        middle = len(lst)//2
        median = (lst[middle] + lst[middle - 1])/2
    else:
        middle = len(lst)//2
        median = lst[middle]
    # print(f'median: {median}')
    return median

def oneline_fasta(filename: str):
    '''
    Writes new fasta file (oneline_<old file name>) where sequences are only on one line each
    New filename is the old file name prepended with 'oneline_'
    '''
    oneline_filename = f'oneline_{filename}'

    with open(filename, 'r') as f, open(oneline_filename, 'w') as oneline:
        for i, line in enumerate(f):
            line = line.strip()
            if i == 0:
                oneline.write(f'{line}\n')
            elif line[0] == '>':
                oneline.write(f'\n{line}\n')
            else:
                oneline.write(line)
        oneline.write('\n')
    return oneline_filename

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")
    assert qual_score("A") == 32.0, "wrong average phred score for 'A'"
    assert qual_score("AC") == 33.0, "wrong average phred score for 'AC'"
    assert qual_score("@@##") == 16.5, "wrong average phred score for '@@##'"
    assert qual_score("EEEEAAA!") == 30.0, "wrong average phred score for 'EEEEAAA!'"
    assert qual_score("$") == 3.0, "wrong average phred score for '$'"
    print('qual_score is working')
    assert validate_base_seq('actgnACTGN'), "validate_base_seq not working for 'actgnACTGN'"
    assert not validate_base_seq('actgnACTGNxirnoap'), "validate_base_seq not working for 'actgnACTGNxirnoap'"
    assert not validate_base_seq('@*$&!()(#F)'), "validate_base_seq not working for '@*$&!()(#F)'"
    print('validate_base_seq is working')
    assert calc_median([1,2,5]) == 2, 'calc_median not working for odd numbered lists'
    assert calc_median([1,2]) == 1.5, 'calc_median not working for even numbered lists'
    assert calc_median([0,0,0,0,1,2,10]) == 0, 'calc_median not working for longer odd lists'
    assert calc_median([0,0,0,0,1,2,10,11]) == 0.5, 'calc_median not working for longer even lists'
    print('calc_median is working')
    test1, test2, test3 = oneline_fasta('test1.fa'), oneline_fasta('test2.fa'), oneline_fasta('test3.fa')
    with open(test1, 'r') as t1, open(test2, 'r') as t2, open(test3, 'r') as t3:
        t1_lines = t1.readlines()
        assert len(t1_lines) == 2, 'oneline_fasta produced wrong length for fasta with one sequence input'
        assert t1_lines[-1][-1] == '\n', 'oneline_fasta '
        assert len(t2.readlines()) == 4, 'oneline_fasta produced wrong length for fasta with multiple sequences input'
        assert len(t3.readlines()) == 2, 'oneline_fasta produced wrong length fasta which is already one line'
    print('oneline_fasta is working')
    assert gc_content('ACTCCGGAAA') == 0.5, 'gc_content not working'
    assert gc_content('actCCGGAAa') == 0.5, 'gc_content not working for mixed capitalization'
    assert gc_content('C') == 1, 'gc_content not working for fully g/c strings'
    print('gc_content is working')
