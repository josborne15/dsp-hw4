import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def count_kmers(k, seq):
    """
    Summary line: Generates the possible number of kmers of size k for a sequence

    Extended Description: Includes a section that will test to make sure k and seq
        are the correct values.If k or seq is incorrect, will return a message about
        what is wrong with the parameter.

    Parameters:
    k (int): the size of the kmer you want to calculate the number of possible kmers for

    seq (str): the sequence you want to calculate the number of possible kmers of size k for

    Return:
    if k or seq are the incorrect values:
        str: gives a message about what is incorrect about the inputed parameters
    if k and seq are correct values:
        int: gives the possible numbers of size k for seq

    """

###Checking to make sure k & seq are correct values
    if isinstance(k, str):
        return('k must be a integer')
    if (k <1):
        return('Not a valid k. k > 0')
    if seq=='':
        return('Sequence must have a least a length of 1')
    if isinstance(k, float):
        return('k must be a whole number')
    if isinstance(seq, int):
        return('Sequence must be a string')
    elif isinstance(seq, float):
        return('Sequence must be a string')
    if seq.isalpha()==False:
        return('Sequence must contain only letters.')
    if k > len(seq):
        return('k cannot be longer than the length of the sequence')
    guide=['A', 'a', 'C', 'c', 'T', 't', 'g', 'G']
    for letter in range(0, len(seq)):
        if seq[letter] not in guide:
            return('Sequence can only contain A, a, C, c, T, t, G, g as values.')


    if (len(seq) > 4**k):
        return(4**k)
    else:
        return(len(seq)-k+1)

def possible_kmers(sequence):
    """
    Summary line: Creates a pandas dataframe containing all possible k and the observed
    and expected kmers.

    Extended Description: Includes a section that will test to make sure sequence is
    the correct value. If sequence is incorrect, will return a message about
    what is wrong with the parameter. Will output a dataframe with 3 columns: k, observed kmers,
    and possible kmers.

    Parameters:
    sequence (str): the sequence you want to determine the values of k, and the observed
    and expected kmer values for

    Return:
    if seq is the incorrect value:
        str: gives a message about what is incorrect about the parameter
    if k and str are correct values:
        pd.DataFrame: a pandas dataframe having 3 columns: k, observed kmers, possible kmers

    """

###Checking to make sure sequence is correct value
    if sequence=='':
        return('Sequence must have a least a length of 1')
    if isinstance(sequence, int):
        return('Sequence must be a string')
    elif isinstance(sequence, float):
        return('Sequence must be a string')
    if sequence.isalpha()==False:
        return('Sequence must contain only letters.')
    guide=['A', 'a', 'C', 'c', 'T', 't', 'g', 'G']
    for letter in range(0, len(sequence)):
        if sequence[letter] not in guide:
            return('Sequence can only contain A, a, C, c, T, t, G, g as values.')

    possible = pd.DataFrame(columns=['k','observed kmers','possible kmers'])

    for i in range(1, len(sequence)+1):
        list2=[]
        for j in range(0, len(sequence)):
            if j+i > len(sequence):
                break
            word=sequence[j:j+i]
            if word not in list2:
                list2.append(word)
        obs_kmers=len(list2)

        if (len(sequence) > 4**i):
            temp=pd.DataFrame({'k':[i], 'observed kmers' :[obs_kmers], 'possible kmers' : [4**i]})
            possible=possible.append(temp, ignore_index=True)
        else:
            temp=pd.DataFrame({'k':[i], 'observed kmers' :[obs_kmers], 'possible kmers' : [len(sequence)-i+1]})
            possible=possible.append(temp, ignore_index=True)
    return(possible)



def linguistic_complexity(sequence):
    """
    Summary line: Generates the linguistic complexity for a sequence

    Extended Description: Includes a section that will test to make sure sequence
        is the correct value. If sequence is incorrect, will return a message about
        what is wrong with the parameter. The linguistic complexity is the proportion
        of kmers that are observed compared to the total number that are possibe.

    Parameters:
    sequence (str): the sequence you want to calculate the number linguistic complexity of

    Return:
    if sequence is the incorrect value:
        str: gives a message about what is incorrect about the parameter
    if sequence is correct value:
        float: gives the proportion of observed kmers to possible kmers (between 0 and 1)

    """
###Checking to make sure sequence is correct value
    if sequence=='':
        return('Sequence must have a least a length of 1')
    if isinstance(sequence, int):
        return('Sequence must be a string')
    elif isinstance(sequence, float):
        return('Sequence must be a string')
    if sequence.isalpha()==False:
        return('Sequence must contain only letters.')
    guide=['A', 'a', 'C', 'c', 'T', 't', 'g', 'G']
    for letter in range(0, len(sequence)):
        if sequence[letter] not in guide:
            return('Sequence can only contain A, a, C, c, T, t, G, g as values.')

    df=possible_kmers(sequence)
    ling_com=sum(df['observed kmers'])/sum(df['possible kmers'])
    return(ling_com)


def kmer_graph(sequence):
    """
    Summary line: Produces a graph of the proportion of each kmer observed

    Extended Description: Includes a section that will test to make sure sequence
        us the correct value.If sequence is correct, will return a message about
        what is wrong with the parameter.

    Parameters:
    sequence (str): the sequence you want to produce the graph of proportion of
    observed kmers

    Return:
    if sequence is the incorrect value:
        str: gives a message about what is incorrect about the parameter
    if sequence is correct value:
        plt: a bar plot of the size of k vs. the observed kmers for the sequence

    """

###Checking to make sure sequence is correct value
    if sequence=='':
        return('Sequence must have a least a length of 1')
    if isinstance(sequence, int):
        return('Sequence must be a string')
    elif isinstance(sequence, float):
        return('Sequence must be a string')
    if sequence.isalpha()==False:
        return('Sequence must contain only letters.')
    guide=['A', 'a', 'C', 'c', 'T', 't', 'g', 'G']
    for letter in range(0, len(sequence)):
        if sequence[letter] not in guide:
            return('Sequence can only contain A, a, C, c, T, t, G, g as values.')

    df=possible_kmers(sequence)
    freq=df['observed kmers']/df['possible kmers']
    x=np.arange(len(df['k']))
    plt.bar(x, freq)
    plt.xticks(x, df['k'])
    plt.ylabel('Frequency Observed')
    plt.xlabel('Size of kmers')
    plt.title('Frequency Observed vs. Size of kmer')
    #plt.show()
    return()

def main():
    #Input code based off of:
    #https://stackoverflow.com/questions/22939211/
       # what-is-the-proper-way-to-take-a-directory-path-as-user-input
    user_input = input("Enter the path of your file: ")

    assert os.path.exists(user_input), "I did not find the file at, "+str(user_input)
    f = open(user_input,'r+')
    sequence=f.read()
    sequence=list(sequence.split())
    f.close()
    for s in range(0, len(sequence)):
        print('For', sequence[s])
        k = input('Enter your kmer size:')
        k= int(k)

        kmer_count=count_kmers(k, sequence[s])
        if isinstance(kmer_count, str):
            print('Error: count_kmers function.', kmer_count)
            exit()

        dataframe=possible_kmers(sequence[s])
        if isinstance(dataframe, str):
            print('Error: possible_kmers function. ', dataframe)
            exit()

        lc=linguistic_complexity(sequence[s])
        if isinstance(lc, float):
            print('The linguistic complexity is for ', sequence[s], ' is ', lc)
            #only outputting this because that's what is sounded like in the prompt
        else:
            print('Error: linguistic_complexity function. ', lc)
            exit()

        g=kmer_graph(sequence[s])
        if isinstance(g, str):
            print('Error: kmer_graph function. ', g)
            exit()

if __name__ == '__main__':
    main()
