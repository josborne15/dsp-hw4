import pytest
import sys

from JO_HW4 import *


def test_count_kmers():
    assert count_kmers(-1, 'ATTTGGATT' ) == 'Not a valid k. k > 0'
    assert count_kmers(1,'') == 'Sequence must have a least a length of 1'
    assert count_kmers(2.5, 'ATTTGGATT') == 'k must be a whole number'
    assert count_kmers('T', 'THH') == 'k must be a integer'
    assert count_kmers(3, 5654) == 'Sequence must be a string'
    assert count_kmers(3, 2.5) == 'Sequence must be a string'
    assert count_kmers(3, '-') == 'Sequence must contain only letters.'
    assert count_kmers(10, 'ATTTGGATT') == 'k cannot be longer than the length of the sequence'
    assert count_kmers(1, 'GATCF') == 'Sequence can only contain A, a, C, c, T, t, G, g as values.'

def test_possible_kmers():
    assert possible_kmers('') == 'Sequence must have a least a length of 1'
    assert possible_kmers(5654) == 'Sequence must be a string'
    assert possible_kmers(2.5) == 'Sequence must be a string'
    assert possible_kmers('-') == 'Sequence must contain only letters.'
    assert possible_kmers('GATCF') == 'Sequence can only contain A, a, C, c, T, t, G, g as values.'

def test_linguistic_complexity():
    assert linguistic_complexity('') == 'Sequence must have a least a length of 1'
    assert linguistic_complexity(5654) == 'Sequence must be a string'
    assert linguistic_complexity(2.5) == 'Sequence must be a string'
    assert linguistic_complexity('-') == 'Sequence must contain only letters.'
    assert linguistic_complexity('GATCF') == 'Sequence can only contain A, a, C, c, T, t, G, g as values.'

def test_kmer_graph():
    assert kmer_graph('') == 'Sequence must have a least a length of 1'
    assert kmer_graph(5654) == 'Sequence must be a string'
    assert kmer_graph(2.5) == 'Sequence must be a string'
    assert kmer_graph('-') == 'Sequence must contain only letters.'
    assert kmer_graph('GATCF') == 'Sequence can only contain A, a, C, c, T, t, G, g as values.'
