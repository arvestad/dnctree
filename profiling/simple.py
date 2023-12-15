import cProfile
import pstats
from pstats import SortKey

from alv.io import read_alignment
import dnctree

def dnc_command():
    #    seqs, _ = read_alignment('six.phy', 'guess', 'phylip', None, None)
    seqs, _ = read_alignment('../testdata/s83_L500.phylip', 'guess', 'phylip', None, None)

    msa = dnctree.MSA(seqs)
    t, _ = dnctree.divide_n_conquer_tree(msa, base_case_size=20)
    #    print(str(t))

if __name__ == "__main__":
    profiling_output = 'profiling_stats.txt'
    cProfile.run('dnc_command()', profiling_output)

    p = pstats.Stats(profiling_output)
    p.strip_dirs().sort_stats(SortKey.TIME).print_stats(10)    

