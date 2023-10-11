#! /usr/bin/env python3
import argparse

lba_length = 1.0

def even_steven(depth, branchlength):
    t, _ = even_steven_helper(depth, branchlength, idx=0, n_kids=3)
    return t + ';'


def even_steven_helper(depth, branchlength, idx, n_kids):
    assert(depth > 0)

    if depth == 1:
        return create_node(idx, branchlength), idx+1
    else:
        subtrees = list()
        for i in range(n_kids):
            t, idx = even_steven_helper(depth-1, branchlength, idx, 2)
            subtrees.append(t)

        subtree_string = ','.join(subtrees)
        tree_string = f'({subtree_string}):{branchlength}'
        return tree_string, idx

def caterpillar(n, branchlength):
    t1=create_node(0, branchlength)
    t2=create_node(1, branchlength)
    tn = caterpillar_helper(n-2, 2, branchlength)

    return f'({t1},{t2},{tn});'

def caterpillar_helper(n, idx, branchlength):
    if n == 1:
        return create_node(idx, branchlength)
    else:
        leaf = create_node(idx, branchlength)
        t = caterpillar_helper(n-1, idx+1, branchlength)
        return f'({leaf},{t}):{branchlength}'

def create_node(idx, branchlength):
    return f'L{idx}:{branchlength}'


def lba_model_tree(depth, long_branchlength, regular_branchlength):
    '''
    Create a tree that could be susceptible to long branch attraction (LBA).
    '''
    idx = 0
    t1, idx = lba_model_helper(depth-1, long_branchlength, regular_branchlength, idx)
    t2, idx = even_steven_helper(depth-1, regular_branchlength, idx, 2)
    t3, idx = lba_model_helper(depth-1, long_branchlength, regular_branchlength, idx)
    t4, idx = even_steven_helper(depth-1, regular_branchlength, idx, 2)

    L = long_branchlength
    t = f'(({t1},{t2}):{regular_branchlength},{t3}, {t4});'
    return t

def lba_model_helper(depth, long_branchlength, regular_branchlength, idx):
    assert depth > 1
    t1, idx = even_steven_helper(depth-1, regular_branchlength, idx, 2)
    t2, idx = even_steven_helper(depth-1, regular_branchlength, idx, 2)
    return f'({t1},{t2}):{long_branchlength}', idx


def main():
    ap = argparse.ArgumentParser(description='Create non-random test trees for dnctree')
    ap.add_argument('-b', '--branchlength', type=float, metavar='b', default=0.1, help='Branchlength to use on all edges.')
    ap.add_argument('-d', '--depth', type=int, help='Bifurcating tree. Specify #levels in the tree')
    ap.add_argument('-c', '--caterpillar', type=int, metavar='n', help='Caterpillar-shaped tree. Specify #leaves.')
    ap.add_argument('-l', '--lba-tree', type=int, metavar='n', help='Tree inducing(?) long branch attraction. Specify depth of the four subtrees.')
    ap.add_argument('-L', type=float, metavar='b', default=lba_length, help=f'Length of a long branch in the LBA test case. Default is {lba_length}.')
    args = ap.parse_args()

    if args.depth:
        assert args.depth > 0
        t = even_steven(args.depth, branchlength=args.branchlength)
    elif args.caterpillar:
        assert args.caterpillar > 2
        t = caterpillar(args.caterpillar, branchlength=args.branchlength)
    elif args.lba_tree:
        t = lba_model_tree(args.lba_tree, args.L, args.branchlength)
    else:
        ap.exit('No. You have to say what kind of tree to generate.')
    print(t)


if __name__ == '__main__':
    main()
