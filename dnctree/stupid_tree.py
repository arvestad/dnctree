#! /usr/bin/env python3
import argparse

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


def main():
    ap = argparse.ArgumentParser(description='Create non-random test trees for dnctree')
    ap.add_argument('-b', '--branchlength', type=float, metavar='b', default=0.1, help='Branchlength to use on all edges.')
    ap.add_argument('-d', '--depth', type=int, help='Bifurcating tree. Specify #levels in the tree')
    ap.add_argument('-c', '--caterpillar', type=int, metavar='n', help='Caterpillar-shaped tree. Specify #leaves.')
    args = ap.parse_args()

    if args.depth:
        assert args.depth > 0
        t = even_steven(args.depth, branchlength=args.branchlength)
    elif args.caterpillar:
        assert args.caterpillar > 2
        t = caterpillar(args.caterpillar, branchlength=args.branchlength)
    else:
        raise('No. You have to say what kind of tree to generate.')
    print(t)


if __name__ == '__main__':
    main()
