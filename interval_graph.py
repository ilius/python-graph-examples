import igraph as ig
import random
from time import time

overlaps = lambda m1, d1, m2, d2: abs(m1-m2) < d1+d2


def makeNormalRandomIntervals(n, m_mean=0, m_sigma=50, d_mean=20, d_sigma=10):
    l = []
    for i in range(n):
        m = random.normalvariate(m_mean, m_sigma)
        ###
        d = 0
        while d <= 0:
            d = random.normalvariate(d_mean, d_sigma)
        ###
        l.append((m, d))
    return l

def makeIntervalGraph(intervals):
    g = ig.Graph()
    n = len(intervals)
    g.add_vertices(n-1)
    g.vs['name'] = range(n)
    g.vs['interval'] = intervals
    for i in range(1, n):
        mi, di = intervals[i]
        for j in range(i-1):
            mj, dj = intervals[j]
            if overlaps(mi, di, mj, dj):
                g.add_edges([
                    (i, j),
                    (j, i),
                ])
    return g

def splitGraphIntoComponents(g):
    cc = g.biconnected_components()
    return [cc.subgraph(i) for i in range(len(cc))]

def colorGraphByOrder(g, ordered_vertices):
    adjlist = g.get_adjlist()
    for i in ordered_vertices:
        colors = set()
        for j in adjlist[i]:
            c = g.vs[j]['color']
            if c is not None:
                colors.add(c)
        c = 0
        while c in colors:
            c += 1
        g.vs[i]['color'] = c

def colorGraph(g, alg):
    n = g.vcount()
    g.vs['color'] = [None]*n
    if alg=='compact':
        for i in range(n):
            s = 0
            vi = g.vs[i]
            mi, di = vi['interval']
            for j in range(n):
                mj, dj = g.vs[j]['interval']
                s += float(abs(mi-mj))/(di+dj)
            vi['compact'] = s
        ordered_vertices = [
            x[1] for x in sorted(
                [
                    (g.vs['compact'][i], i) for i in range(n)
                ],
            )
        ]
    elif alg=='degree_compact':
        for i in range(n):
            s = 0
            vi = g.vs[i]
            mi, di = vi['interval']
            for j in range(n):
                mj, dj = g.vs[j]['interval']
                s += float(abs(mi-mj))/(di+dj)
            vi['compact'] = s
        ordered_vertices = [
            x[2] for x in sorted(
                [
                    (-g.indegree(i), g.vs['compact'][i], i) for i in range(n)
                ],
            )
        ]
    elif alg=='degree':
        ordered_vertices = [
            x[1] for x in sorted(
                [
                    (-g.indegree(i), i) for i in range(n)
                ],
            )
        ]
    elif alg=='random':
        ordered_vertices = range(n)
        random.shuffle(ordered_vertices)
    else:
        raise ValueError(alg)
    ##
    #print ordered_vertices
    colorGraphByOrder(g, ordered_vertices)

def testComponentCount():
    n = 20
    m_sigma = 100.0
    d_mean = m_sigma * 8.0 / n
    d_sigma = d_mean * 0.5
    ###
    triesN = 20
    compN_sum = 0
    for I in range(triesN):
        intervals = makeNormalRandomIntervals(
            n,
            m_sigma=m_sigma,
            d_mean=d_mean,
            d_sigma=d_sigma,
        )
        ###
        graphs = splitGraphIntoComponents(
            makeIntervalGraph(intervals)
        )
        compN_sum += len(graphs)
    print 'average number of components: %.1f'%(float(compN_sum)/triesN)

def testColoring():
    ## degree is faster than compact and degree_compact, much faster than random
    #alg = 'random'
    alg = 'degree' ## [color:degree + split] is faster than [split + color:degree]
    #alg = 'compact' ## [color:compact] is slower than [split + color:compact]
    #alg = 'degree_compact'
    do_split = 0
    repeat = 10000
    n = 40
    m_sigma = 100.0
    d_mean = m_sigma * 3.0 / n
    d_sigma = d_mean * 0.5
    ###
    t0 = time()
    sumMaxColor = 0.0
    for i in range(repeat):
        intervals = makeNormalRandomIntervals(n, m_sigma=m_sigma, d_mean=d_mean, d_sigma=d_sigma)
        graph = makeIntervalGraph(intervals)
        ###
        if do_split:
            maxStepColor = 0
            for graph in splitGraphIntoComponents(graph):
                colorGraph(graph, alg)
                maxStepColor = max([maxStepColor] + graph.vs['color'])
            sumMaxColor += maxStepColor
        else:
            colorGraph(graph, alg)
            sumMaxColor += max(graph.vs['color'])
            splitGraphIntoComponents(graph)
    dt = time() - t0
    print 'algorithm: %s'%alg
    print '  total coloring time: %.1f'%dt
    print 'average coloring time: %e'%(dt/repeat)
    print 'average colors count: %s'%(float(sumMaxColor)/repeat + 1)




if __name__=='__main__':
    #testComponentCount()
    testColoring()

