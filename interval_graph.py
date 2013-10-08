import igraph as ig
import random

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

def colorGraph(g):
    n = g.vcount()
    al = g.get_adjlist()
    g.vs['color'] = [None]*n
    for i in range(n):
        s = 0
        vi = g.vs[i]
        mi, di = vi['interval']
        for j in range(n):
            mj, dj = g.vs[j]['interval']
            s += float(abs(mi-mj))/(di+dj)
        vi['sum'] = s
    ordered_vertices = [
        x[1] for x in sorted(
            [
                (s, i) for i, s in enumerate(g.vs['sum'])
            ],
            ##reverse=True,
        )
    ]
    ##random.shuffle(ordered_vertices)
    #print ordered_vertices
    for i in ordered_vertices:
        colors = set()
        for j in al[i]:
            c = g.vs[j]['color']
            if c is not None:
                colors.add(c)
        c = 0
        while c in colors:
            c += 1
        g.vs[i]['color'] = c

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
    n = 200
    m_sigma = 100.0
    d_mean = m_sigma * 3.0 / n
    d_sigma = d_mean * 0.5
    ###
    intervals = makeNormalRandomIntervals(n, m_sigma=m_sigma, d_mean=d_mean, d_sigma=d_sigma)
    ###
    graphs = splitGraphIntoComponents(
        makeIntervalGraph(intervals)
    )
    print '%s vertices'%n
    print '%s edges'%sum([g.ecount() for g in graphs])
    print '%s components'%len(graphs)
    for i, graph in enumerate(graphs):
        colorGraph(graph)
    colors_count = max([
        max(g.vs['color']) for g in graphs
    ]) + 1
    print '%s colors'%colors_count




if __name__=='__main__':
    #testComponentCount()
    testColoring()

