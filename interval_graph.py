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

def colorGraphByOrder(g, ordered_vertices, adjlist):
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
    adjlist = g.get_adjlist()
    g.vs['color'] = [None]*n
    ordered_vertices = range(n)
    random.shuffle(ordered_vertices)
    if alg=='random':
        pass
    elif alg=='degree':
        ordered_vertices = [
            x[1] for x in sorted(
                [
                    (-g.indegree(i), i) for i in range(n)
                ],
            )
        ]
    elif alg=='-degree':
        ordered_vertices = [
            x[1] for x in sorted(
                [
                    (g.indegree(i), i) for i in range(n)
                ],
            )
        ]
    elif alg=='sumdegree':
        ordered_vertices = [
            x[1] for x in sorted(
                [
                    (
                        -sum([
                            g.indegree(j) for j in adjlist[i]
                        ]),
                        i,
                    )
                    for i in range(n)
                ]
            )
        ]
    elif alg=='compact':
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
    else:
        raise ValueError(alg)
    ##
    #print ordered_vertices
    colorGraphByOrder(g, ordered_vertices, adjlist)

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
    ## degree is the best so far
    ## sumdegree takes twice the time (than degree) for only 0.5% smaller colors count
    ## degree take 10% or 20% more time than random, for 4.7% smaller colors count
    ## -degree is even worst than random, by 7.9% bigger colors count
    algList = [
        'random',
        'degree', ## [color:degree + split] is faster than [split + color:degree]
        #'-degree',
        #'sumdegree',
        #'compact', ## [color:compact] is slower than [split + color:compact]
        #'degree_compact',
    ]
    do_split = 0
    repeat = 5000
    n = 50
    m_sigma = 100.0
    d_mean = m_sigma * 4.0 / n
    d_sigma = d_mean * 0.7
    ###
    data = dict([
        (alg, {
            'sumMaxColor': 0,
        }) for alg in algList
    ])
    for i in range(repeat):
        intervals = makeNormalRandomIntervals(n, m_sigma=m_sigma, d_mean=d_mean, d_sigma=d_sigma)
        graph = makeIntervalGraph(intervals)
        ###
        for i, alg in enumerate(algList):
            data[alg]['time'] = time()
            if do_split:
                maxStepColor = 0
                for g in splitGraphIntoComponents(graph):
                    colorGraph(g, alg)
                    maxStepColor = max([maxStepColor] + g.vs['color'])
                data[alg]['sumMaxColor'] += maxStepColor
            else:
                colorGraph(graph, alg)
                data[alg]['sumMaxColor'] += max(graph.vs['color'])
                #splitGraphIntoComponents(graph)
            data[alg]['time'] = time() - data[alg]['time']
    alg0 = algList[0]
    for alg in algList:
        print 'algorithm: %s'%alg
        ###
        timeStr = '%.2e'%(data[alg]['time']/repeat)
        if alg != alg0:
            percent = int((data[alg]['time']/data[alg0]['time'] - 1) * 100)
            timeStr += ' (' + ('+' if percent > 0 else '-')  + '%d%%)'%abs(percent)
        ###
        print 'average coloring time: %s'%timeStr
        color = float(data[alg]['sumMaxColor'])/repeat + 1
        colorStr = '%.2f'%color
        if alg != alg0:
            color0 = float(data[alg0]['sumMaxColor'])/repeat + 1
            #print color, color0, color/color0
            percent = abs((color/color0 - 1) * 100.0)
            colorStr += ' (' + ('+' if color >= color0 else '-')  + '%.1f%%)'%abs(percent)
        print 'average colors count:  %s'%colorStr
        print '---------------------------'




if __name__=='__main__':
    #testComponentCount()
    testColoring()

