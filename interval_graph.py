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
    #g.add_edges([(i, i) for i in range(n)])
    for i in range(1, n):
        mi, di = intervals[i]
        for j in range(i):
            mj, dj = intervals[j]
            if overlaps(mi, di, mj, dj):
                g.add_edges([
                    (i, j),
                    (j, i),
                ])
    return g

def assignComponents(g):
    g.vs['component'] = [None]*g.vcount()
    subgraphs = g.decompose()
    for compI, sg in enumerate(subgraphs):
        for v in sg.vs:
            origV = g.vs[v['name']]
            origV['component'] = compI
    return len(subgraphs)


def hslToRgb(h, s, l):
    ## 0.0 <= h <= 360.0
    ## 0.0 <= s <= 1.0
    ## 0.0 <= l <= 1.0
    if l < 0.5:
        q = l * (1.0+s)
    else:
        q = l + s - l*s
    p = 2*l - q
    hk = h/360.0
    tr = (hk+1.0/3) % 1
    tg = hk % 1
    tb = (hk-1.0/3) % 1
    rgb = []
    for tc in (tr, tg, tb):
        if tc < 1.0/6:
            c = p + (q-p)*6*tc
        elif 1.0/6 <= tc < 1.0/2:
            c = q
        elif 1.0/2 <= tc < 2.0/3:
            c = p + (q-p)*6*(2.0/3-tc)
        else:
            c = p
        rgb.append(int(c*255))
    return tuple(rgb)

rgbToHtmlColor = lambda r, g, b: '#%.2x%.2x%.2x'%(r, g, b)

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


def testComponents():
    n = 50
    m_sigma = 100.0
    d_mean = m_sigma * 1.0 / n
    d_sigma = d_mean * 0.5
    ###
    triesN = 5
    compN_sum = 0
    for I in range(triesN):
        intervals = makeNormalRandomIntervals(
            n,
            m_sigma=m_sigma,
            d_mean=d_mean,
            d_sigma=d_sigma,
        )
        ###
        graph = makeIntervalGraph(intervals)
        c_graphs = graph.decompose()
        compN_sum += len(c_graphs)
    print 'average number of components: %.1f'%(float(compN_sum)/triesN)

def drawComponents():
    n = 16
    m_sigma = 100.0
    d_mean = m_sigma * 2.0 / n
    d_sigma = d_mean * 0.5
    ###
    triesN = 1
    compN_sum = 0
    for I in range(triesN):
        intervals = makeNormalRandomIntervals(
            n,
            m_sigma=m_sigma,
            d_mean=d_mean,
            d_sigma=d_sigma,
        )
        ###
        graph = makeIntervalGraph(intervals)
        ###
        compN = assignComponents(graph)
        ###
        colors = []
        for v in graph.vs:
            i = v['component']
            color = rgbToHtmlColor(*hslToRgb(
                360.0 * i / compN,
                1.0,
                0.5,
            ))
            colors.append(color)
        graph.write_svg(
            'graph-%s.svg'%I,
            'fruchterman_reingold',
            ## fruchterman_reingold
            ## grid_fruchterman_reingold
            ## graphopt
            ## circle
            ## sphere
            ## reingold_tilford_circular
            ## kamada_kawai
            ## lgl
            colors=colors,
        )


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
                for g in graph.decompose():
                    colorGraph(g, alg)
                    maxStepColor = max([maxStepColor] + g.vs['color'])
                data[alg]['sumMaxColor'] += maxStepColor
            else:
                colorGraph(graph, alg)
                data[alg]['sumMaxColor'] += max(graph.vs['color'])
                #graph.decompose()
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
    #testComponents()
    drawComponents()
    #testColoring()

