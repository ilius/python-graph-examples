#!/usr/bin/env python3

import random
from time import time

import igraph as ig


def overlaps(m1, d1, m2, d2):
	return abs(m1 - m2) < d1 + d2


def makeNormalRandomIntervals(n, m_mean=0, m_sigma=1.0, d_mean=0.1, d_sigma=0.01):
	lst = []
	for _i in range(n):
		m = random.normalvariate(m_mean, m_sigma)
		###
		d = 0
		while d <= 0:
			d = random.normalvariate(d_mean, d_sigma)
		###
		lst.append((m, d))
	return lst


def makeIntervalGraph(intervals):
	g = ig.Graph(1)
	n = len(intervals)
	g.add_vertices(n - 1)
	g.vs["name"] = list(range(n))
	g.vs["interval"] = intervals
	for i in range(1, n):
		mi, di = intervals[i]
		for j in range(i):
			mj, dj = intervals[j]
			if overlaps(mi, di, mj, dj):
				g.add_edges(
					[
						(i, j),
					],
				)
	return g


def assignComponents(g):
	g.vs["component"] = [None] * g.vcount()
	subgraphs = g.decompose()
	for compI, sg in enumerate(subgraphs):
		for v in sg.vs:
			origV = g.vs[v["name"]]
			origV["component"] = compI
	return len(subgraphs)


def hslToRgb(hue, sat, lum):
	# 0.0 <= h <= 360.0
	# 0.0 <= s <= 1.0
	# 0.0 <= l <= 1.0
	q = lum * (1.0 + sat) if lum < 0.5 else lum + sat - lum * sat
	p = 2 * lum - q
	hk = hue / 360.0
	tr = (hk + 1.0 / 3) % 1
	tg = hk % 1
	tb = (hk - 1.0 / 3) % 1
	rgb = []
	for tc in (tr, tg, tb):
		if tc < 1.0 / 6:
			c = p + (q - p) * 6 * tc
		elif 1.0 / 6 <= tc < 1.0 / 2:
			c = q
		elif 1.0 / 2 <= tc < 2.0 / 3:
			c = p + (q - p) * 6 * (2.0 / 3 - tc)
		else:
			c = p
		rgb.append(int(c * 255))
	return tuple(rgb)


def rgbToHtmlColor(r, g, b):
	return f"#{r:2x}{g:2x}{b:2x}"


def colorGraphByOrder(g, ordered_vertices, adjlist):
	for i in ordered_vertices:
		colors = set()
		for j in adjlist[i]:
			c = g.vs[j]["color"]
			if c is not None:
				colors.add(c)
		c = 0
		while c in colors:
			c += 1
		g.vs[i]["color"] = c


def colorGraph(g, alg):
	n = g.vcount()
	adjlist = g.get_adjlist()
	g.vs["color"] = [None] * n
	if alg == "random":
		ordered_vertices = list(range(n))
		random.shuffle(ordered_vertices)
	elif alg == "degree":
		ordered_vertices = [
			x[-1]
			for x in sorted(
				[(-g.degree(i), i) for i in range(n)],
			)
		]
	elif alg == "indegree":
		ordered_vertices = [
			x[-1]
			for x in sorted(
				[(-g.indegree(i), i) for i in range(n)],
			)
		]
	elif alg == "outdegree":
		ordered_vertices = [
			x[-1]
			for x in sorted(
				[(-g.outdegree(i), i) for i in range(n)],
			)
		]
	elif alg == "degree_len":
		ordered_vertices = [
			x[-1]
			for x in sorted(
				[(-g.degree(i), -g.vs[i]["interval"][1], i) for i in range(n)],
			)
		]
	elif alg == "len":
		ordered_vertices = [
			x[-1]
			for x in sorted(
				[(-g.vs[i]["interval"][1], i) for i in range(n)],
			)
		]
	elif alg == "-degree":
		ordered_vertices = [
			x[-1]
			for x in sorted(
				[(g.degree(i), i) for i in range(n)],
			)
		]
	elif alg == "sumdegree":
		ordered_vertices = [
			x[-1]
			for x in sorted(
				[
					(
						-sum([g.degree(j) for j in adjlist[i]]),
						i,
					)
					for i in range(n)
				],
			)
		]
	elif alg == "compact":
		for i in range(n):
			s = 0
			vi = g.vs[i]
			mi, di = vi["interval"]
			for j in range(n):
				mj, dj = g.vs[j]["interval"]
				s += float(abs(mi - mj)) / (di + dj)
			vi["compact"] = s
		ordered_vertices = [
			x[-1]
			for x in sorted(
				[(g.vs["compact"][i], i) for i in range(n)],
			)
		]
	elif alg == "degree_compact":
		for i in range(n):
			s = 0
			vi = g.vs[i]
			mi, di = vi["interval"]
			for j in range(n):
				mj, dj = g.vs[j]["interval"]
				s += float(abs(mi - mj)) / (di + dj)
			vi["compact"] = s
		ordered_vertices = [
			x[-1]
			for x in sorted(
				[(-g.degree(i), g.vs["compact"][i], i) for i in range(n)],
			)
		]
	else:
		raise ValueError(alg)
	##
	# print ordered_vertices
	colorGraphByOrder(g, ordered_vertices, adjlist)


def colorGraphSL(g):
	n = g.vcount()
	adjlist = g.get_adjlist()
	g.vs["color"] = [None] * n
	for row in sorted(
		[(-g.degree(i), i) for i in range(n)],
	):
		i = row[1]
		colors = set()
		for j in adjlist[i]:
			c = g.vs[j]["color"]
			if c is not None:
				colors.add(c)
		c = 0
		while c in colors:
			c += 1
		g.vs[i]["color"] = c


def testComponents():
	n = 12
	d_mean = 2.0 / n
	d_sigma = d_mean * 0.5
	###
	triesN = 1
	compN_sum = 0
	for _I in range(triesN):
		intervals = makeNormalRandomIntervals(
			n,
			d_mean=d_mean,
			d_sigma=d_sigma,
		)
		###
		graph = makeIntervalGraph(intervals)
		print(
			"\n".join(
				[
					"".join(
						[str(x) for x in row],
					)
					for row in graph.get_adjacency()
				],
			),
		)
		c_graphs = graph.decompose()
		compN_sum += len(c_graphs)
	print(f"average number of components: {float(compN_sum) / triesN:.1f}")


def drawComponents():
	n = 16
	d_mean = 2.0 / n
	d_sigma = d_mean * 0.5
	###
	triesN = 1
	for index in range(triesN):
		intervals = makeNormalRandomIntervals(
			n,
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
			i = v["component"]
			color = rgbToHtmlColor(
				*hslToRgb(
					360.0 * i / compN,
					1.0,
					0.5,
				),
			)
			colors.append(color)
		graph.write_svg(
			f"graph-{index}.svg",
			"fruchterman_reingold",
			# fruchterman_reingold
			# grid_fruchterman_reingold
			# graphopt
			# circle
			# sphere
			# reingold_tilford_circular
			# kamada_kawai
			# lgl
			colors=colors,
		)


def testColoring(showStat=1, doDraw=0):
	# degree is the best so far
	# sumdegree takes twice the time (than degree) for only 0.5% smaller colors count
	# degree take 10% or 20% more time than random, for 4.7% smaller colors count
	# -degree is even worst than random, by 7.9% bigger colors count
	algList = [
		"random",
		"degree",  # 3.0% less colors count
		# [color:degree + split] is faster than [split + color:degree]
		# 'indegree',## 7% slower than degree, the same colors count
		# 'outdegree',## 5% slower than degree, the same colors count
		# '-degree',
		# 'sumdegree',
		# 'compact', ## [color:compact] is slower than [split + color:compact]
		# 'degree_compact',
		# 'degree_len',## 9% slower than degree, the same colors count
		# 'len',## 0% to 2% faster than degree, 0.7% to 1.2% more colors count
	]
	do_split = 0
	repeat = 1000
	n = 20
	d_mean = 4.0 / n
	d_sigma = d_mean * 0.7
	###
	data = {
		alg: {
			"sumMaxColor": 0,
			"sumTime": 0,
		}
		for alg in algList
	}
	for stepI in range(repeat):
		intervals = makeNormalRandomIntervals(n, d_mean=d_mean, d_sigma=d_sigma)
		graph = makeIntervalGraph(intervals)
		###
		for _algI, alg in enumerate(algList):
			t0 = time()
			if do_split:
				maxStepColor = 0
				for g in graph.decompose():
					colorGraph(g, alg)
					maxStepColor = max([maxStepColor] + g.vs["color"])
				data[alg]["sumMaxColor"] += maxStepColor
			else:
				colorGraph(graph, alg)
				maxStepColor = max(graph.vs["color"])
				data[alg]["sumMaxColor"] += maxStepColor
				# graph.decompose()
			data[alg]["sumTime"] += time() - t0
			###
			if doDraw and graph.vcount() > 1:
				hcolors = [
                    rgbToHtmlColor(
                        *hslToRgb(
                            360.0 * v["color"] / (maxStepColor + 1),
                            1.0,
                            0.5,
                        ),
                    )
                    for v in graph.vs
                ]
				graph.write_svg(
					f"graph-{stepI}.svg",
					"circle",
					# fruchterman_reingold
					# grid_fruchterman_reingold
					# graphopt
					# circle
					# sphere
					# reingold_tilford_circular
					# kamada_kawai
					# lgl
					colors=hcolors,
					labels="color",
				)
	if showStat:
		alg0 = algList[0]
		for alg in algList:
			print(f"algorithm: {alg}")
			###
			timeStr = f"{data[alg]['sumTime'] / repeat:.2e}"
			if alg != alg0:
				percent = int((data[alg]["sumTime"] / data[alg0]["sumTime"] - 1) * 100)
				timeStr += " (" + ("+" if percent > 0 else "-") + "%d%%)" % abs(percent)
			###
			print(f"average coloring time: {timeStr}")
			color = float(data[alg]["sumMaxColor"]) / repeat + 1
			colorStr = f"{color:.2f}"
			if alg != alg0:
				color0 = float(data[alg0]["sumMaxColor"]) / repeat + 1
				# print color, color0, color/color0
				percent = abs((color / color0 - 1) * 100.0)
				colorStr += (
					" (" + ("+" if color >= color0 else "-") + f"{abs(percent):.1f}%)"
				)
			print(f"average colors count:  {colorStr}")
			print("---------------------------")


def testGraphCreation():
	repeat = 1000
	n = 20
	d_mean = 2.0 / n
	d_sigma = d_mean * 0.2
	###
	sumTime = 0
	for _stepI in range(repeat):
		intervals = makeNormalRandomIntervals(n, d_mean=d_mean, d_sigma=d_sigma)
		t0 = time()
		makeIntervalGraph(intervals)
		sumTime += time() - t0
	print(f"average time: {sumTime / repeat:e}")


def testColoringSL(showStat=1, doDraw=0):
	# the same 'degree' alg
	do_split = 1
	repeat = 1000
	n = 20
	d_mean = 2.0 / n
	d_sigma = d_mean * 0.2
	###
	sumMaxColor = 0
	sumTime = 0
	for stepI in range(repeat):
		intervals = makeNormalRandomIntervals(n, d_mean=d_mean, d_sigma=d_sigma)
		graph = makeIntervalGraph(intervals)
		###
		t0 = time()
		if do_split:
			maxStepColor = 0
			for g in graph.decompose():
				colorGraphSL(g)
				maxStepColor = max([maxStepColor] + g.vs["color"])
		else:
			colorGraphSL(graph)
			maxStepColor = max(graph.vs["color"])
			# graph.decompose()
		sumMaxColor += maxStepColor
		sumTime += time() - t0
		###
		if doDraw and graph.vcount() > 1:
			hcolors = [
                rgbToHtmlColor(
                    *hslToRgb(
                        360.0 * v["color"] / (maxStepColor + 1),
                        1.0,
                        0.5,
                    ),
                )
                for v in graph.vs
            ]
			graph.write_svg(
				f"graph-{stepI}.svg",
				"circle",
				# fruchterman_reingold
				# grid_fruchterman_reingold
				# graphopt
				# circle
				# sphere
				# reingold_tilford_circular
				# kamada_kawai
				# lgl
				colors=hcolors,
				labels="color",
			)
	if showStat:
		print(f"average coloring time: {sumTime / repeat:e}")
		print(f"average colors count:  {float(sumMaxColor) / repeat + 1}")


if __name__ == "__main__":
	testGraphCreation()
	testComponents()
	drawComponents()
	testColoring()
	testColoringSL()
