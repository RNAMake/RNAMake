
f = open("test_gu.dat")
lines = f.readlines()
f.close()

points = []
minp, maxp = 1000, -1000

for l in lines:
    d = float(l)
    points.append(d)
    if d > maxp:
        maxp = d
    if d < minp:
        minp = d

bins = [0 for x in range(0,101)]
incr = (maxp - minp) / 100.0
current = minp
i = 0
for p in points:
    rounded = round((p - minp) / incr)
    bins[int(rounded)] += 1

for i in range(0,101):
    print minp+i*incr,bins[i]



print maxp, minp
