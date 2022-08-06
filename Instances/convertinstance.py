import sys


with open("Instances/" + sys.argv[1] + "/" + sys.argv[2]) as file:
    lines = file.readlines()

current = lines[0].rstrip().split()

J = int(current[0])
I = int(current[1])

s = []
f = []

for j in range(1, J+1):
    current = lines[j].rstrip().split()
    s.append(float(current[0]))
    f.append(float(current[1]))


d = [float(x) for x in lines[J+1].rstrip().split()]

c = []

for j in range(0, J):
    c.append([float(x) for x in lines[J+2 + j].rstrip().split()])

# print(J)
# print(I)
# print(s)
# print(f)
# print(c)

estimations = [[0 for j1 in range(0, J)] for j2 in range(0, J)]

for j1 in range(0, J):
    for j2 in range(0, j):
        minimum = c[j1][0] + c[j2][0]
        for i in range(1, I):
            if (c[j1][i] + c[j2][i] < minimum):
                minimum = c[j1][i] + c[j2][i]
        estimations[j1][j2] = minimum
        estimations[j2][j1] = minimum
    estimations[j1][j1] = 2 * max(c[j1])

# print(estimations)

nbvoisins = int(sys.argv[3])
K = int(sys.argv[4])

V = []

for j in range(0, J):
    V.append(sorted(range(J), key = lambda site: estimations[j][site])[:nbvoisins])

# print(J)
# print(V)


instance = open("Instances/Reliable/" + str(nbvoisins) + "_" + str(K) + "_"  + sys.argv[1] + "_" + sys.argv[2] + ".txt", "w")

instance.write("J = " + str(J) + " ;" + '\n')
instance.write("I = " + str(I) + " ;" + '\n')
instance.write("K = " + str(K) + " ;" + '\n')
instance.write("c = " + str(s) + " ;" + '\n')
instance.write("d = " + str(d) + " ;" + '\n')
instance.write("b = " + str(f) + " ;" + '\n')
instance.write("a = " + str(c) + " ;" + '\n')
instance.write("V = " + str(V) + " ;" + '\n') 