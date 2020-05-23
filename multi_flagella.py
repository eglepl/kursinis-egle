import numpy
import random
import matplotlib.pyplot as plt


t_end = 10
N = 1
M = 2
zone_a = 1000

M_ift = 3
alpha = 1 # how much rate drops on increse of length
gamma = 1 # length per element
zeta = 1 #elements to dissociate from flagella end.

beta = 0.01 # ift container, takes percent

lmb_1 = 100
lmb_2 = 100
lmb_3 = 100
lmb_4 = 100
lmb_p = 100
lmb_m = 100

mu = 0

s_p = 1
s_m = 1

class IFT:
    def __init__(self, x, d, c, z):
        self.x = x
        self.d = d
        self.c = c
        self.z = z

    def copy(self):
        return IFT(self.x, self.d, self.c, self.z)

    def __str__(self):
        return "{x:" + str(self.x) + ", d:"+ str(self.d) + ", c:" + str(self.c) + ", z:" + str(self.z)  + "}"

class S:
    def __init__(self):
        self.A = zone_a
        self.B = [1 for j in range(M)]
        self.U = [1 for j in range(M)]
        self.F = [IFT("E", -1, 0, 0) for i in range(N)]

    def copy(self):
        scp = S()
        scp.A = self.A
        scp.B = list(self.B)
        scp.U = list(self.U)
        scp.F = [f.copy() for f in self.F]

        return scp

    def __str__(self):
        return "A:" + str(self.A) + " B:" + ", ".join([str(b) for b in self.B]) + " U:" + ", ".join([str(u) for u in self.U]) + " F:" + " ,".join([str(f) for f in self.F])

s = S()

def L(s, j):
    return gamma * s.U[j]

def attach_to_flagella_zone_a(s, i, j):
    s.F[i].x = 0
    s.F[i].d = +1
    e = min(M_ift, int(s.A * beta))
    s.A += - e
    s.F[i].c = e
    s.F[i].z = j

def move_plus(s, i):
    s.F[i].x += s_p

def detach_to_zone_b(s, i):
    s.F[i].x = "E"
    s.F[i].d = +1
    s.U[s.F[i].z] += s.F[i].c
    s.F[i].c = 0

def attach_to_flagella_zone_b(s, i):
    s.F[i].x = L(s, s.F[i].z)
    s.F[i].d = -1
    e = min(M_ift, int(s.B[s.F[i].z] * beta))
    s.B[s.F[i].z] -= e
    s.F[i].c = e

def move_minus(s, i):
    s.F[i].x -= s_m

def detach_to_zone_a(s, i):
    s.F[i].x = "E"
    s.F[i].d = -1
    s.A += s.F[i].c
    s.F[i].c = 0

def flagella_end_dissociation(s, j):
    s.U[j] -= zeta
    if(s.U[j] < 1):
        s.U[j] = 1
    for f in s.F:
        if (f.x > L(s, j)):
            f.x = L(s, j)


def choice(p):
        r = random.random()
        psum = 0
        for i in range(0, len(p)):
            psum += p[i]
            if psum >= r:
                return i
        return len(p) - 1

t = 0

T = []
Y = []

# save initial state and time
T.append(t)
Y.append(s.copy())

while(t < t_end):
    V = list()
    for i in range(N):
        if (s.F[i].x == "E" and s.F[i].d == -1):
            for j in range(M):
                jj = j
                V.append([float(lmb_1) / (alpha * s.U[j]), lambda: attach_to_flagella_zone_a(s, i, jj)])
        elif (s.F[i].x == "E" and s.F[i].d == +1):
            V.append([lmb_3, lambda: attach_to_flagella_zone_b(s, i)])
        elif (s.F[i].x <= L(s, s.F[i].z) - s_p) and s.F[i].d == +1:
            V.append([lmb_p, lambda: move_plus(s, i)])
        elif (s.F[i].x > L(s, s.F[i].z) - s_p) and s.F[i].d == +1:
            V.append([lmb_2, lambda: detach_to_zone_b(s, i)])
        elif (s.F[i].x >= s_m) and s.F[i].d == -1:
            V.append([lmb_m, lambda: move_minus(s, i)])
        elif (s.F[i].x < s_m) and s.F[i].d == -1:
            V.append([lmb_4, lambda: detach_to_zone_a(s, i)])
        else:
            print("something wrong!!!")
    for j in range(M):
        V.append([mu, lambda: flagella_end_dissociation(s, j)])

    a0 = 0
    for v in V:
        a0 += v[0]

    print(V)
    print(s)

    dt = numpy.random.exponential(1.0/a0)
    t += dt
    if (t <= t_end):
        p = [v[0] / a0 for v in V]
        idx = choice(p)
        fun = V[idx][1]
        fun()

    T.append(t)
    Y.append(s.copy())

#print([y.F[0].x for y in Y])

# Plot the simulation data.
plt.figure()
# Create a 5% (0.05) and 10% (0.1) padding in the
# x and y directions respectively.
plt.margins(0.05, 0.1)
plt.grid(True)
#plt.plot(T, [y.A for y in Y], label="A")
plt.plot(T, [y.U[0] for y in Y], label="U1")
plt.plot(T, [y.U[1] for y in Y], label="U2")
#plt.plot(T, [y.U[2] for y in Y], label="U3")
plt.legend()
plt.xlabel('Time')
plt.ylabel('Elements')
plt.show()