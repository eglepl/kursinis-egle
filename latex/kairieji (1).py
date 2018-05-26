#!/usr/bin/env python


def f(x):
  return x*x


def integralas(a,b,N,funkcija):
  h = (b-a)*1.0/N
  x_taskai = [i*h for i in range(N)]
  y_taskai = [funkcija(x) for x in x_taskai]
  return sum([ y_taskai[ i ] * h for i in range(N) ])

a = 0.0
b = 1.0
N1 = 100
N2 = 1000


print "N1=" + str(N1) + " integralas = " + str(integralas(a,b,N1, f))
print "N2=" + str(N2) + " integralas = " + str(integralas(a,b,N2, f))
print " integralas = " + str(1.0/3)

