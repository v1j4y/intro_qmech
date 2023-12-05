#import "template.typ": project
#import "@preview/physica:0.8.1": *

#show: project.with(
  title: "Introduction to Model hamiltonians",
  authors: (
    (name: "Vijay Gopal CHILKURI", 
    email: "vijay-gopal.chilkuri@uni-amu.fr", 
    affiliation: "University of Aix-Marseille", 
    postal: "Avenue Escadrille Niemen - 130013", 
    phone: "+33413945595"),
  ),
  abstract: [A simple introduction to model
  hamiltonians for quantum chemists.]
)

#set math.equation(numbering: "1.")


= Introduction

Model hamiltonians are of premordial importance 
for understanding chemical and physical behavior
of molecules and materials.
Here, we shall briefly describe the various models
and their formulation in as simple terms as possible.

= Derivation of the Schrodinger equation

The schrodinger equation can be derived using the
path integral formulation as shown by Feynman.@feynman1948

== Lagrangian

In order to demonstrate the derivation by Feynman,
one needs to first define the notion of the
lagrangian @eq:deriv1.

$
L = T - V
$ <eq:deriv1>

Where, $T$ is the kinetic energy and $V$ is the 
potential energy. In order to better understand 
the lagrangian and its relation to Newton's 
equations of motion, in @eq:deriv2, @eq:deriv3 we derive the equations
of motion in lagrange formulation and its connection
to the usual newtons equations of motion.

$
L(r, dot(r)) &= T - V \
L(r, dot(r)) &= sum_i frac(1,2)m dot(r)_i^2 - V(r_1,...,r_n) \
$ <eq:deriv2>

where, $r$ is the position and $dot(r)=frac(d r,d t)$ 
is the velocity. 
Using this definition of the lagrangian, we
can derive the so called Euler-Lagrange equation
which is equivalent to Newton's equation @eq:deriv3.

$
frac(d, d t)(pdv(L,dot(r)_i))  - pdv(L,r_i)&= 0
$ <eq:deriv3>

Where, the second term on the left of @eq:deriv3 is the derivative
of the potential i.e. the force (@eq:deriv4).

$
pdv(L,r_i) = pdv(V(r_1,...,r_n),r_i) = F_i
$ <eq:deriv4>

and the first term of @eq:deriv3 is the acceleration (@eq:deriv5).

$
pdv(,t) pdv(L,dot(r)_i) = m pdv(dot(r)_i,t) = m a_i
$ <eq:deriv5>

where, $a_i = dot.double(r)_i$. Therefore @eq:deriv3 
is equivalent to Newton's equation (@eq:deriv6).

$
F_i = m a_i
$ <eq:deriv6>

== Action

The action is defined as the integral of
the Lagrangian along a specific path between
two points, $A$ at time $t_a$ to point $B$ in time $t_b$ @eq:action1.

$
S[r(t)] = integral_(t_a)^(t_b) L(r(t),dot(r)(t)) d t
$ <eq:action1>

The action is an important quantity and describes 
the weight and phase of each path. Using the 
action, we can derive the @eq:deriv3. This can 
be done using the principle of least action 
which says that the path that survives is the one
that minimises the action @eq:action2.

$
integral_(t_a)^(t_b) delta L d t &= 0 \
  delta S &= 0
$ <eq:action2>

The derivation of @eq:deriv3 follows from the 
above @eq:action2 once it is simplified using
integration by parts (@eq:action3).

$

integral_(t_a)^(t_b) delta L d t &=
integral_(t_a)^(t_b) sum_i^n ( pdv(L,r_i) + frac(d ,d t)(pdv(L,dot(r)_i)) -  frac(d, d t) pdv(L,dot(r)_i) ) d t\
 &= sum_i^n [ pdv(L,dot(q))delta q_j ]_(t_a)^(t_b)
 + integral_(t_a)^(t_b) sum_i^n ( pdv(L,r_i) - frac(d, d t)(pdv(L,dot(r)_i))  ) d t
$ <eq:action3>

Therefore, if $integral_(t_a)^(t_b) delta L d t = 0$,
the left hand side of @eq:action3 is $0$. All terms of
@eq:action3 including the value of the integral. This
implies @eq:deriv3.

== Postulates of Feynman

Feynman put forth two postulates to derive
the schrodinger equation.
The first postulate says that the total
action is the sum of the actions of individual
paths, i.e. @eq:postul1.

$
S = sum_i S[r(t)_i]
$ <eq:postul1>

The second postulate says that the wavefunction
$phi$ can be be expressed as an exponential
function of the position $r(t)$ and its first
deriavtive $dot(r)(t)$, i.e. @eq:postul2.

$
phi(x_k,t) = lim_(epsilon arrow.r 0) integral_R 
exp( frac(i,planck.reduce) sum_i S[r(t)_i]) ...frac(d x_(i-1), A) frac(d x_(i-2), A)...
$ <eq:postul2>

where the integral is over the region $R$ which
contains all the paths.

== Derivation

The equation of motion describes the evolution
of the wavefunction $phi(x_(k+1),t)$ from
time $t$ to time $t+epsilon$ (@eq:deriv7).

$
phi(x_(k+1), t + epsilon) = lim_(epsilon arrow.r 0) integral_R 
exp( frac(i,planck.reduce) sum_i S[r(t)_i]) ...frac(d x_(i), A) frac(d x_(i-1), A)...
$ <eq:deriv7>

Using the definition of $phi(x_k, t)$ given in @eq:postul2,
we can use it to obtain the wavefunction at time
$t+epsilon$ (@eq:deriv8).

$
phi(x_(k+1),t+epsilon) =  [ integral_R S[x_(k+1),x_k]] phi(x_k, t)frac(d x_k, A)
$ <eq:deriv8>

The integral in @eq:deriv8 can be interpreted as 
the hamiltonian once we substitute the action (@eq:deriv9).

$
S(x_(k+1),x_k) = frac(m epsilon,2) (frac(x_(k+1)-x_k,epsilon))^2 - epsilon V(x_(k+1))
$ <eq:deriv9>

now the @eq:deriv8 becomes,

$

phi(x_(k+1),t+epsilon) =  [ integral frac(m epsilon,2) (frac(x_(k+1)-x_k,epsilon))^2 - epsilon V(x_(k+1))] phi(x_k, t)frac(d x_k, A)
$ <eq:deriv10>

Expanding the wavefunction $phi(x_(k+1),t)$ around $x_k$ using the
taylor series gives,

$
phi(x_(k+1),t+epsilon) =\  exp(frac(-i epsilon V, planck.reduce)) times integral exp(frac(i epsilon xi^2,2 planck.reduce epsilon))[psi(x, t) - xi pdv(psi(x,t),x)  + frac(xi^2,2) pdv(psi(x,t),x,2) - ...] (d xi)/A
$ <eq:deriv11>

where, $x_(k+1) - x_k = xi$. Expanding the left hand
also around $xi$ gives.

$
phi(x_(k+1),t) + epsilon pdv(phi(x,t),t)= \ exp(frac(-i epsilon V, planck.reduce)) times integral exp(frac(i epsilon xi^2,2 planck.reduce epsilon))[psi(x, t) - xi pdv(psi(x,t),x)  + frac(xi^2,2) pdv(psi(x,t),x,2) - ...] (d xi)/A
$ <eq:deriv12>

The factors in the integrand on the right of
@eq:deriv12 which contain $xi$, $xi^3$ etc are
zero because they are odd integrals (@eq:deriv13).

$

phi(x_(k+1),t) + epsilon pdv(phi(x,t),t)= \ exp(frac(-i epsilon V, planck.reduce)) times 
sqrt(2 pi planck.reduce i /m)/A [  psi(x, t) + frac(planck.reduce epsilon i,2 m) pdv(psi(x,t),x,2) + ...] 
$ <eq:deriv13>

Finally, equating the terms of same order in $epsilon$,
we get @eq:deriv14

$
-planck.reduce/i pdv(psi,t) &= 1/(2 m) (planck.reduce/i pdv(,x) )^2 psi + V(x) psi \
-planck.reduce/i pdv(psi,t) &= H psi
$ <eq:deriv14>

The above equation can be compared to the 
time dependent schrodinger equation. The
time independent form describes stationary wavefunctions
which is given as @eq:deriv15.

$
H psi = lambda psi
$ <eq:deriv15>

= 1D particle in a box

The problem of a particle in a box can be defined
as shown in @eq:part1

$
frac(d^2 psi, d x^2) = lambda psi
$ <eq:part1>

= Second quantization

= Huckel hamiltonian

= Hubbard hamiltonian

= Double exchange hamiltonian

// Bibliography section
#pagebreak(weak: true)
#set page(header: [])
#bibliography("biblio.bib")