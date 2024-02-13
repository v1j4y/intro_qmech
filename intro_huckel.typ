#import "template.typ": project
#import "@preview/physica:0.8.1": *
#import "@preview/cetz:0.1.2"
#import "@preview/showybox:2.0.1": showybox

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
potential energy. 

//#showybox(
//  title: "Stokes' theorem",
//  frame: (
//    border-color: blue,
//    title-color: blue.lighten(30%),
//    body-color: blue.lighten(95%),
//    footer-color: blue.lighten(80%)
//  ),
//  footer: "Information extracted from a well-known public encyclopedia"
//)[
//  Let $Sigma$ be a smooth oriented surface in $RR^3$ with boundary $diff Sigma equiv Gamma$. If a vector field $bold(F)(x,y,z)=(F_x (x,y,z), F_y (x,y,z), F_z (x,y,z))$ is defined and has continuous first order partial derivatives in a region containing $Sigma$, then
//
//  $ integral.double_Sigma (bold(nabla) times bold(F)) dot bold(Sigma) = integral.cont_(diff Sigma) bold(F) dot dif bold(Gamma) $
//]
//
//// First showybox
//#showybox(
//  frame: (
//    border-color: red.darken(50%),
//    title-color: red.lighten(60%),
//    body-color: red.lighten(80%)
//  ),
//  title-style: (
//    color: black,
//    weight: "regular",
//    align: center
//  ),
//  shadow: (
//    offset: 3pt,
//  ),
//  title: "Red-ish showybox with separated sections!",
//  lorem(20),
//  lorem(12)
//)

In order to better understand 
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
the schrodinger equation. @feynman1948

#showybox(
  footer-style: (
    sep-thickness: 0pt,
    align: right,
    color: black
  ),
  title: "Feynman's postulates of Quantum Mechanics: Postulate I",
  footer: [
    Note that the sum is over all possible paths $r_i (t)$ which 
    are possible to take from point $r_0$ to $r_1$ in time
    $t_0$ to $t_1$.
  ]
)[
The first postulate says that the total
action is the sum of the actions of individual
paths, i.e. @eq:postul1.

$
S = sum_i S[r_i (t)]
$ <eq:postul1>
]

#showybox(
  footer-style: (
    sep-thickness: 0pt,
    align: right,
    color: black
  ),
  title: "Feynman's postulates of Quantum Mechanics: Postulate II",
  footer: [
where the integral is over the region $R$ which
contains all the paths.
  ]
)[

The second postulate says that the wavefunction
$phi$ can be be expressed as an exponential
function of the position $r(t)$ and its first
deriavtive $dot(r)(t)$, i.e. @eq:postul2.

$
phi(x_k,t) = lim_(epsilon arrow.r 0) integral_R 
exp( epsilon frac(i,hbar) sum_i S[r_i (t)]) ...frac(d x_(i-1), A) frac(d x_(i-2), A)...
$ <eq:postul2>
]

== Derivation

The equation of motion describes the evolution
of the wavefunction $phi(x_(k+1),t)$ from
time $t$ to time $t+epsilon$ (@eq:deriv7).

$
phi(x_(k+1), t + epsilon) = lim_(epsilon arrow.r 0) integral_R 
exp( frac(i,hbar) sum_i S[r(t)_i]) ...frac(d x_(i), A) frac(d x_(i-1), A)...
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
phi(x_(k+1),t+epsilon) =\  exp(frac(-i epsilon V, hbar)) times integral exp(frac(i epsilon xi^2,2 hbar epsilon))[psi(x, t) - xi pdv(psi(x,t),x)  + frac(xi^2,2) pdv(psi(x,t),x,2) - ...] (d xi)/A
$ <eq:deriv11>

where, $x_(k+1) - x_k = xi$. Expanding the left hand
also around $xi$ gives.

$
phi(x_(k+1),t) + epsilon pdv(phi(x,t),t)= \ exp(frac(-i epsilon V, hbar)) times integral exp(frac(i epsilon xi^2,2 hbar epsilon))[psi(x, t) - xi pdv(psi(x,t),x)  + frac(xi^2,2) pdv(psi(x,t),x,2) - ...] (d xi)/A
$ <eq:deriv12>

The factors in the integrand on the right of
@eq:deriv12 which contain $xi$, $xi^3$ etc are
zero because they are odd integrals (@eq:deriv13).

$

phi(x_(k+1),t) + epsilon pdv(phi(x,t),t)= \ exp(frac(-i epsilon V, hbar)) times 
sqrt(2 pi hbar i /m)/A [  psi(x, t) + frac(hbar epsilon i,2 m) pdv(psi(x,t),x,2) + ...] 
$ <eq:deriv13>

Finally, equating the terms of same order in $epsilon$,
we get @eq:deriv14

$
-hbar/i pdv(psi,t) &= 1/(2 m) (hbar/i pdv(,x) )^2 psi + V(x) psi \
-hbar/i pdv(psi,t) &= H psi
$ <eq:deriv14>

The above equation can be compared to the 
time dependent schrodinger equation. The
time independent form describes stationary wavefunctions
which is given as @eq:deriv15.

$
H psi = lambda psi
$ <eq:deriv15>

= Discretization of the Hamiltonian

The position operator can be written as
shown in @eq:dis1

$
hat(q) phi(x) = x phi(x)
$ <eq:dis1>

Similarly, the momentum operator can be 
written as given in @eq:dis2.

$
hat(p) phi(x) = -i pdv(phi(x),x) " " " " (hbar"= 1 in a.u")
$ <eq:dis2>

Both these operators can be discretized on
a uniform grid of a fixed number of points 
separated by a distance $d$ as shown in
@eq:dis3.

$
hat(q) phi(x_i) &= x_i phi(x) \
hat(p) phi(x_i) &= -i frac(phi(x_(i+1)) - phi(x_(i-1)), 2d)
$ <eq:dis3>

The hamiltonian is then given by the
square of the momentum operator $p$ 
@eq:dis4.

$
hat(H) = hat(p)^2/(2m) = -1/2 pdv(,x,2)
$ <eq:dis4>

Note that the hamiltonian is real and 
depends on the coordinate $x$ in this 
one-dimensional example.

$
hat(H) = -frac(phi(x_(i+1)) + phi(x_(i-1)) - 2 phi(x_i),2 d^2)
$ <eq:dis5>

This follows from the fact that the second
derivative can be discretized using the
approximation $(f(x_(i+1)) + f(x_(i-1)) - 2 f(x_i)) / d^2$.


= 1D particle in a box

The problem of a particle in a box can be defined
as shown in @fig:partin1dbox. The position of the particle inside the box
can be defined via the wavefunction $psi(x)$. 
The particle is inside a box with infinitely large walls. 
Therefore,
the probability of finding the particle on the wall
is zero. These constitute the boundary conditions
for for the wavefunction of the particle $psi(x)$, i.e. @eq:part1.

$
psi(0) &= 0 \
psi(L) &= 0 
$ <eq:part1>

#figure(
cetz.canvas({
  import cetz.draw: *
  // Your drawing code goes here
  line((0, 0), (4.0, 0), stroke: (thickness: 4pt), name: "base")
  place-marks(
    line((0.0, -0.07), (0.0, 4.0), stroke: (thickness: 4pt)),
    (mark: ">", size: .7, pos: 1),
    fill: black
  )
  place-marks(
    line((4.0, -0.07), (4.0, 4.0), stroke: (thickness: 4pt)),
    (mark: ">", size: .7, pos: 1),
    fill: black
  )
  circle((2,1), radius: 3pt, fill: red)
  content( (-1,2), "V=∞")
  content( ( 5,2), "V=∞")
  content( ( 2,2), "V=0")
  line((0, -1.0), (4.0, -1.0), stroke: (thickness: 4pt), name: "base", mark: (begin: "<"))
  content( (2,-1.5), "L")
  circle((2,0), radius: 3pt, fill: black)
  content( (2,-0.5), "O")
}),
caption: [Particle in an infinite potential well.]
) <fig:partin1dbox>


== Finite difference equations

The schrodinger equation is given as shown in @eq:part2. 

$
hat(H)psi(x) = -1/2 grad^2 psi(x) = -1/2 frac(d^2 psi,d x^2) = lambda psi(x)
$ <eq:part2>

Here, $lambda$ represents the eigenvalues and $psi$ the 
eigenvectors of the hamiltonian. The @eq:part2 is 
a second order differential equation also known
as the one dimensional Laplace equation or 
the Poisson equation. This equation can 
be solved using numerical integration techniques.
Using Taylor series expansion of @eq:part2, one can 
obtain the finite difference formulae to evaluate the
derivative at point $x_i$.

$
lr(frac(d psi,d x) bar)_(x_i) &= lim_(epsilon arrow.r 0) ( psi(x_i + epsilon/2) - psi(x_i - epsilon/2) ) / (epsilon) \
lr(frac(d psi,d x) bar)_(x_i+epsilon/2) &= lim_(epsilon arrow.r 0) ( psi(x_i + epsilon) - psi(x_i) ) / (epsilon) \
lr(frac(d psi,d x) bar)_(x_i-epsilon/2) &= lim_(epsilon arrow.r 0) ( psi(x_i) - psi(x_i - epsilon) ) / (epsilon) \
lr(frac(d^2 psi,d x^2) bar)_(x_i) &= lim_(epsilon arrow.r 0) ( psi'(x_i + epsilon/2) - psi'(x_i - epsilon/2) ) / (epsilon) \
lr(frac(d^2 psi,d x^2) bar)_(x_i) &= lim_(epsilon arrow.r 0) ( psi(x_i + epsilon) - 2psi(x_i) + psi(x_i - epsilon) ) / (epsilon^2) 
$ <eq:part3>

The above operator on the right hand side of @eq:part3 
(called $T$ ) can be used to write the finite difference
form of the schrodinger equation @eq:part2.
This finite difference form shown in @eq:part4 can be
used to find the eigenvalues and eigenvectors of the 
schrodinger equation by dividing the segment into 
a finite numebr of uniformly distributed points.

$
T psi (x) = epsilon^2 lambda psi(x)
$ <eq:part4>

The matrix form of T is shown in @eq:part5 below 
where one can clearly see the tridiagonal form of the
Laplace operator.

#let matright(..aa) = math.mat(..aa.pos().map(row => row.map(y => {y; [$&$]})))

$
T = 1/(2 epsilon^2)  
mat(
  delim: "[",
  -1,2, -1, . , ..., ., ., .; 
  . , -1,2, -1, ..., ., ., .;
  ... , ..., ..., ...,..., ..., ..., ...;
  . , ., .,  ., ...,2, -1, .;
  . , ., .,  ., ..., -1,2, -1;
)
$ <eq:part5>

The parts not shown in the matrix above are all zeros.
The finite difference form of the schrodinger equation
can then be written as @eq:part6.

$
T psi(x) &= lambda psi(x) \
(T - lambda) psi(x) &= 0 
$ <eq:part6>

The matrix above can be diagonalized using the Lanczos
or other algorithms to obtain the eigenvalues and eigenvectors 
for the problem of particle in a box.

== Eigenvalues and eigenvectors

The eigenvectors with $n=16$ is shown in @fig:eigen1.

#figure(
   image("laplace_n16.png", width: 350pt),
   caption: [The solution of the Laplace equation with $n=16$ points.],
 ) <fig:eigen1>

As one can see, the boundary values are not 
consistent with the boundary conditions defined for the 
problem in @eq:part1. 
This is due to the finite step size (i.e. $epsilon$) and
depends on the number of points chosen for the discretization.
Increasing the number of points from $n=16$ to $n=1024$
gives a much better agreement to the bounday values defined above.

#figure(
   image("laplace_n1024.png", width: 350pt),
   caption: [The solution of the Laplace equation with $n=1024$ points.],
 ) <fig:eigen2>

Since the solutions $psi_i$ are eigenfunction of the
laplacian operator, they are by definition orthonormalized.

$
braket(psi_i , psi_i) &= 1, forall i \
braket(psi_i , psi_j) &= 0, forall i,j
$ <eq:part7>

Since we are here in real space, i.e. real coordinates $x$, 
the overlap (also known as the measure) is defined as simply
the integral @eq:part8.

$
braket(psi(x), psi(x)) = integral_0^1 psi(x)^(dagger) psi(x) d x
$ <eq:part8>

where $psi(x)^(dagger)$ is the complex conjugate of $psi(x)$.
Here, since the laplacian $T$ is hermitian (i.e. $T^(dagger) = T$),
the eigenvalues are real.

In the present case, where we have used a numerical method
to perform the integration, we can also perform the integral
numerically as shown in @eq:part9.

$
integral_0^1 psi(x)^dagger psi(x) d x = sum_k^N psi_i^dagger (k)psi_i (k)  = 1
$ <eq:part9>

where $N$ is the total number of points. Similarly,
the orthogonality constraint says that @eq:part10 holds.

$
integral_0^1 psi(x)^dagger psi(x) d x = sum_k^N psi_i^dagger (k)psi_j (k)  = 0
$ <eq:part10>

As can be easily verified using basic linear algebra.

== Probability density

The eigenvectors obtained in the previous section $psi(x)$,
can be used to plot the probability density for the various
states as given by @eq:part11.

$
rho(x) = psi(x)^dagger psi(x)
$ <eq:part11>

Numerically, this can be written as the dot product
between the wavefunction as shown in @eq:part12.

$
rho_i = sum_k psi_i (k) psi_i (k)
$ <eq:part12>

where we have assumed that the wavefunction $psi(x)$ is real.
The probability density can be shown to be positive everywhere
and defines the nodes of the state, i.e. the regions where
the probability of finding the particle is zero as shown in 
@fig:part1.

#figure(
   image("density_n1024.png", width: 350pt),
   caption: [The solution of the Laplace equation with $n=1024$ points.],
  ) <fig:part1>

The ground state has zero nodes, the first excited state has 
exactly one node, the second excited state has two nodes etc.
In general, the nodes of the function increase with increasing
energy compared to the ground state.

== Analytic solution

Analytic solution to the problem of 
particle in a one-dimensional box is of the
form shown in @eq:ana1.

$
psi_k (x) = A_k cos(k x) + B_k sin(k x)
$ <eq:ana1>

where, $k$ refers to the state in question.
The unknowns, $A_k$ and $B_k$ can be 
found using the boundray conditions and
the @eq:part2. Using the boundray conditions,
we get @eq:ana2.

$
psi(-L/2) &= psi(L/2) = 0 \
psi_k (-L/2) &= A_k cos(k L/2) - B_k sin(k L/2) = 0 \
psi_k ( L/2) &= A_k cos(k L/2) + B_k sin(k L/2) = 0
$  <eq:ana2>

which gives the final analytical solution as 
shown in @eq:ana3.

$ psi(x) := cases(
  sqrt(2/L) sin((n pi )/L x) "if" n "is even",
  sqrt(2/L) cos((n pi )/L x) "if" n "is odd",
) $ <eq:ana3>

//== Dipole moment and transition dipole moment
//
//The dipole moment operator for a given state
//is given as shown below @eq:part13.
//
//$
//angle.l x_i angle.r = integral_0^1 rho(x) x d x
//$ <eq:part13>
//
//The transition dipole moment between two states
//$i$ and $j$ is defined as shown in @eq:part14.
//
//$
//angle.l x_(i j) angle.r = integral_0^1 psi_i(x) x psi_j(x) d x
//$ <eq:part14>
//
//The table @fig:table1 shows the 
//dipole moment for the first state and the transition dipole 
//between the ground and first excited state.
//
// #figure(
//   table( columns: (auto, auto), inset: 10pt, align: center, 
// [$(i,j)$], [$angle.l x angle.r_(Psi^2)$],
// [(1,1)], [ 0.000489237   ],
// [(1,2)], [-0.180655      ],
// [(1,3)], [ 0.0           ],
// [(2,2)], [ 0.000489237   ],
// [(2,3)], [ 0.195108      ],
// ),
//   caption: [The dipole moment and transition dipole moment
//between the pairs of states shown in the first column.],
// ) <fig:table1>
 


= Second quantization

= Huckel hamiltonian

= Hubbard hamiltonian

= Double exchange hamiltonian

// Bibliography section
#pagebreak(weak: true)
#set page(header: [])
#bibliography("biblio.bib")
