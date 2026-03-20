International Journal of Modern Physics C
Vol. 23, No. 12 (2012) 1250080 (24 pages)
#.c World Scienti¯c Publishing Company
DOI: 10.1142/S0129183112500805

KLEIN TUNNELING IN THE PRESENCE
OF RANDOM IMPURITIES

S. PALPACELLI

Numidia s.r.l., via Giacomo Peroni
130, 00131 Roma, Italy
silviapalpacelli@gmail.com

M. MENDOZA* and H. J. HERRMANN†
ETH Z €urich, Computational Physics for Engineering Materials
Institute for Building Materials, Schafmattstrasse 6
HIF, CH-8093 Z €urich, Switzerland
*mmendoza@ethz.ch
†hjherrmann@ethz.ch

S. SUCCI

Istituto per le Applicazioni del Calcolo C.N.R.
Via dei Taurini, 19, 00185, Roma, Italy
and
Freiburg Institute for Advanced Studies, Albertstrasse
19, D-79104 Freiburg, Germany
succi@iac.cnr.it

Received 7 August 2012
Accepted 31 August 2012
Published 19 November 2012

In this paper, we study Klein tunneling in random media. To this purpose, we simulate the
propagation of a relativistic Gaussian wave packet through a disordered medium with randomly
distributed potential barriers (impurities). The simulations, based on a relativistic quantum
lattice Boltzmann (QLB) method, permit to compute the transmission coe±cient across the
sample, thereby providing an estimate for the conductivity (or permeability) as a function of
impurity concentration and strength of the potentials. It is found that the conductivity loss due
to impurities is signi¯cantly higher for wave packets of massive particles, as compared to
massless ones. A general expression for the loss of conductivity as a function of the impurity
percentage is presented and successfully compared with the Kozeny(cid:1)Carman law for disordered
media in classical °uid-dynamics.

Keywords: Klein paradox; disorder media; quantum lattice Boltzmann.

PACS Nos.: 66.35.+a, 03.65.Pm, 72.80.Vp.

1250080-1

S. Palpacelli et al.

1. Introduction

As opposed to classical quantum mechanics, where fermions tunneling into a barrier
are exponentially damped, relativistic scattering was shown by Klein in 19291 to
follow a very unexpected behavior: If the potential is of the order of the fermion mass
or higher, the barrier becomes virtually transparent to the fermions. The capability
of quantum wave functions to undergo zero re°ection from a potential barrier much
higher than the energy of the wave function itself is a property that relies exclusively
upon the spinorial nature of the Dirac wave function. It stands in stark contrast with
the corresponding nonrelativistic behavior, which predicts an exponential decay of
the transmission coe±cient with the di®erence V0 (cid:1) E, V0 being the height of the
barrier and E the wave function energy. Based on an analytical solution of the
it has been shown that,
scattering problem for a monochromatic plane wave,
depending on a series of geometrical and energy parameters, special angles of inci-
dence (resonant angles) provide literally zero re°ectivity2,3: the plane wave goes
completely across the barrier. Furthermore, in the presence of a random media,
transport laws similar to the ones ruling °uid motion in diluted porous media, may be
expected to apply. We refer here, e.g. to the Carman(cid:1)Kozeny law,4(cid:1)8 which relates
the permeability of a porous medium (conductivity of the medium) to the solid
concentration (impurity density).

The conductivity of two-dimensional massless fermions in disordered media has
become the object of intense studies in the literature.9(cid:1)13 The Klein(cid:1)paradox (KP)
has also been addressed in the recent literature,13(cid:1)20 but mostly in connection with
plane wave formulations. Some works have also considered localized Gaussian wave
packets,21,22 but only in connection with single and multiple ordered barriers. The
contribution of the present work to this subject encompasses the three following
directions: (i) Investigate the KP for the case of Gaussian wave packets, rather than
plane waves, for disordered samples; (ii) Discuss the viability of semi-classical
descriptions of fermionic excitations in disordered media, based on quantitative
analogies with °ows in porous media; (iii) Expose the QLB method as a new
computational tool for fermionic transport in disordered media, which might bear a
special interest for prospective implementations on parallel computers.

We wish to remark that while Klein tunneling leads to spectacular e®ects for
plane waves across one-dimensional barriers, its import for the case of two-dimen-
sional random media might turn out to be relatively milder, for two reasons. First,
because in con¯ned nonperiodic geometries and the presence of impurities,23 plane
waves might be hard to realize, hence the motivation for inspecting localized wave
functions, such as prototypical Gaussian wave packets. However, since Gaussian
wave packets, or any other localized wave function, necessarily include nonresonant
frequencies, which su®er partial re°ection, they are expected to experience less
spectacular KP e®ects than plane waves. Second, whenever the wave packet extent
exceeds the impurity size, it can split and turn around the obstacle like a classical
°uid, thereby undergoing partial transmission, with no need for any quantum

1250080-2

Klein Tunneling in the Presence of Random Impurities

tunneling through the barrier. Based on the above, it appears of interest to explore to
what extent Klein tunneling is indeed a®ecting fermionic transport within generic
disordered media.

The paper is organized as follows: ¯rst, we introduce a brief description of the
QLB method24; second, we study the case of Klein tunneling of a Gaussian wave
packet through a rectangular potential barrier. Subsequently, we present numerical
solutions of the Dirac equation in the presence of random impurities, thereby pro-
viding an estimate for the e®ects of the impurity concentration on the conductivity of
the medium, for both cases, massless and massive Dirac fermions. The simulations
are performed using a QLB model, which is also introduced as a new tool to study
transport phenomena through disordered media. Finally, we discuss and summarize
the results.

2. The Quantum Lattice Boltzmann Method

The QLB method24 is a quantum-kinetic technique originally devised to solve non-
relativistic single-body Schr€odinger equation and other related quantum problems,
such as the Gross(cid:1)Pitaevski equation for Bose(cid:1)Einstein condensates.25 Only re-
cently, it has been shown to provide a second-order accurate solver for relativistic
wave scattering and propagation.26 To forestall any confusion, we wish to clarify that
QLB is targeted to (genuine or e®ective) one-body quantum wave functions, and it
does not refer to the dynamics of statistical ensembles, such as those described by the
quantum Boltzmann equation. The name Boltzmann is attached to QLB simply as
an inheritance from its classical counterpart, the lattice Boltzmann method, which is
indeed a minimal form of Boltzmann kinetic equation to solve classical °uid-dynamic
problems from a statistical mechanics (kinetic theory) perspective.39 Since the
method is relatively new in the relativistic context, for the sake of self-containedness,
we revisit here its main technical aspects. For full details, see Refs. 27 and 28.

QLB equation was initially derived from a formal parallel between the kinetic
lattice Boltzmann equation and the single-particle Dirac equation.24,29,a For our
purpose, it proves expedient to transform the standard form of the Dirac equation
into the Majorana form, in which all matrices are real,30

½@t þ cð(cid:1)(cid:1)x@x þ (cid:2)@y (cid:1) (cid:1)z@zÞ þ i!c(cid:1)y (cid:1) igI(cid:2)  ¼ 0:

ð1Þ

This form is obtained by multiplying the standard Dirac equation on the left- and
right-hand side by the involution matrix U ¼ 2(cid:1)1=2ð(cid:1)y þ (cid:2)Þ. In the above, c is the
light speed, I is the identity operator and !c ¼ mc2=~ is the Compton frequency for a
particle of mass m, with ~ the reduced Planck's constant. The wave function   is a

a We wish to forestall any potential confusion, possibly arising from the denomination \Boltzmann" in the
QLB framework. Here \Boltzmann" stems from the formal analogy between the single-particle Dirac
equation and the Boltzmann equation of classical statistical mechanics. In this respect, QLB should be
kept distinct from the quantum Boltzmann equation, typically used to address collective quantum
transport phenomena.

1250080-3

S. Palpacelli et al.

complex four-spinors and (cid:1) and (cid:2) are the standard Dirac matrices. The last term
couples the wave function to an applied scalar potential V ðx; y; zÞ via the coe±cient
g ¼ qV =~, where q is the electric charge.30 Note that since the spin states mix-up
during propagation (spinning particles), there is no basis in which all three matrices
(cid:1)x, (cid:1)y, (cid:1)z can be simultaneously diagonalized.

Let us consider a one-dimensional version of the Dirac equation. In particular, let

Z be a unitary matrix, diagonalizing the streaming matrix (cid:1)(cid:1)z:

Z ¼

1
p
ﬃﬃﬃ
2

0

B
B
@

0 (cid:1)1
1
0
1

1
0
0 (cid:1)1 0
1
0
1
0
1
0

1

C
C
A:

ð2Þ

Applying the matrix Z to Eq. (1), the streaming matrix along z is diagonalized and
the collision matrix is also transformed accordingly

½@t þ cZ(cid:1)1ð(cid:1)(cid:1)zÞZ@z þ Z(cid:1)1ð(cid:1)c(cid:1)x@x þ c(cid:2)@y þ i!c(cid:1)y (cid:1) igIÞZ(cid:2)Z(cid:1)1  ¼ 0:

ð3Þ

Neglecting any dependence of   on the x and y coordinates, Eq. (3) may be written as
a pair of one-dimensional Dirac equations

@tu1;2 þ c@zu1;2 ¼ !cd2;1 þ igu1;2;
@td1;2 (cid:1) c@zd1;2 ¼ (cid:1)!cu2;1 þ igd1;2;

ð4Þ

for the variables ðu1; d2Þ and ðu2; d1Þ that represent the rotated wave function
Z(cid:1)1  ¼ ðu1; u2; d1; d2ÞT . The components u and d propagate up and down the z axis,
respectively, and the subscripts indicate the spin up (1) and spin down (2) states,
respectively. The system of Eq. (4) may be treated as a Boltzmann equation for a pair
of complex distribution functions u1;2 and d1;2.24 Equation (4) may thus be dis-
cretized using the same approach as in lattice Boltzmann method, i.e. by integrating
along the characteristic light-cones dz ¼ (cid:3)cdt.

The resulting system of algebraic equations reads as follows:

bu1;2 (cid:1) u1;2 ¼

bd1;2 (cid:1) d1;2 ¼ (cid:1)

1
2 ~mðd2;1 þ bd2;1Þ þ
1
2 ~mðu2;1 þ bu2;1Þ þ

1
2 i~gðu1;2 þ bu1;2Þ;
1
2 i~gðd1;2 þ bd1;2Þ;

ð5Þ

where the hat superscript (^) indicates that the wave function is evaluated at the
end-point of the corresponding streaming step, namely

bu1;2 ¼ u1;2ðz þ (cid:1)z; t þ (cid:1)tÞ;
bd1;2 ¼ d1;2ðz (cid:1) (cid:1)z; t þ (cid:1)tÞ;

u1;2 ¼ u1;2ðz; tÞ;
d1;2 ¼ d1;2ðz; tÞ:

ð6Þ

The dimensionless Compton frequency is ~m ¼ !c(cid:1)t and the dimensionless scalar
potential is ~g ¼ gðz; tÞ(cid:1)t.

1250080-4

Klein Tunneling in the Presence of Random Impurities

The pair of Eqs. (5) can be solved algebraically, delivering explicit expressions for

bu1;2 and bd1;2:

bu1;2 ¼ au1;2 þ bd2;1;
bd1;2 ¼ ad1;2 (cid:1) bu2;1;

ð7Þ

where the coe±cients a and b are

a ¼

1 (cid:1) (cid:2)=4
1 þ (cid:2)=4 (cid:1) i~g

;

b ¼

~m
1 þ (cid:2)=4 (cid:1) i~g

; (cid:2) ¼ ~m2 (cid:1) ~g2:

These coe±cients satisfy jaj2 þ jbj2 ¼ 1, so that the right-hand side of Eq. (7) cor-
responds to multiplying the rotated wave function Z(cid:1)1  ¼ ðu1; u2; d1; d2ÞT by the
unitary collision matrix

0

B
B
B
@

0

0 b
a
b 0
0
a
0 (cid:1)b a 0
0 a
0
(cid:1)b

1

C
C
C
A

:

Q ¼

ð8Þ

The streaming step propagates u1;2 upwards and d1;2 downwards, along the light-
cones given by (cid:1)z ¼ (cid:3)c(cid:1)t. Note that this unitary operation is numerically exact,
without round-o® error, because the distribution function is integrally transferred
from the source to the destination site, and no fractional transport is involved. Since
both streaming and collisions step are unitary, the overall QLB scheme evolves the
discrete wave function through a sequence of unitary operations for any value of
the discrete time step (cid:1)t. In addition, since streaming proceeds upwind only (no
centered spatial di®erences) along the discrete light-cones associated with each
component (cid:3)i, the QLB dispersion relation is automatically free from fermion-
doubling.31

This, together with the excellent e±ciency of the method,b especially on parallel
computers,32,33 makes QLB an appealing candidate for computational studies of
fermionic transport through disordered media.

The scheme extends to multiple dimensions through an operator splitting tech-
nique. Within this method, the three-dimensional Dirac equation splits into the sum
of three numbers of one-dimensional equations, each involving spatial derivatives
along one single direction. Each of the three stages representing evolution by a
time step dt is accomplished by rotating   to diagonalize the relevant streaming
matrix, taking one time step of the existing one-dimensional QLB scheme described
above, and rotating   back to its original basis. The algorithm is thus composed of
the following three steps: (1) Rotate   with X(cid:1)1, collide with X(cid:1)1 bQX, stream along
x, rotate back with X; (2) Rotate   with Y (cid:1)1, collide with Y (cid:1)1 bQY , stream along y,

b The model updates about ¯ve million lattice sites per second on a standard PC (Intel Processor of
2.4 GHz)

1250080-5

S. Palpacelli et al.

rotate back with Y ; (3) Rotate   with Z(cid:1)1, collide with Z(cid:1)1 bQZ, stream along z,
rotate back with Z . This form emphasizes the symmetry between the three steps, but
since the streaming matrix along y is already diagonal in the Majorana form, Y ¼ I is
the identity matrix. The matrix X reads as follows:

X ¼

1
p
ﬃﬃﬃ
2

0

B
B
@

1

C
C
A

(cid:1)1 0 1
0
1
0

0
1 0 (cid:1)1
0
0 1
1
1 0

ð9Þ

and the Z matrix is given in Eq. (2) above.

The collision-term splits into three parts, each of which is combined with the
corresponding streaming step. The collision matrix thus coincides, up to a unitary
transformation, with the collision matrix for the one-dimensional QLB scheme, with
a time step 1=3dt.28 In particular, bQ is given by:

0

B
B
B
B
@

bQ ¼

0

0 (cid:1)^b
^a
^b
0
^a
0
0 (cid:1)^b ^a
^b
0
0

^a

0

1

C
C
C
C
A

;

ð10Þ

where the coe±cients

^a ¼

1 (cid:1) (cid:2)3=4
1 þ (cid:2)3=4 (cid:1) i~g3

;

^b ¼

~m3
1 þ (cid:2)3=4 (cid:1) i~g3

;

are written in terms of the rescaled dimensionless Compton and potential frequencies

(cid:2)3 ¼ ~m 2

3 (cid:1) ~g

2
3;

~m3 ¼

1
3

!cdt;

~g3 ¼

1
3 gdt:

The pattern of þ and (cid:1) signs in the ^b terms on the o®-diagonal of bQ follows the same
pattern as the (cid:1)y matrix. The rotated matrices X(cid:1)1 bQX and Z(cid:1)1 bQZ have the same
sign pattern as Q, but bQ does not.

Summarizing, QLB provides a unitary, explicit algorithm for quantum wave
functions in which information propagates along classical trajectories represented by
a sequence of three one-dimensional light-cones, thereby avoiding any mixing of the
spinorial components during the streaming step. Although detailed comparisons
with other techniques remain to be developed, there are reasons to believe that such
simpli¯cation may result in enhanced computational e±ciency, especially with
parallel computers in mind. Finally, we wish to point out that the same algorithm
describes both relativistic and nonrelativistic quantum wave packets, depending on
the value of the mass m and the characteristic strength of the potential energy.

1250080-6

Klein Tunneling in the Presence of Random Impurities

3. Relativistic Gaussian Wave Packets

Our simulations will be performed in two spatial dimensions, (for more details see
Ref. 28). The propagation of a plane wave through a rectangular potential barrier
was discussed in Ref. 34. However due to the fact that it only applies to mono-
chromatic plane waves, i.e. in¯nitely extended states which may not necessarily be
realized under all experimental conditions, it is of interest to explore to what extent
are such results a®ected by the ¯nite extent of the wave function. KP e®ects using
Gaussian wave packets have been addressed in recent numerical studies.21,22,35,36
Here, we build a semi-analytical approach and test it against numerical simulations.c

For simplicity, we consider a Gaussian wave packet of the form:

 lðx; yÞ ¼

Al

ð4(cid:3)(cid:4)2Þ1=2 e(cid:1)r2=4(cid:4)2

eiðkxxþkyyÞ;

l ¼ 1; 2;

ð11Þ

where r2 ¼ x2 þ y2, A1 ¼ 1=A, A2 ¼ ei(cid:5)=A with A ¼
and (cid:5) is the angle of
the momentum vector. The rectangular box potential of height V0 and width D is
de¯ned as follows:

p

ﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃ
1 þ A 2
A 2
2

(cid:2)

V ðxÞ ¼ V0;
0;

if 0 < x < D;
elsewhere:

ð12Þ

Given the linearity of the Dirac equation and the fact that wave packets are
constituted by a Gaussian superposition of plane waves, it is natural to express the
transmission coe±cient of a Gaussian wave packet of size (cid:4) through the following
convolution:

(cid:3)

T(cid:4)ðkx; kyÞ ¼

Z

G

Sf

k (cid:1) k0
(cid:4)k

(cid:4)

T ðk

0
x; k

0
yÞdk

0
xdk

0
y;

ð13Þ

¼ k 2

F , with k 2

where Sf ¼ (cid:3)k 2
y, denotes the Fermi area and G is a Gaussian
kernel of width (cid:4)k ¼ 1=(cid:4) in momentum space. The function T ðkx; kyÞ is the trans-
mission coe±cient of a plane wave with vector k (cid:4) ðkx; kyÞ, which according to
Ref. 34, can be calculated as T ¼ 1 (cid:1) jrj2 with

x þ k 2

F

r ¼

2iei(cid:5)ðss0Þ(cid:1)1 sinðqxDÞðsin (cid:5) (cid:1) ss0 sin (cid:6)Þ
½e(cid:1)iqxD cosð(cid:5) þ (cid:6)Þ þ eiqxD cosð(cid:5) (cid:1) (cid:6)Þ(cid:2) (cid:1) 2i sinðqxDÞ

;

ð14Þ

where (cid:5) is the incidence angle, q 2
tion angle, s ¼ signðEÞ, s0 ¼ signðE (cid:1) V0Þ and E is the energy.

x ¼ ðE (cid:1) V0Þ2=~2c2, (cid:6) ¼ tan(cid:1)1ðky=qxÞ is the refrac-

Since the transmission coe±cient for a plane wave only depends on the wave
number ky, and due to the fact that the x component of the wave vector experiences a
perfect transmission, as a ¯rst-order approximation, we perform the convolution in

c We have also tested our QLB scheme against the numerical results given in Ref. 22, and found satis-
factory agreement.

1250080-7

S. Palpacelli et al.

just one dimension, ky, i.e.

T(cid:4)ðkÞ ¼

Z

kF

(cid:3)

G

(cid:1)kF

k (cid:1) k0
(cid:4)k

(cid:4)

T ðk0Þdk0;

ð15Þ

where we have de¯ned k (cid:4) ky. By setting k0 ¼ k þ q, and expanding T ðk þ qÞ around
q ¼ 0 to second-order, Eq. (15) delivers

T(cid:4)ðkÞ (cid:5) T ðkÞ þ

(cid:4) 2
2 T 00ðkÞ þ Oð(cid:4) 2
kÞ;
k

ð16Þ

where T 00 is the second derivative of T with respect to k. The above expression means
that resonant peaks (T 0ðkrÞ ¼ 0; T 00ðkrÞ < 0) are smoothed out whenever the ¯lter
width (cid:4)k, is su±ciently high, or, more precisely, (cid:4) 2
> jT 00ðkÞj=2T ðkÞ. This smoothing
is the e®ect of nonresonant wave numbers. Given that (cid:4) ¼ 1=(cid:4)k, one could readily
estimate the minimal width (cid:4) above which the secondary resonant peak would no
longer be seen by the Gaussian wave packet. However, the asymptotic expansion
given by Eq. (16) fails to represent the actual transmission coe±cient of the Gaussian
wave packet near the secondary resonant peak, the reason being that, around that
peak, a second-order expansion is grossly inaccurate because (cid:4)2T 00 (cid:5) 1 and higher
orders will be even less accurate. As a result, the convolution integral, Eq. (15), needs
to be computed.

k

3.1. Computing the convolution

To gain a quantitative sense of the dependence of the transmission coe±cient
of the Gaussian wave packet with the spatial spread (cid:4), we have numerically
computed the convolution integral of Eq. (15), for the following values (cid:4)=D ¼ 0:15;
0:31; 0:46; 0:92; 1:85, where D ¼ 100 nm is the width of the potential barrier. The
parameters are the same as in Ref. 34, namely E ¼ 0:08 eV, V0 ¼ 0:2 eV and
D ¼ 100 nm. The results are shown in Fig. 1. From this ¯gure, it is seen that, for
(cid:5) ¼ 0, T ðkrÞ ¼ T ðkF cosð(cid:5)rÞÞ goes from 1 to 0.7348, slightly over a 25% reduction.
The same ¯gure also shows that around the secondary resonance (at (cid:5) ¼ 2(cid:3)=9),
narrow wave packets with (cid:4)=D < 0:46 feature T (cid:5) 0:5, with no sign of the secondary
resonant peak. On the other hand, the secondary peak is seen to re-emerge
for (cid:4)=D > 0:92, i.e. when (cid:4) is of the order of 100 nm, comparable with the barrier
width. With (cid:4)=D ¼ 1:85, the secondary peak is recovered, but only to about 80%.
Note that, for (cid:5) ¼ (cid:3)=2, the transmission coe±cient is not zero, which is a conse-
quence of the approximation made to obtain Eq. (15) from Eq. (13). However, as
shown in Sec. 3.2, the numerical simulation of the transmission coe±cient using
QLB, shows generally a pretty satisfactory agreement with the approximation of
Eq. (15).

In order to use the plane wave approximation, one needs to ensure that the
condition (cid:4) > D is ful¯lled, which sounds pretty plausible. However, this condition is

1250080-8

Klein Tunneling in the Presence of Random Impurities

1

0.8

0.6

0.4

0.2

T

0
−2

−1

T plane wave
T filtered σ/D=0.15
T filtered σ/D=0.31
T filtered σ/D=0.46
T filtered σ/D=0.92
T filtered σ/D=1.85

0
φ

1

2

Fig. 1.
(Color online) The transmission coe±cient of a Gaussian wave packet, as computed with the
analytical convolution, Eq. (15), as a function of the incidence angle (cid:5) for (cid:4)=D ¼ 0:15; 0:31; 0:46; 0:92; 1:85.
The blue line corresponds to the un¯ltered case, (cid:4) ! 1, corresponding to a plane wave.

strongly dependent on the angle of incidence. In particular, it is far more stringent for
oblique than for head-on ((cid:5) ¼ 0) incidence. Indeed, for (cid:5) ¼ 0, (cid:4)=D (cid:5) 0:5 yields a
substantial T ¼ 0:9 for (cid:5) ¼ 0, while at (cid:5) ¼ 2(cid:3)=9, we obtain a mere T (cid:5) 0:4. At
(cid:4)=D (cid:5) 2, perfect transmission, T ¼ 1, is practically recovered at (cid:5) ¼ 0, while for
(cid:5) ¼ 2(cid:3)=9, T (cid:5) 0:8, i.e. about 80% of full transmission.

We conclude that, for head-on incidence ((cid:5) (cid:5) 0), the transmission coe±cient of
Gaussian packets is still similar to the one of plane waves, as soon their extent
becomes comparable to the barrier width. On the other hand, the secondary reso-
nance, at oblique incidence, is highly a®ected by the ¯nite-size of the wave packet,
and full recovery of perfect transmission seems to require wave packet extents sig-
ni¯cantly larger than the barrier width.

3.2. Numerical simulations

The analytical expression of Eq. (15) has been compared against direct numerical
simulation of the Dirac equation, using the QLB method. In order to back-up the
previous ¯ndings, we have computed full numerical solutions of the Dirac equation
using a QLB solver. The simulations are performed on 10242; 20482 and 40962 grids,
depending on the size of the Gaussian packet.

The physical parameters are taken from Ref. 34, i.e. E ¼ 0:080 eV, V0 ¼ 0:200 eV

and D ¼ 100 nm.

Space and time units are de¯ned by the lattice spacing in space and time, (cid:1)x and
(cid:1)t ¼ (cid:1)x=c, respectively, whereas the energy unit is taken as (cid:1)E ¼ ~=(cid:1)t. The lat-
tice spacing is chosen in such a way to resolve well the barrier width, in our case ~D ¼
D=(cid:1)x ¼ 52, corresponding to (cid:1)x (cid:5) 1:92 nm. This gives (cid:1)t ¼ (cid:1)x=c (cid:5) 1:92 (cid:6) 10(cid:1)15 s
and (cid:1)E (cid:5) 0:35 eV. The following sequence of wave packets spreading, (cid:4) ¼ 24; 48; 96

1250080-9

σ/D=0.46

T filtered
T QLB

S. Palpacelli et al.

1

0.8

0.6

0.4

0.2

T

−1

0
φ

1

2

T filtered
T QLB

0
−2

1

σ/D=0.92

T

0.8

0.6

0.4

0.2

0
−2

−1

σ/D=1.85

1

0.8

0.6

0.4

0.2

T

0
−2

−1

0
φ

0
φ

1

2

T filtered
T QLB

1

2

Fig. 2.
(Color online) The transmission coe±cient of a Gaussian wave packet as a function of the
incidence angle (cid:5) for (cid:4) ¼ 24, 48 and 96 (in lattice units), corresponding to (cid:4)=D ¼ 0:46; 0:92; 1:85, as
computed via convolution (solid line) and by QLB simulations (line with dots).

has been simulated, with D ¼ 52, all in lattice units. The results of the QLB simu-
lations appear substantially in line with the prediction of the convolution integral,
i.e. they clearly show the disappearance of the secondary peak for (cid:4)=D < 0:46, and
its progressive reappearance above this threshold (see Fig. 2). Note that, di®erent

1250080-10

Klein Tunneling in the Presence of Random Impurities

(Color online) Snapshots of the wave packet density at various instants, t ¼ 0; 420; 1050 (lattice
Fig. 3.
units), for the case (cid:5) ¼ 0 (left) and (cid:5) ¼ 2(cid:3)=9 (middle) and (cid:5) ¼ (cid:3)=3 (right) for (cid:4)=D ¼ 1:85. In the middle,
as one can see, after signi¯cant distortion in the intermediate stage of the evolution, the wave packet
manages to be transmitted across the barrier to a substantial extent (T ¼ 0:76). On the other hand, at the
right, the packet is mostly bounced-back by the barrier, with transmission coe±cient as low as T ¼ 0:13.
For visualization purposes, the color bar scale has been modi¯ed independently for each ¯gure.

from the solution of the convolution integral, Eq. (15), the transmission coe±cient
measured by the simulation is zero for (cid:5) ¼ (cid:3)=2, as should be expected, but the
appearance of the second resonant peak is still retained.

In Fig. 3, we show typical snapshots of the wave packets for the cases (cid:5) ¼ 0, 2(cid:3)=9
and (cid:3)=3, for (cid:4)=D ¼ 1:85. The snapshots clearly show that, in the case (cid:5) ¼ 0, the
wave packet crosses the barrier totally unperturbed, with literally no distortion at
any stage of the evolution. In the case of oblique resonant propagation, the packet
still manages to cross the barrier to a large extent, (T ¼ 0:76), with signi¯cant
distortions in the intermediate stages of the evolution, leaving 24% of the packet
behind. Finally, in the case of oblique nonresonant propagation, (cid:5) ¼ (cid:3)=3, the packet
is mostly bounced-back by the barrier, with a transmission coe±cient as low as
T ¼ 0:13.

4. Klein{Paradox in Random Media
To gain insight into the macroscopic properties of a medium a®ected by the presence
of impurities, it is of interest to investigate the propagation of relativistic wave
packets within a disordered sample. To analyze these transport phenomena, we
simulate the propagation of a relativistic Gaussian wave packet through a two-
dimensional domain composed of three regions: an inlet region, where the wave

1250080-11

S. Palpacelli et al.

INLET

IMPURITIES

OUTLET

Fig. 4.
Gaussian wave packet through a porous medium.

(Color online) Sketch of the domain setting used in our simulations of the propagation of a

packet is positioned at the initial time t ¼ 0; the impurity region, i.e. the central part
of the domain where randomly distributed barriers (impurities) are located; and the
outlet region, which is the ¯nal region, where measurements of the transmitted wave
packets are taken. Since Klein tunneling is a single-particle e®ect, in the sequel we
shall neglect any fermion(cid:1)fermion interaction and consider the interaction of a single
wave packet with a random distribution of impurities. Given that the wave packet
loses momentum on the impurities, the permeability (conductivity) of the sample can
be de¯ned and measured, in close analogy with a classical °uid through porous
media.d

The impurity concentration is given by C ¼ Nd2=A, where N is the number of
square impurities of cross-section d2, distributed over a total area A ¼ Ly (cid:6) Lz. Here,
Ly and Lz represent the vertical and horizontal lengths of the simulation zone,
respectively. The impurities are randomly located, and have the same size, shape and
potential value. For the present simulations, d ¼ 8 and C is varied in the range
0.001(cid:1)0.05. In Fig. 4, the computational domain is sketched, periodic boundary
conditions are imposed at top and bottom boundaries, while a bounce-back condition
is enforced at the inlet, and an open boundary condition is imposed at the outlet
(so that the transmitted wave packet is not re°ected back). This implies that
the transmission coe±cient will always converge to the unit value asymptotically
in time.

We use a square lattice of size 2048 (cid:6) 512 cells, such that the regions
½0; 512Þ (cid:6) 512, ½512; 1536Þ (cid:6) 512 and ½1536; 2048(cid:2) (cid:6) 512 correspond to the inlet, im-
purity and outlet regions, respectively. The lattice spacing is chosen in such a way as
to properly resolve the smallest lengths in the problem, namely the obstacle diameter
D, as well as the wave packet extent (cid:7). The cell size is chosen to be (cid:1)x ¼ 0:96 nm,
corresponding to (cid:4) ¼ 48 lattice spacings for the spreading of the initial Gaussian
wave packet and d ¼ 8 for the obstacle side. This yields an energy E ¼ 0:117
(80 meV in physical units).

In our study, we use two values for the mass of the particles, m ¼ 0 and m ¼ 0:1
(mc2 ¼ 0:1 ~=(cid:1)t (cid:5) 0:07 eV in physical units), and vary the impurity potential and
the concentration. We note that the energy of the massive Gaussian wave packet is

d We have also performed simulations where massless and massive particles carry the same energy, and
come to the very same conclusions: massive wave packets show smaller transmittivity than massless ones.

1250080-12

Klein Tunneling in the Presence of Random Impurities

p

ﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃ
k 2
0 þ (cid:4)(cid:1)2 þ m2

given by E ¼
and (cid:4) (cid:8) 1, the width contribution can be safely neglected, to obtain

(in lattice units). Since throughout this study k0 (cid:7) 1

q

ﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃ
k 2
0 þ m2

E (cid:5)

:

ð17Þ

p

ﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃﬃ
1 þ m2=k 2
0

Since we are principally interested in the momentum losses against the impurities, in
the sequel we shall consider the case where the massless and massive particles have
the same momentum. Based on the expression (17), this means that the massive

particle is a factor

more energetic than the massless one.

Five barrier heights are considered, namely V ¼ 25, 50, 100, 200 and 285 meV.
Note that, while the ¯rst two lie below E, hence can be overcome classically, the
others can only be traversed head-on via quantum tunneling. It should be further
observed, though, that since the wave packet is wider than the single impurity, i.e.
(cid:4) > d, even in the case E < V , the wave packet can split and turn around the
obstacle like a \classical" °uid. Our results can be classi¯ed according to the energy
of the particles, the potential of the barrier and their mass as follows: weak poten-
tials, V < E (cid:1) mc2; intermediate potentials, E (cid:1) mc2 < V < E þ mc2; and strong
potentials, V > E þ mc2. The transmission coe±cient T ðtÞ is obtained by computing
T ðtÞ ¼
(cid:8)ðz; y; tÞdzdy, where (cid:8) is the wave packet density de¯ned as
(cid:8) ¼ ju1j2 þ ju2j2 þ jd1j2 þ jd2j2, with   ¼ ðu1; u2; d1; d2ÞT being the Dirac quad-
rispinor introduced in Eq. (1). For the two-dimensional system discussed in this
paper, the Dirac equation only involves two streaming matrices along z and y, re-
spectively, so that we can set (cid:1)x ¼ 0. In the absence of magnetic ¯eld, the up-down
components of the quadrispinor collapse, and the Dirac equation for the resulting bi-
spinor can be described in terms of 2 (cid:6) 2 Pauli matrices.37

z>zoutlet

R

4.1. Wave packet mass m = 0
In this ¯rst set of simulations, we ¯x m ¼ 0, and vary the impurity concentration,
C and the strength of the impurity potential, V . In Fig. 5, we ¯x the value of V
and we compare T while varying the impurity percentage, including the reference
value for the pure sample C ¼ 0. From this ¯gure, we observe that the wave
packet takes longer to regroup for high impurity concentration and high impurity
potential. This is a natural consequence of the randomness induced in the wave
function by the disordered media. However, in all cases, the complete wave packet
is reconstructed after some time, with no stagnant regions left behind. This can
be related to the momentum loss due to the presence of the impurities, and
therefore the motion of the wave packet experiences a corresponding slow-down.
in order to recover the complete wave function, the simulations
Note that,
have been performed in a longer domain. Otherwise the right-moving wave packet
would leave the outlet region too early while the left-mover is still in the domain.
In order to provide a measurement of momentum dissipation, i.e. the loss of
conductivity due to impurities, we compute the momentum transmission

1250080-13

S. Palpacelli et al.

1

0.8

0.6

0.4

0.2

0
0

1

0.8

0.6

0.4

0.2

0
0

1

0.8

0.6

0.4

0.2

T

T

T

0
0

1

V = 50 meV

C=0%
C=0.1%
C=0.5%
C=1%
C=5%

5000

Time step

10000

15000

V = 100 meV

0.5

1
Time step

1.5

V = 285 meV

C=0%
C=0.1%
C=0.5%
C=1%
C=5%

2
x 104

C=0%
C=0.1%
C=0.5%
C=1%
C=5%

2
3
Time step

4

5
x 104

Fig. 5.
(Color online) Transmission coe±cient as a function of time for the impurity potential set at
V ¼ 50, 100 and 285 meV while varying the impurity percentage C (C ¼ 0:1, 0.5, 1 and 5%) for m ¼ 0.

coe±cient as follows:

TJzðtÞ ¼

Z

z>zout

Jzðz; y; tÞdzdy;

ð18Þ

1250080-14

Klein Tunneling in the Presence of Random Impurities

where

Jz ¼ Ã†AzÃ þ Ã†A

†
zÃ;

ð19Þ

is the z-component of the current density with Az the streaming matrix along z
and Ã ¼ ðu1; u2; d1; d2ÞT the Dirac quadrispinor.

In Fig. 6, we ¯x the value of V and compare TJz, while varying the impurity
percentage. The subscript Jz denotes the transmission coe±cient due to the
z-component of the current density, Jz. As a reference, we also plot TJzðtÞ when the
impurity percentage is set to C ¼ 0. From Fig. 6 we can observe that, unlike the
density, the momentum transmission coe±cient does not saturate at unity (its value
in the inlet region at the beginning of the simulation), because momentum is irre-
versibly lost in the impurity region. Furthermore, as expected, the momentum loss
increases with increasing impurity potential and concentration.

As a characteristic quantity associated with the dynamics of the transmission
coe±cient T , in Fig. 7, we report the escape time, t0:9, i.e. the time at which the
transmission coe±cient reaches 90%, (i.e. at 90% of the wave packet is transmitted
through the obstacle region). As above, we plot t0:9 as a function of the impurity
percentage for two values of V . We notice that for high impurity concentration the
Gaussian wave packet takes longer to cross the impurity barrier. The same e®ect
occurs when the impurity potential is increased. At low impurity concentration,
C ¼ 0:001, the e®ect of the potential barrier is relatively minor, but, as the con-
centration is increased, the escape time grows approximately linearly with the barrier
voltage.

In Fig. 8, we show some representative snapshots of the ¯rst 1800 time steps of
the simulation, for impurity percentage C ¼ 0:5% and V ¼ 50 meV. Here, we can see
the way how the wave packet is scattered by the impurities, generating a plane front,
as a result of the fragmentation of the wave function due to the random obstacles.

4.2. Wave packet mass m = 0.1
Next, we repeat the same simulations for the case of massive particles, with m ¼ 0:1.
Note that, since mc2=E ¼ 0:83, the rest energy is a signi¯cant fraction of the kinetic
energy, and therefore the wave function comes in the form of a superposition of two
wave packets, both moving at the speed of light, along opposite directions and
mixing through the nonzero mass term.

In Fig. 9, we ¯x the value of V and compare T , while varying the impurity
concentration C . As a reference, we also plot T with C ¼ 0. From the results, we
observe that the wave packet takes longer to cross the impurity region than for the
case of m ¼ 0 (the time it takes to reach a unit value of the transmission coe±cient is
longer). This is due to the slow-down of the wave function as compared to the speed
of light, because of the nonzero particle mass. Note the peak in the transmission
coe±cient, once the wave packet exits from the impurity region. This is due to the
fact that TJz takes negative values in the late stage of the evolution, indicating the

1250080-15

S. Palpacelli et al.

z
J

T

1

0.8

0.6

0.4

0.2

0

C=0%
C=0.1%
C=0.5%
C=1%
C=5%

V = 50 meV

1000

2000

3000

4000

5000

6000

Time step

C=0%
C=0.1%
C=0.5%
C=1%
C=5%

V = 100 meV

z
J

T

1

0.8

0.6

0.4

0.2

0

1000

2000

3000

4000

5000

6000

Time step

C=0%
C=0.1%
C=0.5%
C=1%
C=5%

V = 285 meV

z
J

T

1

0.8

0.6

0.4

0.2

0

1000

2000

3000

4000

5000

6000

Time step

(Color online) Momentum transmission coe±cient TJz as a function of time for the impurity
Fig. 6.
potential set at V ¼ 50, 100 and 285 meV while varying the impurity percentage C (C ¼ 0:1, 0.5, 1 and
5%) for m ¼ 0.

prevalence of the left-moving component of the wave packet once the right-moving
one has left the domain.

We compute the momentum transmission coe±cient using Eqs. (18) and (19). In
Fig. 10, we ¯x the value of V and compare TJz while varying the impurity

1250080-16

Klein Tunneling in the Presence of Random Impurities

x 104

V=50 meV
V=100 meV
V=200 meV
V=285 meV

3.5

3

2.5

2

1.5

1

0.5

p
e
t
s

e
m
T

i

0
0

1

2

3

4

5

C

(Color online) Time at which 90% of the wave packet has been transmitted, t0:9, as a function of
Fig. 7.
the impurity percentage for ¯xed values of V and m ¼ 0. The potential barriers are as follows: V ¼ 50, 100,
200 and 285 meV. The impurity percentage values are C ¼ 0:1%, 0.5%, 1% and 5%.

(Color online) Wave packet density (cid:8) at times ¼ 0, 900, 1500 and 1800 (lattice units) for the

Fig. 8.
simulation performed with impurity percentage C ¼ 0:5% and V ¼ 50 meV.

percentage. As a reference, we also plot TJzðtÞ when the impurity percentage is set to
zero. Note that, as expected, due to the inertia when the mass is increased, the curve
of the momentum transmission becomes wider than for the case of zero mass,
re°ecting the fact that the wave packet takes longer to move across the impurity
region. In addition, the maximum momentum is smaller than for the case of zero

1250080-17

S. Palpacelli et al.

T

1

0.8

0.6

0.4

0.2

0
0

1

T

0.8

0.6

0.4

0.2

0
0

1

T

0.8

0.6

0.4

0.2

V = 50 meV

C=0%
C=0.1%
C=0.5%
C=1%
C=5%

1

2
3
Time step

4

5
x 104

V=100 meV

C=0%
C=0.1%
C=0.5%
C=1%
C=5%

1

2
3
Time step

4

5
x 104

V=285 meV

C=0%
C=0.1%
C=0.5%
C=1%
C=5%

2
3
Time step

4

5
x 104

0
0

1

Fig. 9.
(Color online) Transmission coe±cient as a function of time for the impurity potential set at
V ¼ 50, 100 and 285 meV while varying the impurity percentage (C ¼ 0:1, 0.5, 1 and 5%) for m ¼ 0:1.

mass, which indicates higher momentum losses. Thus, a nonzero mass of the (quasi)-
particles, results in higher momentum losses. Also to be noted, are the negative
values of TJz in the late stage of the evolution, indicating the presence of a left-
moving component, most likely due to a spurious re°ection at the outlet boundary.

1250080-18

Klein Tunneling in the Presence of Random Impurities

C=0%
C=0.1%
C=0.5%
C=1%
C=5%

z
J

T

1

0.8

0.6

0.4

0.2

0

−0.2

V = 50 meV

1000

2000

3000

4000

5000

6000

Time step

1

0.8

0.6

0.4

0.2

0

z
J

T

−0.2

1000

z
J

T

1

0.8

0.6

0.4

0.2

0

−0.2

C=0%
C=0.1%
C=0.5%
C=1%
C=5%

V = 100 meV

2000

3000

4000

5000

6000

Time step

C=0%
C=0.1%
C=0.5%
C=1%
C=5%

V = 285 meV

1000

2000

3000

4000

5000

6000

Time step

(Color online) Momentum transmission coe±cient TJz as a function of time for the impurity
Fig. 10.
potential set at V ¼ 50, 100 and 285 meV while varying the impurity percentage (C = 0:1, 0.5, 1 and 5%)
for m ¼ 0:1.

In Fig. 11, we show selected snapshots from the ¯rst 1800 time steps of the
simulation for impurity percentage C ¼ 0:5% and V ¼ 50 meV. From this ¯gure, we
observe that a portion of the wave packet gets \trapped," moving at lower speed,
within the impurity medium, while another portion manages to move out faster.

1250080-19

S. Palpacelli et al.

Fig. 11.
simulation performed with impurity percentage ¼ 0:5% and V ¼ 50 meV and with m ¼ 0:1.

(Color online) Wave packet density (cid:8) at times 0, 900, 1500 and 1800 (lattice units) for the

4.3. Momentum transmission coe±cient TJz

In order to summarize the results obtained in the previous sections, we inspect the
maximum of the transmission coe±cient TJz in Figs. 6 and 10, as a function of the
impurity potential and concentration, for three di®erent values of mass, m ¼
0; 0:05; 0:1 (see Fig. 12). These data summarize the loss of momentum, hence resis-
tivity, due to the random impurities, formally measured by:

(cid:9)ðC; V Þ ¼ maxðTJz

ðC; V ÞÞ:

ð20Þ

From these ¯gures, we observe that at high impurity concentration, C ¼ 0:05 and a
barrier V ¼ 100 meV, the relativistic wave packet loses about 50% of its momentum, as
compared with the case of a pure sample (C ¼ 0). At the same concentration, a massive
wave packet with m ¼ 0:1, would lose more than 80%, indicating a signi¯cant drop of
transmissivity due to inertia. At low impurity level, C ¼ 0:001, both massless and
massive wave packets show a mild reduction of transmittivity, below 10%.

Let us now de¯ne the following \transmittance":

(cid:4)ðC; V Þ (cid:4)

(cid:9)
1 (cid:1) (cid:9)

:

ð21Þ

This de¯nition allows one to draw a quantitative parallel with the concept of per-
meability of a classical °uid moving through a porous medium. That is, when the
transmittance is unity, the conductivity goes formally to in¯nity, whereas zero
transmittance connotes zero conductivity.

1250080-20

Klein Tunneling in the Presence of Random Impurities

1

0.9

0.8

0.7

0.6

0.5

0.4

0.3

0.2

0.1

)

V

,

C

(
η

V = 285 meV

V = 25 meV

V = 50 meV

V = 100 meV

V = 200 meV

0

0.02

0.04

C

0.06

0.08

0.1

1

0.9

0.8

0.7

0.6

0.5

0.4

0.3

0.2

0.1

0

1

0.8

0.6

0.4

0.2

0
0

)

V

,

C

(
η

)

V

,

C

(
η

V = 25 meV

V = 50 meV

V = 100 meV

V = 200 meV

V = 285 meV

0.02

0.04

C

0.06

0.08

0.1

V = 25 meV

V = 50 meV

V = 285 meV

V = 100 meV

V = 200 meV

0.02

0.04

C

0.06

0.08

0.1

(Color online) Maximum value of TJz as a function of the impurity percentage for each value
Fig. 12.
of the impurity potential V ¼ 50(cid:1)285. For three values of the mass, m ¼ 0 (top), 0.05 (middle), 0.1
(bottom).

1250080-21

S. Palpacelli et al.

Table 1. Set of parameters that has been obtained by ¯tting the
numerical results for (cid:4) using Eq. (21).

V (meV)

25

50

100

200

285

m ¼ 0

m ¼ 0:05

m ¼ 0:1

A
n
(cid:4)0
A
n
(cid:4)0
A
n
(cid:4)0

1.09
0.8
0.51

0.68
0.84
0

0.27
0.89
0

0.26
0.85
0.23

0.16
0.88
0

0.053
0.96
0

0.046
0.98
0.17

0.03
0.99
0

0.011
1.04
0

0.017
0.97
0

0.009
1.01
0

0.0053
1.01
0

0.0097
0.94
0

0.005
1.01
0

0.0039
1.00
0

Using Eqs. (20) and (21), we have found that the numerical results are satisfac-

torily ¯tted by the following analytical expression:
ð1 (cid:1) CÞnþ1
Cn

(cid:4)ðC; V Þ ¼ A

þ (cid:4)0;

ð22Þ

where A; n; (cid:4)0 are ¯tting parameters, which depend on the strength of the potential and
the mass of the particles. In Fig. 12, we report the results of the ¯tting (solid line),
showing good agreement with the numerical data. We have plotted (cid:9) instead of (cid:4), in
order to avoid the divergence at C ¼ 0. The values of the parameters can be found in
Table 1. From this table, we appreciate that the residual (cid:4)0, is zero when the mass is
di®erent from zero, which points to this minimum permeability (conductivity) as to a
property of massless particles. On the other hand, massive particles show a closer ad-
herence to the Kozeny(cid:1)Carman law, in the context of classical °uid-dynamics,4(cid:1)8 where
no residual conductivity is observed at C ¼ 1. Also, note that for low potential barriers,
the exponent is around n (cid:5) 0:85, while for intermediate and strong potentials it is near
n (cid:5) 1, i.e. the value it takes for classical °uid-dynamics in a dilute disordered medium.
Thus, for strong potentials, the classical analogy shows satisfactory results, while for
intermediate and weak potentials, it presents deviations, typically of the order of 15%.
Finally, we observe that the case m ¼ 0 shows a signi¯cantly higher transmission than
the corresponding data with m > 0, which is due to the higher momentum losses in the
impurity region. It appears plausible to interpret the surplus of relativistic conductivity,
especially for the three cases with energy E < V , as an indirect manifestation of Klein
tunneling. Indeed, Klein tunneling requires V > mc2, hence it becomes less e®ective with
increasing mass. Di®erently rephrased, in the massive case the Dirac equation naturally
reduces to the Schr€odinger equation, which does not support any Klein tunneling.

5. Conclusions and Discussion

In this paper, we have performed a numerical study of a relativistic Gaussian
wave packet propagating through a disordered medium, which we modeled as a set of
randomly located potential barriers.

1250080-22

Klein Tunneling in the Presence of Random Impurities

From the numerical results, we conclude that for high concentration of impurities,
the wave packet presents higher losses in momentum. Furthermore, for a given
impurity concentration, by increasing the potential of each impurity, we also ¯nd a
loss of momentum. Systems with massive excitations are also studied. A nonzero
mass is found to produce higher losses of momentum in the impurity region. The
actual numerical values show that at high impurity concentration, C ¼ 0:05, the
wave packet loses more than half of its momentum with barriers of 100 meV and up
to 85% with V ¼ 285 meV. At low concentrations, C ¼ 0:001, however, the losses are
much milder, going from about 5(cid:1)20%, for V ¼ 100(cid:1)285 meV, respectively.

These data can be regrouped into an analytical expression, which bears a strong
similarity with the permeability of porous media, as a function of the porosity. We
have estimated the value of the conductivity from the transmission coe±cient and
¯tted it by using the Carman(cid:1)Kozeny law for porous media, relating the perme-
ability with the concentration of impurities. We have found that this analogy works
pretty well for the massive case, which shows no residual conductivity and a scaling
exponent pretty close to unity. On the other hand, the massless case shows a residual
conductivity. Moreover, for weak and intermediate potential strengths, the exponent
is not unity, corresponding to a fractional Kozeny(cid:1)Carman law. On the other hand,
for strong potentials, the exponent one is recovered to a good accuracy, bringing the
results closer to the analogy with classical °uids.38 The applicability of this classical
analogy indicates that, at least for the parameter set investigated in this paper,
quantum tunneling is not the dominant transport mechanism, as compared to the
semi-classical dynamics of the wave function, which can turn around the obstacles in
a similar way as a classical °uid would do.

Finally, as a byproduct, we have introduced a new single-particle tool to model
fermionic transport through disordered media, namely the QLB method. The present
QLB solves the one-particle Dirac equation, but future extensions to include e®ective
many-body interactions can also be conceived. QLB shares a remarkable computa-
tional e±ciency, especially on parallel computers, and easy handling of complex
geometries with its well-established classical LB counterpart. As a result, it is
hoped and expected that the present model can make a contribution to the
computational study of transport phenomena in any physical system governed by
the Dirac equation.

Acknowledgments

The authors are grateful for the ¯nancial support of the Eidgen€ossische Technische
Hochschule Zürich (ETHZ) under Grant No. 06 11-1.

References

1. O. Klein, Z. Phys. A, Hadrons Nuclei 53, 157 (1929).
2. K. S. Novoselov et al., Nature 438, 197 (2005).
3. K. Novoselov et al., Science 306, 666 (2004).

1250080-23

S. Palpacelli et al.

4. R. R. Rumer and P. A. Drinker, Proc. Amer. Soc. Civil Engrg., J. Hydraulic Div. 92, 155

(1966).

5. J. Bear, Dynamics of Fluids in Porous Media (American Elsevier Publising Company,

1972).

6. J. Kozeny, Sitzungsber Akad. Wiss. 136, 271 (1927).
7. P. C. Carman, Trans. Inst. Chem. Engrg. 15, 150 (1937).
8. P. C. Carman, Flows of gases through porous media (Academic Press, 1956).
9. K. Nomura and A. H. MacDonald, Phys. Rev. Lett. 98, 076602 (2007).
10. V. M. Galitski et al., Phys. Rev. B 76, 245405 (2007).
11. E. Rossi and S. Das Sarma, Phys. Rev. Lett. 101, 166803 (2008).
12. M. Polini et al., Phys. Rev. B 78, 115426 (2008).
13. S. Das Sarma, S. Adam, E. H. Hwang and E. Rossi, Rev. Mod. Phys. 83, 407 (2011).
14. A. F. Young and P. Kim, Nat. Phys. 5, 222 (2009).
15. N. Stander, B. Huard and D. Goldhaber-Gordon, Phys. Rev. Lett. 102, 026807 (2009).
16. V. V. Cheianov, V. Fal'ko and B. L. Altshuler, Science 315, 1252 (2007).
17. C. W. J. Beenakker, R. A. Sepkhanov, A. R. Akhmerov and J. Tworzydlo, Phys. Rev.

Lett. 102, 146804 (2009).

18. S. Ghosh and M. Sharma, J. Phys., Condens. Matter 21, 292204 (2009).
19. F. M. Zhang, Y. He and X. Chen, Appl. Phys. Lett. 94, 212105 (2009).
20. J. M. Pereira, V. Mlinar, F. M. Peeters and P. Vasilopoulos, Phys. Rev. B 74, 045424

(2006).

21. M. R. Setare and D. Jahani, Physica B, Condens. Matter 405, 1433 (2010).
22. K. Y. Rakhimov, A. Chaves, G. A. Farias, F. M. Peeters, J. Phys., Condens. Matter 23,

275801 (2011).

23. M. S. Jang, H. Kim, H. A. Atwater and W. A. Goddard, Appl. Phys. Lett. 97, 043504

(2010).

24. S. Succi and R. Benzi, Physica D 69, 327 (1993).
25. S. Palpacelli, S. Succi and R. Spigler, Phys. Rev. E 76, 036712 (2007).
26. P. Dellar and D. Lapitski, Phil. Trans. Roy. Soc. A 369, 2155 (2011).
27. P. J. Dellar, D. Lapitski, S. Palpacelli and S. Succi, Phys. Rev. E 83, 046706 (2011).
28. S. Palpacelli and S. Succi, Commun. Comput. Phys. 4, 980 (2008).
29. R. Benzi, S. Succi and M. Vergassola, Phys. Rep. 222, 145 (1992).
30. V. B. Berestetskii, L. P. Pitaevskii and E. Lifshitz, Quantum Electrodynamics, 2nd edn.

(Butterworths-Heinemann, Oxford, 1982).

31. J. Tworzydlo et al., Phys. Rev. B 78, 235438 (2008).
32. M. Bernaschi et al., Comput. Phys. Commun. 180, 1495 (2009).
33. M. Bernaschi et al., Proc. Int. Conf., High Performance Computing, Networking, Storage
and Analysis (SC), November 12(cid:1)18, 2011, Gordon Bell Award Honorable Mention.

34. M. I. Katsnelson, K. S. Novoselov and A. K. Geim, Nat. Phys. 2, 620 (2006).
35. F. M. Peeters, J. Phys., Condens. Matter 23, 275801 (2011).
36. J. M. Pereira Jr., F. M. Peeters, A. Chaves and G. A. Farias, Semicond. Sci. Technol. 25,

033002 (2010).

37. G. Baym, Lectures on Quantum Mechanics (American Westview Press, 1990).
38. R. D'Agosta and M. Di Ventra, J. Phys., Condens. Matter 18, 11059 (2006).
39. S. Succi, Eur. Phys. J. B. 64, 471 (2008).

1250080-24



