COMMUNICATIONSINCOMPUTATIONALPHYSICS Commun.Comput.Phys.
Vol.4,No.5,pp.980-1007 November2008
R A
EVIEW RTICLE
The Quantum Lattice Boltzmann Equation: Recent
Developments†
Silvia Palpacelli1, andSauro Succi2
∗
1 DipartimentodiMatematica,Universita` RomaTre,LargoSanLeonardoMurialdo
1,00146,Roma,Italy.
2IstitutoApplicazionidelCalcolo,VialePoliclinico137,00161Roma,Italy.
Received26February2008;Accepted(inrevisedversion)25June2008
Availableonline8July2008
Abstract. The derivation of the quantum lattice Boltzmann model is reviewed with
special emphasis on recent developments of the model, namely, the extension to a
multi-dimensionalformulationandtheapplicationtothe computationoftheground
stateoftheGross-Pitaevskiiequation(GPE).Numericalresultsforthelinearandnon-
linear Schro¨dinger equation and for the ground state solution of the GPE are also
presented and validated against analytical results or other classical schemes such as
Crank-Nicholson.
PACS:02.70.-c,03.65-w,03.67.Lx
Key words: Quantum lattice Boltzmann, multi-dimensions, imaginary-time model, linear and
non-linearSchro¨dingerequation,adiabaticlimit.
Contents
1 Introduction 981
2 FormalparallelbetweenLBEandDiracequation 982
3 QuantumlatticeBoltzmannequation 985
4 One-dimensionalquantumlatticeBoltzmannmodel 986
5 Extensiontotwoandthreespatialdimensions 987
6 AddingapotentialtotheqLBmodel 989
7 Imaginary-timequantumlatticeBoltzmannmodel 990
†DedicatedtoProfessorXiantuHeontheoccasionofhis70thbirthday.
∗Correspondingauthor. Emailaddresses: palpacel@mat.uniroma3.it (S.Palpacelli),succi@iac.rm.cnr.it
(S.Succi)
http://www.global-sci.com/ 980 (cid:13)c2008Global-SciencePress

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007 981
8 Numericalresults 993
9 Conclusionsandoutlook 1005
1 Introduction
Lattice Boltzmann models (LBMs) have become a competitive numerical tool for simu-
lating fluid flows over a wide range of complex physical problems [1–7]. LBMs were
initially derived from lattice gas cellular automata (LGCA). The basic idea of LGCA is
tosimulatethemacroscopicbehaviorofafluidflowbyimplementinganextremelysim-
plified model of the microscopic interactions between particles. LBMs were developed,
startingfromLGCA,intheattempttoovercometheirmajordrawbacks: statisticalnoise,
increasingcomplexityofthecollisionoperator(forthreedimensionalproblems)andhigh
viscosity(duetosmallnumberofcollisions)[1–3]. Nowadays,LBMhasconsolidatedinto
a powerfulalternative to more classical computational fluid dynamics models based on
thediscretizationoftheNavier-Stokesequationsofcontinuummechanics.
However,LBMand,ingeneral,thelatticekineticapproachhasbeenmostlyusedwith
classical(non-quantum)fluid. Nonetheless,withthetheorizationofquantumcomputers,
some authors have extended the lattice kinetic approach to quantum mechanics [8–16].
In fact, as it was first suggestedby Feynman [17], the most natural application of quan-
tum computers would be quantum mechanics [18]. The lattice kinetic approach is very
interesting in this respect, because it was shown that the so-called quantum lattice gas
cellularautomata(QLGCA)[11]canbeusedtosimulatesystemsofnonrelativisticquan-
tumparticleswithexponentialspeedupinthenumberofparticles[8].
Besidestheirhypotheticalandfutureapplicationtoquantumcomputing,theselattice
kineticmethodsforquantummechanicsareinterestingnumericalschemes,whichcanbe
implementedonclassicalcomputersretainingtheusualattractivefeaturesofLGCAand
LBM:simplicity,computationalspeed,straightforwardparallelimplementation.
Inthis paper,wewill focusontheso-called quantumlattice Boltzmann (qLB)model
proposedbySucciandBenzi[16,19]. TheqLBmodelwasinitiallyderivedfromaformal
parallel between the kinetic lattice Boltzmann equation (LBE) and the relativistic Dirac
equation. It was then shown that the non-relativistic Schro¨dingerequation ensuesfrom
the Dirac equation under an adiabatic assumption that is formally similar to the one
whichtakestheBoltzmannequationtotheNavier-Stokesequationsinkinetictheory[16].
The basic idea of the qLB model is to associate the wave functions composing the
Dirac quadrispinorwiththediscretedistributionfunctionsoftheLBE.Inonespatialdi-
mension,thisanalogyisnaturalandthequadrispinorcomponentscanbeassimilatedto
quantum particles of different typespropagating with velocities c and colliding when
±
theymeetat the same space-time location. However,in multi-dimensional formulation,
the analogy is no longer straightforward. This is mainly due to the fact that the Dirac
streaming operator is not diagonal along all the spatial directions (i.e., Dirac matrices
can not be simultaneously diagonalized). We could roughly say that, unlike classical

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007
982
particles,quantumparticlesofdifferenttypesmixupwhilepropagating(“spinningpar-
ticles”). To cope with this problem, a new step has to be included besides the classical
collision and streaming steps: a so-called “rotation” step. The rotation stepis neededto
alignthespinalongeachdirectionofpropagation[16].
Recently, such multi-dimensional version of the model has been implemented and
numerically validated [20]. Moreover, an imaginary-time versionofthemodelhas been
proposed to compute the ground state solution of the Gross-Pitaevskii equation (GPE)
[21]. In this paper, we will review the theoretical derivation of the qLB model and its
most recent developmentsand applications. Numerical results for the two-dimensional
linear Schro¨dingerequation,one-and two-dimensionalnonlinear Schro¨dingerequation
(namelytheGPE)andforthegroundstatesolutionoftheGPEarealsopresented.
| 2 Formal | parallel | between |     | LBE | and | Dirac | equation |     |
| -------- | -------- | ------- | --- | --- | --- | ----- | -------- | --- |
The quantum lattice Boltmzmann equation (qLBE) was initially derived from a formal
parallel between the kinetic lattice Boltzmann equation and the Dirac equation [16,22].
Thisassociationwassuggestedbytheinterestinganalogiesbetweenquantummechanics
and fluid mechanics which were pointed out from the early days of the formulation of
quantumtheory[23]. Forexample,itiswellknownthatthenon-relativisticSchro¨dinger
equationcan bewrittenin fluidformbysimplydefiningthequantumfluiddensityand
2
momentum as ρ= ψ and J ρu =(h¯/m)ρ∂ θ, where the complex wave function ψ
|     | |   | |   | a ≡ | a   |     | a   |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
is represented as ψ=ρ1/2exp(iθ). A similar analogy is also valid for the relativistic
Dirac equation. In this case, the quantum fluid can be seen as a mixture of particles
of four different types, since the Dirac equation describes the time evolution of a com-
|     | ψ=(ψ |     |     | )T. |     |     |     |     |
| --- | ---- | --- | --- | --- | --- | --- | --- | --- |
plex quadrispinor 1 ,ψ 2 ,ψ 3 ,ψ 4 As for classical particle motion, these quantum
particlespropagateinspaceandcollidewhenthey“meet”,but,differentlyfromclassical
particles,theygetmixedduringthestreamingstepbecausetheDiracstreamingoperator
isnotdiagonal[22].
For all their intellectual charm, it is now commonly accepted that they are only for-
mal similarities. However, they can be extremely useful for modeling purposes to for-
mulatenon-relativisticquantummechanicsintermsoffirst-order(relativistic)numerical
schemes. As we mentioned, the qLBE is based on an analogy between the LBE and the
| Diracequation. | Toclarifythispoint,weconsiderthekineticLBE |      |          |     |         |     |       |       |
| -------------- | ------------------------------------------ | ---- | -------- | --- | ------- | --- | ----- | ----- |
|                |                                            | (x+v | ∆t,t+∆t) |     | (x,t)=A |     | eq    |       |
|                |                                            | f    |          |     | f       | (f  | f ),  | (2.1) |
|                |                                            | i    | i        |     | − i     | ij  | j − j |       |
,barethediscretedistributionfunctionsalongthelatticespeedsv eq
| where f fori=1, |     |     |     |     |     |     |     | , f |
| --------------- | --- | --- | --- | --- | --- | --- | --- | --- |
| i               | ··· |     |     |     |     |     |     | i i |
aretheequilibriumdistributionfunctionsand A isthescatteringmatrix. Thisequation
ij
canbethoughofasadiscretizationofthefollowingsetofpartialdifferentialequations:
|     |     |     | ∂ f | +v ∂ | f =A | (f f eq ). |     | (2.2) |
| --- | --- | --- | --- | ---- | ---- | ---------- | --- | ----- |
|     |     |     | t i | ia   | a i  | ij j       |     |       |
|     |     |     |     |      |      | − j        |     |       |

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007
983
Itcan beshownthatEq.(2.2)isthefirst-orderequationresultingfromthemultiscale ex-
pansionprocedurecommonlyadoptedtostudythemacrodynamics ofLBE[4]. Eq.(2.2)
isformallyverysimilartotheDiracequation
mc2
+cαa∂
|     | ∂ ψ |      | ψ =i | β     | ψ,  | (2.3) |
| --- | --- | ---- | ---- | ----- | --- | ----- |
|     | t   | i ij | a j  | h¯ ij | j   |       |
whereαa fora=x,y,z and βarethestandardDiracmatrices,cisthespeedoflightandm
isthemassoftheparticle.
By projecting Eq. (2.2) upon the eigenvectors of the scattering matrix A , a set of
ij
hyperbolicequationsforthehydrodynamicfields
| ρ=∑ |     | =∑   |       | =∑  |         |     |
| --- | --- | ---- | ----- | --- | ------- | --- |
|     | f , | ρu a | v f , | P   | v v f   |     |
|     | i   |      | ia i  | ab  | ia ib i |     |
|     | i   |      | i     |     | i       |     |
are derived. From LBE theory, it is known that Navier-Stokes equations ensuefrom the
hyperbolic systemderiving from Eq. (2.2) by means of an adiabatic assumption. In this
context, adiabatic assumption means that the shear tensor S =∑ Q f (with Q =
|     |     |     |     |     | ab i iab i | iab |
| --- | --- | --- | --- | --- | ---------- | --- |
(v2/D)δ v
v ia v ib ab , where v is the norm of the lattice vectors i and D is the number of
−
spatialdimensions)isadiabaticallyenslavedtoitsequilibriumvalueinthelowKnudsen
numberlimit:
eq
|                                                   |     | ∂ S     | λ(S | S ), |        |     |
| ------------------------------------------------- | --- | ------- | --- | ---- | ------ | --- |
|                                                   | |   | t ab |≪ | ab  | − ab |        |     |
| whereλistheleadingeigenvalueofthescatteringmatrix |     |         |     |      | A [4]. |     |
ij
TheSchro¨dingerequationcanbederivedfromtheDiracequationinaformallyequiv-
β=v/c
alent adiabatic assumption valid in the non-relativistic limit 1, where v is the
≪
particle speed. To show this point, we consider, for the sake of simplicity, the one di-
mensionalversionofEq.(2.3)writtenintheMajoranaform[24],whereallthestreaming
matricesarerealvalued. Thisreads:
|     | ∂ u | +c∂ | u =ω  | d     | ,   |     |
| --- | --- | --- | ----- | ----- | --- | --- |
|     | t   | 1,2 | z 1,2 | c 2,1 |     |     |
(2.4)
=
|     | ∂ t d | 1,2 c∂ | z d 1,2 | ω c u | 2,1 , |     |
| --- | ----- | ------ | ------- | ----- | ----- | --- |
|     |       | −      | −       |       |       |     |
where u 1,2 and d 1,2 are the four wave functions composing the Dirac quadrispinor and
ω =mc2/h¯ istheComptonfrequency. Letusdefinethesymmetric/antisymmetricmodes
c
accordingtotheunitarytransformation
1
|     | φ   | =    | (u  | id ). |     |     |
| --- | --- | ---- | --- | ----- | --- | --- |
|     |     | 1±,2 | 1,2 | 2,1   |     |     |
|     |     | √2   | ±   |       |     |     |
StartingfromEq.(2.4),itiseasytocheckthatφ 1±,2 fulfillthefollowingequations:
|     | ∂ φ+ | +c∂ | φ =    | iω φ+ | ,   |       |
| --- | ---- | --- | ------ | ----- | --- | ----- |
|     | t    | 1,2 | z 1−,2 | c     | 1,2 |       |
|     |      |     | −      |       |     | (2.5) |
+
|     | ∂ t φ | 1−,2 +c∂ | z φ =iω | c φ 1−,2 | .   |     |
| --- | ----- | -------- | ------- | -------- | --- | --- |
1 ,2

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007
984
Up to now, the system is still symmetric and time and space derivative are in balance
(bothfirstorder). Asinkinetictheory,weneedtobreakthesymmetryofthishyperbolic
systemand writeit in adissipativeform, wherespatialand timederivative arenomore
in balance. Thesymmetryis brokenbychoosinga specifictime directionvia theenergy
| exp(iω | t)  |     |     |     |     |
| ------ | --- | --- | --- | --- | --- |
phase-shiftφ 1±,2→ φ 1±,2 c [25]. Withthisdefinitionofφ 1±,2 ,Eq.(2.5)transforminto:
+
|     | ∂ φ +c∂ | φ 1−,2 =0, |     |     |     |
| --- | ------- | ---------- | --- | --- | --- |
|     | t 1 ,2  | z          |     |     |     |
(2.6)
|     | ∂ φ +c∂ | φ+ =2iω | φ .    |     |     |
| --- | ------- | ------- | ------ | --- | --- |
|     | t 1−,2  | z       | c 1−,2 |     |     |
1,2
Inthenon-relativisticlimit, β=v/c 1, thefollowingadiabaticassumptionholds:
≪
|     | ∂ φ | 2ω       | φ .   |     |     |
| --- | --- | -------- | ----- | --- | --- |
|     | t   | 1−,2|≪ c | 1−,2| |     |     |
|     | |   | |        |       |     |     |
FromthesecondequationofEq.(2.6),byneglectingthetimederivative,weobtain
c
|     | φ        | ∂   | φ .  |     | (2.7) |
| --- | -------- | --- | ---- | --- | ----- |
|     | 1−,2∼2iω | z   | 1−,2 |     |       |
c
Inserting Eq. (2.7) into the first equation of Eq. (2.6), we finally obtain the Schro¨dinger
equationforafreeparticleofmassm,
h¯2
|     | ih¯∂ φ+ | =   | ∂2φ+ . |     |     |
| --- | ------- | --- | ------ | --- | --- |
|     | t       | 1,2 | z 1,2  |     |     |
−2m
The fast modes φ can be though of as “ghost” variables of the dynamics, in the sense
1−,2
thattheyare neededtopreservethecorrectsymmetries,althoughtheydonot“emerge”
at the macroscopic scale. To inspect the behavior of φ with respect to φ+ , we rewrite
1−,2
1,2
Eq.(2.5)intermsoftheenergyandmomentumoperatorsofquantummechanicsih¯∂ E,
t
→
ih¯∂ p :
− z → z
|     | +     | =mc2φ     | +    |     |     |
| --- | ----- | --------- | ---- | --- | --- |
|     | Eφ    | cp φ 1−,2 | ,    |     |     |
|     | 1 ,2− | z         | 1 ,2 |     |     |
(2.8)
φ+ = mc2φ
|     | Eφ 1−,2− | cp z | 1−,2 . |     |     |
| --- | -------- | ---- | ------ | --- | --- |
1,2 −
+mc2
In order to take the non-relativistic limit, we make the usual replacement E E ′
→
with E mc2 for β 0 [26]. This corresponds to the energy shift and the adiabatic
′
≪ →
assumption. Hence,Eq.(2.8)becomes
+
| E φ     | cp φ =0, |     |     |     |     |
| ------- | -------- | --- | --- | --- | --- |
| ′ 1 ,2− | z 1−,2   |     |     |     |     |
(2.9)
| +2mc2)φ |          | + 2mc2φ |          | +       |     |
| ------- | -------- | ------- | -------- | ------- | --- |
| (E ′    | 1−,2− cp | z φ     | 1−,2− cp | z φ =0. |     |
|         |          | 1 ,2∼   |          | 1 ,2    |     |
FromthesecondofEq.(2.9),weobtain
φ
|     |     | 1−,2|= 1v | β   |     |     |
| --- | --- | --------- | --- | --- | --- |
|     | |   | =         | .   |     |     |
φ+
2 c 2
| 1,2|

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007
985
FromastandardFourieranalysisofEq.(2.6),i.e.,taking
|     | φ+=ρ+exp[i(kz |     | t)], | =ρ −exp[i(kz |     | t)], |     |
| --- | ------------- | --- | ---- | ------------ | --- | ---- | --- |
|     |               |     | ω +  | φ −          |     | ω    |     |
|     |               | −   |      |              | −   | −    |     |
it can be checked that ω /ω 1/β2. Hence, the amplitude of φ goes to zero for
|     |     | +   |     |     |     | 1−,2 |     |
| --- | --- | --- | --- | --- | --- | ---- | --- |
|     |     | − ∼ |     |     |     |      |     |
β 0, while their frequencies increase as 1/β2. In conclusion, φ are small-amplitude,
1−,2
→
fast-oscillatingwavefunctions. Theseghostfieldsaresimilartotheghostfieldsemerging
fromthehydrodynamicLBE.However,thereisafundamentaldifference: ghostfieldsfor
hydrodynamic LBE tend to die out in a short time because of the real valued relaxation
coefficient,whileφ keeposcillatingsincetheirrelaxationcoefficientispurelyimaginary
1−,2
(time-reversibledynamics).
| 3 Quantum | lattice | Boltzmann | equation |     |     |     |     |
| --------- | ------- | --------- | -------- | --- | --- | --- | --- |
Intheprevioussection,wepointedoutanintriguinganalogybetweentheDiracequation
and the LBE: the Schro¨dinger equation can be obtained from the Dirac equation in the
samewayastheNavier-StokesequationsisderivedfromtheLBE.Thisinvitesaquanti-
tative correspondencebetweenLBEand Dirac equation. To this end, let us considerthe
| three-dimensionalDiracequationinMajoranaform.               |     |            |      |         | Thisreads: |     |       |
| ----------------------------------------------------------- | --- | ---------- | ---- | ------- | ---------- | --- | ----- |
|                                                             |     | cαx∂       | cαz∂ |         | αyψ.       |     |       |
|                                                             |     | ∂ t x +cβ∂ | y    | z ψ=    | iω c       |     | (3.1) |
|                                                             |     | −          | −    |         | −          |     |       |
| TheformalparallelbetweenEq.(3.1)andEq.(2.2)isasfollows[16]: |     | (cid:0)    |      | (cid:1) |            |     |       |
| f ψ,                                                        |     |            |      |         |            |     |       |
| • i → i                                                     |     |            |      |         |            |     |       |
L c( αx,β, αz),
v i
| • → ≡ | −   | −   |     |     |     |     |     |
| ----- | --- | --- | --- | --- | --- | --- | --- |
y
| A iω | α .  |     |     |     |     |     |     |
| ---- | ---- | --- | --- | --- | --- | --- | --- |
| ij   | c ij |     |     |     |     |     |     |
• →−
Someobservationsareinorder:
1. the distribution functions f are real valued, whereas ψ are complex wave func-
|     |     | i   |     |     | i   |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- |
tions;
2. the number of distribution functions f is related to the lattice used to discretize
i
the phase space, whereas ψ is composed by exactly four complex wave functions
independentlyfromthelattice;
3. the LBE streaming operator is diagonal along all the spatial directions, while the
Diracstreamingoperatorisnot,becauseitisnotpossibletosimultaneouslydiago-
| nalizethematrices |     | L=c( αx,β, | αz). |     |     |     |     |
| ----------------- | --- | ---------- | ---- | --- | --- | --- | --- |
|                   |     | −          | −    |     |     |     |     |
Themainproblemofthisapproachisclearlygivenbypointthreeabove. However,there
is a way out: to diagonalize each matrix of L separately in a sequence. The basic idea
is to use an operator splitting technique and hence consider three equivalent formula-
tions of the same equation, each having a diagonal streaming operator along x, y and
z, respectively. In practice, we write the one-dimensional Dirac equation with, say, the

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007
986
z-streaming operator in diagonal form. As we shall see, in this one-dimensional formu-
lation, a full correspondence with LBE is achieved. Thus, collision and streaming are
performedalongzbyusingtheone-dimensionalqLBE.Thenthesystemisrotated(i.e. a
unitary transformation Y is applied) so that the y-streaming matrix is diagonal and the
qLBEis usedalong y. Finally, the equationis transformedagain so that the x-streaming
matrix is diagonal and theqLBE is applied once again. After the threesequentialappli-
cationsoftheqLBE,theDiracquadrispinoristransformedback.
It is evident, from the sketch of this procedure, that the model is built upon the one
dimensional version of Eq. (3.1). Hence, in the following section we revise the one-
dimensionalversionofthemodel.
| 4 One-dimensional | quantum | lattice Boltzmann | model |     |     |
| ----------------- | ------- | ----------------- | ----- | --- | --- |
Let us consider Eq. (2.4), the one-dimensional Dirac equation in Majorana form. As ob-
servedin[16],thisisadiscreteBoltzmannequationforacoupleofcomplexbispinoru
1,2
and d . The propagation step consists on streaming u and d along z with speeds
| 1,2 |     |     | 1,2 1,2 |     |     |
| --- | --- | --- | ------- | --- | --- |
c, respectively, while thecollision stepis performedaccording tothescatteringmatrix
±
oftherighthandsideofEq.(2.4).
ThequantumlatticeBoltzmann(qLB)modelisobtainedbyintegratingEq.(2.4)along
the characteristics of u and d respectively and approximating the right hand side
1,2 1,2
Byassuming∆z=c∆t
| integralbyusingthetrapezoidalrule. |     |     | (light-conerule),wehave |     |     |
| ---------------------------------- | --- | --- | ----------------------- | --- | --- |
m
|     | uˆ u    | = (d +dˆ | ),  |     |     |
| --- | ------- | -------- | --- | --- | --- |
|     | 1,2 1,2 | 2,1 2,1  |     |     |     |
|     | −       | 2        |     |     |     |
(4.1)
em
dˆ
|     | 1,2 d 1,2 | = (u 2,1 +uˆ | 2,1 ), |     |     |
| --- | --------- | ------------ | ------ | --- | --- |
|     | −         | − 2          |        |     |     |
e
whereuˆ =u (z+∆z,t+∆t), dˆ =d (z ∆z,t+∆t), u =u (z,t),d =d (z,t) and
| 1,2 1,2 | 1,2 | 1,2 | 1,2 1,2 | 1,2 1,2 |     |
| ------- | --- | --- | ------- | ------- | --- |
−
m=ω ∆t is the lattice Compton frequency. Note that, in the scheme of Eq. (4.1), the
c
integration of the lhs is exact, a numerical error is introduced only by the discretization
oef the rhs integral. The system of Eq. (4.1) can be algebraically solved for uˆ and dˆ
|     |     |     |     | 1,2 | 1,2 |
| --- | --- | --- | --- | --- | --- |
yieldingtheqLBschemeinexplicitform
|     | uˆ 1,2 | =au 1,2 +bd 2,1 , |     |     |     |
| --- | ------ | ----------------- | --- | --- | --- |
(4.2)
|     | dˆ  | =ad bu , |     |     |     |
| --- | --- | -------- | --- | --- | --- |
|     | 1,2 | 1,2 2,1  |     |     |     |
−
where
m2/4
|     | 1      |        | m   |     |       |
| --- | ------ | ------ | --- | --- | ----- |
|     | a= −   | , b=   | .   |     | (4.3) |
|     | 1+m2/4 | 1+m2/4 |     |     |       |
|     | e      |        | e   |     |       |
The scheme of Eq. (4.2) is a lattice Boltzmann equation in matrix form [3], where the
|     | e   |     | e   |     |     |
| --- | --- | --- | --- | --- | --- |

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007
987
collisionstepisperformedbyapplyingtheunitarycollisionmatrix
|     |     |     |     | a   | 0 0 | b   |     |     |       |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ----- |
|     |     |     |     | 0   | a b | 0   |     |     |       |
|     |     |     |     | Q= |     | .  |     |     | (4.4) |
|     |     |     |     | 0   | b a | 0   |     |     |       |
−
|     |     |     |     |  b | 0 0 | a  |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     |     |     | −  |     |    |     |     |     |
|     |     |     |     |    |     |    |     |     |     |
Notethat,duetounitarityofthecollisionmatrixQ,theqLBmethodoffersunconditioned
stability with the size of the time step/mesh spacing (making sure that the light-cone
relationisfulfilled). However,itsaccuracyissubjectedtotheconditionω ∆t=∆z/λ 1,
|     |     |     |     |     |     |     |     | c   | B   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
≤
where λ =c/ω is the De Broglie wavelength of the particle. Since the time stepscales
|     | B   | c   |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
linearlywiththemeshspacing,thegridresolutioncanbeincreasedwithoutsufferingthe
timestepcollapsetypicalofclassical explicitschemeswheretheCFLstabilitycondition,
∆t <(2m/h¯)∆z2,holds. Ontheotherhand,alackofadiabaticitycouldoccurforω ∆t 1,
c
≪
and this effect must be carefully watched, while decreasing the lattice spacing, in order
topreservethevalidityofqLB.
| 5   | Extension | to  | two and | three spatial |     | dimensions |     |     |     |
| --- | --------- | --- | ------- | ------------- | --- | ---------- | --- | --- | --- |
AswementionedinSection3,theextensiontotwoandthreespatialdimensionsrequires
the introduction of a “rotation” step, besides the usual collision and propagation steps.
Thisisduetothefactthatwesplitthestreamingoperatorandapplythreetimestheone-
dimensional qLB model. In particular, assume we start from a formulation of Eq. (3.1),
forwhichthez-streamingmatrixisdiagonal,andweapplythe1D-qLBschemealongz:
|                  |     |           |     | ψ(P ,t+∆t)=S                                        | Qψ(P,t), |     |     |     |     |
| ---------------- | --- | --------- | --- | --------------------------------------------------- | -------- | --- | --- | --- | --- |
|                  |     |           |     | z                                                   | z        |     |     |     |     |
| whereP=(x,y,z),P |     | =P+∆zkˆ,S |     |                                                     |          |     |     |     |     |
|                  |     | z         |     | z isthez-streaminegoperatorandQisthecollisionmatrix |          |     |     |     |     |
(as weshallsee,thecollision matrix isnotexactly equalto Q ofEq.(4.4), thisis duetoa
factor 1/D emerging fromtheoperatorsplittingprocedure). Nowe,we rotatethesystem
so that it is aligned along y, and we apply, to the transformed equation, the 1D-qLB
schemealongy:
|     |     | (P ,t+∆t)=S |     | Qyψ(P ,t+∆t), |      | =Yψ, | Qy=Y | 1QY, |     |
| --- | --- | ----------- | --- | ------------- | ---- | ---- | ---- | ---- | --- |
|     |     | ψ yz        | y   | z             | with | ψ    |      | −    |     |
|     |     | y           |     |               |      | y    |      |      |     |
where P =P+∆yjˆ+∆zkˆ, S eis the streaming operator along y.eThe systeem is rotated
|     | yz  |     |     | y   |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
againsothatitisalignedalong xandthe1D-qLBisappliedforthelasttime:
|     | ψ   | (P ,t+∆t)=S | Qxyψ | (P ,t+∆t), | with | ψ =Xψ | , Qxy=X | 1QyX, |     |
| --- | --- | ----------- | ---- | ---------- | ---- | ----- | ------- | ----- | --- |
|     | xy  | xyz         | x    | y yz       |      | xy    | y       | −     |     |
=P+∆xiˆ+∆yjˆ+e∆zkˆ,S isthestreamingoperatoralongx.eFinally,theeupdated
whereP
|                        | xyz |     |     | x     |     |     |     |     |     |
| ---------------------- | --- | --- | --- | ----- | --- | --- | --- | --- | --- |
|                        |     |     | ψ=Y | 1X 1ψ |     |     |     |     |     |
| valueistransformedback |     |     |     | − − . |     |     |     |     |     |
xy

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007
988
Inthefollowing,wediscussthemathematicaldetailsofthisprocedureintwodimen-
sions,for thethree-dimensionalcase werefer totheoriginal reference[20], howeverthe
extensionisstraightforward.
Thestartingpointisthetwo-dimensionalversionofEq.(3.1)
|                                                |     |     |     |         | cαz∂   |         |     | αyψ, |     |     |     |
| ---------------------------------------------- | --- | --- | --- | ------- | ------ | ------- | --- | ---- | --- | --- | --- |
|                                                |     |     |     | ∂ t     | +cβ∂ y | z       | ψ=  | iω c |     |     |     |
|                                                |     |     |     |         | −      |         | −   |      |     |     |     |
| weapplytothisequationtheunitarytransformationZ |     |     |     | (cid:0) |        | (cid:1) |     |      |     |     |     |
|                                                |     |     |     |         |        | 0 1     | 0   | 1    |     |     |     |
−
|     |     |     |     |     | 1   | 1 0 | 1   | 0   |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     |     |     | Z=  |    |     | −   | ,  |     |     |     |
|     |     |     |     |     | √2  | 0 1 | 0   | 1   |     |     |     |
|     |     |     |     |     | 1  | 0   | 1   | 0  |     |     |     |
|     |     |     |     |     |    |     |     |    |     |     |     |
|     |     |     |     |     |    |     |     |    |     |     |     |
sothatthez-streamingmatrixoperatorbecomesdiagonal. We,thus,obtainthefollowing
equivalentproblem
|     |     |     |     |                     | +cAz∂ +cAy∂ |        |         |     |     |     |       |
| --- | --- | --- | --- | ------------------- | ----------- | ------ | ------- | --- | --- | --- | ----- |
|     |     |     |     | ∂                   |             |        | ψ=ω     | Cψ  |     |     |       |
|     |     |     |     |                     | t z         |        | y       | c   |     |     | (5.1) |
|     |     |     |     | ((cid:0) ψ(z,y,0)=ψ |             | (z,y), |         |     |     |     |       |
|     |     |     |     |                     |             | 0      | (cid:1) |     |     |     |       |
where
|     | 1   | 0 0 | 0   |     | 0   | 0   | 1   | 0   | 0   | 0   | 0 1 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
−
|      | 0   | 1 0 | 0   |         | 0   | 0   | 0   | 1    | 0   | 0   | 1 0  |
| ---- | --- | --- | --- | ------- | --- | --- | --- | ---- | --- | --- | ---- |
| Az= |     |     |     | , Ay= |     |     |     | − , | C= |     | .   |
|      | 0   | 0   | 1 0 |         | 1   | 0   | 0   | 0    | 0   | 1   | 0 0  |
|      |     | −   |     |         | −   |     |     |      |     | −   |      |
|      | 0  | 0 0 | 1  |         |  0 | 1   | 0   | 0   |    | 1 0 | 0 0 |
|      |    |     | −   |        |    | −   |     |     | −  |     |     |
|      |    |     |     |        |    |     |     |     |    |     |     |
Byusingthesequentialsplittingapproach, Eq.(5.1) separatesinto twoone-dimensional
problems
|     |     |     | +cAz∂   |         | )ψn=ωcCψn, |      |        | 1)∆t < | ∆t, |     |       |
| --- | --- | --- | ------- | ------- | ---------- | ---- | ------ | ------ | --- | --- | ----- |
|     |     |     | (∂      |         |            |      | (n     |        | t   |     |       |
|     |     |     | t       |         | z 1 2      | 1    | −      |        | ≤   |     | (5.2) |
|     |     |     | (ψ n[(n | 1)∆t]=ψ | n          | 1[(n | 1)∆t], |        |     |     |       |
|     |     |     | 1       |         | 2          | −    |        |        |     |     |       |
|     |     |     |         | −       |            | −    |        |        |     |     |       |
and
|     |     |     | (∂ +cAy∂ |                | )ψn=ωcCψn, |     | (n  | 1)∆t < | ∆t  |     |       |
| --- | --- | --- | -------- | -------------- | ---------- | --- | --- | ------ | --- | --- | ----- |
|     |     |     | t        |                | y          |     |     |        | t   |     |       |
|     |     |     |          |                | 2 2        | 2   | −   |        | ≤   |     | (5.3) |
|     |     |     | (ψn[(n   | 1)∆t]=ψn(n∆t), |            |     |     |        |     |     |       |
|     |     |     | 2        | −              | 1          |     |     |        |     |     |       |
for n=1,2, ,N. Tostarttheprocedure,weset ψ0=ψ . Theonedimensionalproblems
|     |     |     |     |     |     |     | 2   | 0   |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
···
of Eqs. (5.2) and (5.3) can now be solved by using the 1D-qLB scheme. However, while
Az
is already in diagonal form so that the 1D-qLB scheme can be directly applied, the
same is not true for Ay. Hence, Eq. (5.3) must be transformed (rotation step)in order to
diagonalize Ay andonepossiblechoiceforthetransformationmatrixY is
|     |     |     |     |     |     | 1   | 0 0 | 1   |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
−
|     |     |     |     |     | 1   | 0   | 1 1 | 0   |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     |     |     | Y=  |    | −   |     | .  |     |     |     |
|     |     |     |     |     | √2  | 1   | 0 0 | 1   |     |     |     |
|     |     |     |     |     |    | 0   | 1 1 | 0  |     |     |     |
|     |     |     |     |     |    |     |     |    |     |     |     |
|     |     |     |     |     |    |     |     |    |     |     |     |

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007 989
Note that, the collision matrix of this 2D-qLB is not given by Eq. (4.4), because of the
factor 1/2 in the collision term of Eqs. (5.2) and (5.3). By direct calculations, it can be
easily shown that, in this case, the collision matrix is obtained from Eq. (4.4) by simply
substitutingmwithm/2inthedefinitionofaandb(seeEq.(4.3)).
e e
6 Adding a potential to the qLB model
In this section, we will show how to add a potential interaction into the model. We
will explicitly refer to the one dimensional case for the sake of simplicity. However, the
extensiontomulti-dimensionalcaseisstraightforward.
Consider the one-dimensional Dirac equation with an electrostatic potential interac-
tion
∂ u +c∂ u =ω d +igu ,
t 1,2 z 1,2 c 2,1 1,2
(6.1)
∂ d c∂ d = ω u +igd ,
t 1,2 z 1,2 c 2,1 1,2
− −
where g=qV/h¯ is the space dependent frequency coupling to the external potential V
andqistheparticleelectriccharge. Discretizingasinthefreeparticlecase,weobtain
uˆ =a u +b d ,
1,2 g 1,2 g 2,1
(6.2)
dˆ =a d b u ,
1,2 g 1,2 g 2,1
−
where
1 Ω/4 m
a = − , b = ,
g 1+Ω/4 ig g 1+Ω/4 ig
− −
e
with Ω=m2 g2 and g=g∆t. We note that, when adding a potential, the Schro¨dinger
− e e
equationisstillobtainedin theadiabatic limit ω ω ω +g butwiththeadditional
c c
| − |≤| |
constrainteof“semall”peotentialinteraction g ω
c
. Infact,byfollowingthesameproce-
| |≪
dure outlined in Section 2 to recover the Schro¨dinger equation from Dirac equation, we
obtain
h¯c2 2ω
ih¯∂ φ+ = ∂ c ∂ φ+ qVφ+
t 1,2 −2ω z 2ω +g z 1,2 − 1,2
c (cid:18) c (cid:19)
h¯2
∂2φ+ qVφ+ ,
≈−2m z 1,2− 1,2
where the latter approximation is valid only if the potential g is negligible with respect
toω .
c
6.1 The time-dependentGross-Pitaevskii equation
ItisstraightforwardtoapplytheschemeofEq.(6.2)tothesolutionoftheGross-Pitaevskii
equation(GPE),whereanonlinearself-interactionofthewavefunctionisinvolved[27].

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007
990
Atzerotemperature,thedynamicsofatrappedBose-Einsteincondensate(BEC)isde-
TheGPEforaquantumwavefunctionψ(r,t)withr=(x,y,z)T
scribedbytheGPE[28].
∈
R3 readsas
|     |     | ∂ψ(r,t) |       | h¯2          |        |           |     |
| --- | --- | ------- | ----- | ------------ | ------ | --------- | --- |
|     |     | ih¯     | =     | ∆ r+V (r)+NU | ψ(r,t) | 2 ψ(r,t), |     |
|     |     |         |       | ext          | 0      |           |     |
|     |     | ∂t      |   −2m |              | |      | | !       |     |
where U =4πh¯2a/m is the coupling strength, a is the scattering length and N is the
0
number of particles in the condensate, V (r) is the external trapping potential, usually
ext
Furthermore,thewavefunctionψ(r,t)satisfiesthenormalization
anharmonicpotential.
condition
|     |     |     |     | ψ(r,t) 2dr=1. |     |     | (6.3) |
| --- | --- | --- | --- | ------------- | --- | --- | ----- |
|     |     |     |     | R3| |         |     |     |       |
Z
The three-dimensional GPE can be reduced to two dimensions or even one dimension,
still maintaining the same form, for twoparticular choices oftheharmonic trap [29–32].
Hence,ingeneral,weconsidertheequation
|     |     | ∂ψ(r,t) |     | h¯2    |        |         |       |
| --- | --- | ------- | --- | ------ | ------ | ------- | ----- |
|     |     |         |     | (r)+NU | ψ(r,t) | ψ(r,t), |       |
|     |     | ih¯     | =   | ∆ r+V  |        | 2       | (6.4) |
|     |     | ∂t      | −2m | ext    | d |    | |       |       |
|     |     |         |     |        |        | !       |       |
forr Rd
withd=1,2,3.
∈
To apply the qLB model to this equation, we simply needs to define the following
totalpotential
|       |        |          | V(r,t)=V                                           | (r)+NU | φ+ 2, |     |     |
| ----- | ------ | -------- | -------------------------------------------------- | ------ | ----- | --- | --- |
|       |        |          |                                                    | ext    | d | | |     |     |
|       | φ+ 2   | φ+ 2+ φ+ | 2.                                                 |        |       |     |     |
| where |        |          | Addingthispotentialtothemodel,asshownintheprevious |        |       |     |     |
|       | | | ≡| | 1 | |    | 2 |                                                |        |       |     |     |
section,theslowmodesφ+ satisfytheGPE.Notethat,thecollisionmatrixisstillunitary
1,2
( a 2+ b 2=1),butthisisduetothefactthatV isevaluatedonlyattimetandposition
| g   | g   |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- |
| | | | | | |     |     |     |     |     |     |
r, i.e. it istreatedas aconstantin theintegrationoftherhsofEq.(6.2). Inparticular, the
nonlinearity is evaluated only at the previous time step t, without any further iteration
soastopreserveunitarityandtheconsequentunconditionalstability[33,34].
| 7 Imaginary-time |     |     | quantum | lattice | Boltzmann | model |     |
| ---------------- | --- | --- | ------- | ------- | --------- | ----- | --- |
The qLB-model has been recently used to compute the ground state solution of the
GPE Eq. (6.4) [21]. In this section, we will revise the derivation of the one-dimensional
imaginary-timeqLBschemeanditsapplicationtothecomputationofthegroundstateof
theGPE. Theextensionto themulti-dimensional case follows thesame line already dis-
cussedforthereal-timeqLBschemeandfordetailswerefertotheoriginalreference[21].
| 7.1 The | time-independentGross-Piatevskii |     |     |     | equation |     |     |
| ------- | -------------------------------- | --- | --- | --- | -------- | --- | --- |
ψ(r,t)
In order to find stationary solutions of Eq. (6.4), one usually sets =
exp( iµt/h¯)φ(r),whereµisthechemicalpotential[33,35].
Byinsertingtheabovewave
−

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007
991
functionintoEq.(6.4),thefollowingequationforµandφ(r)isderived
h¯2
|                      |     |     | µφ(r)=    |     |     |      | (r)+NU | φ(r) | 2 φ(r), |       |
| -------------------- | --- | --- | --------- | --- | --- | ---- | ------ | ---- | ------- | ----- |
|                      |     |     |           |     | ∆   | r+V  |        |      |         | (7.1) |
|                      |     |     |           | −2m |     | ext  |        | d |  | |       |       |
|                      |     |     |           |     |     |      |        |      | !       |       |
| withthenormalization |     |     | condition |     |     |      |        |      |         |       |
|                      |     |     |           |     |     | φ(r) | 2=1.   |      |         |       |
|                      |     |     |           |     |     | Rd|  | |      |      |         |       |
Z
This is a constrained nonlinear eigenvalue problem, any eigenvalue µ can be computed
(r)
from its corresponding eigenfunction φ, by simply multiplying Eq. (7.1) by φ ∗ and
| integratingover |     | Rd. |     |     |     |     |     |     |     |     |
| --------------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
The ground state solution is the eigenfunction φ (r) associated with the minimum
g
eigenvalueµandsatisfyingthenormalizationcondition.
Typically, φ g is found by applying a transformation known as Wick rotation, to the
time-dependent GPE. This consists of “rotating” the time axis on the complex plane so
that it becomes purely imaginary [36–39]. Let us introduce an imaginary variable τ de-
fined as τ=it. Rewriting Eq. (6.4) in terms of τ, we obtain a diffusion equation with an
absorption/emissiontermgivenbythepotential
h¯2
|     |     | h¯∂ | ψ(r,τ)= |     | ∆   | V     | (r) | NU ψ(r,τ) | 2 ψ(r,τ). | (7.2) |
| --- | --- | --- | ------- | --- | --- | ----- | --- | --------- | --------- | ----- |
|     |     |     | τ       |     |     | r ext |     | d         |           |       |
|     |     |     |         |     | 2m  | −     | −   | |         | | !       |       |
The problem is reduced to solve Eq. (7.2) under the normalization condition constraint
givenbyEq.(6.3). Infact, for τ +∞,thesolutionofthisconstrainedequationtendsto
→
astationaryprofile,whichturnsouttobethegroundstatewavefunctionoftheGPE.
| 7.2 The | one-dimensionalimaginary-time |     |     |     |     |     |     | qLBscheme |     |     |
| ------- | ----------------------------- | --- | --- | --- | --- | --- | --- | --------- | --- | --- |
The imaginary-time qLB scheme is obtained by applying the Wick rotation to the Dirac
Eq. (6.1) from which qLB starts from. The basic idea is to rewrite the Dirac equation
(writteninMajoranaform)withrespecttotheimaginaryvariable τ=it.
In the one-dimensional case, we consider Eq. (6.1) and apply the Wick rotation, to
yield:
|     |     |     |     | ∂ u | ic∂   | u     | = iω | d +gu | ,   |       |
| --- | --- | --- | --- | --- | ----- | ----- | ---- | ----- | --- | ----- |
|     |     |     |     | τ   | 1,2 − | z 1,2 | −    | c 2,1 | 1,2 | (7.3) |
|     |     |     |     | ∂ d | +ic∂  | d     | =iω  | u +gd | .   |       |
|     |     |     |     | τ   | 1,2   | z 1,2 | c    | 2,1   | 1,2 |       |
Now, let ∆τ=i∆t be the imaginary-time discretization step, while ∆z= ic∆τ=c∆t is,
−
as usual, the spatial step. Integrating Eq. (7.3) along the characteristics of u and d
1,2 1,2
respectivelyandapproximatingtherhsintegralwiththetrapezoidalrule,weobtain
|     |     |     |     |     |     | m   |     | g   |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
+dˆ
|     |     |     | uˆ  | u     | =   | i (d  |     | )+ (uˆ | +u ),   |     |
| --- | --- | --- | --- | ----- | --- | ----- | --- | ------ | ------- | --- |
|     |     |     | 1,2 | − 1,2 | −   | 2 2,1 | 2,1 | 2      | 1,2 1,2 |     |
(7.4)
|     |     |     |     |     | m   | e       |     | g e  |       |     |
| --- | --- | --- | --- | --- | --- | ------- | --- | ---- | ----- | --- |
|     |     |     | dˆ  | d   | =i  | ( u +uˆ | )+  | ( dˆ | +d ), |     |
|     |     |     | 1,2 | 1,2 |     | 2,1     | 2,1 | 1,2  | 1,2   |     |
|     |     |     |     | −   | 2   |         |     | 2    |       |     |
|     |     |     |     |     | e   |         |     | e    |       |     |

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007
992
where uˆ =u (z+∆z,τ+∆τ), dˆ =d (z ∆z,τ+∆τ), u =u (z,τ), d =d (z,τ),
| 1,2 | 1,2 |     | 1,2 | 1,2 |     |     | 1,2 1,2 | 1,2 | 1,2 |
| --- | --- | --- | --- | --- | --- | --- | ------- | --- | --- |
|     |     |     |     | −   |     |     |         |     | dˆ  |
m=ω ∆t and g=g∆t. ThesystemofEq. (7.4) is algebraically solvedfor uˆ and , to
| c   |     |     |     |     |     |     |     | 1,2 | 1,2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
yieldthefollowingexplicitscheme
| e   | e   |     |        |      |     |     |     |     |     |
| --- | --- | --- | ------ | ---- | --- | --- | --- | --- | --- |
|     |     |     |        | =aiu | bid |     |     |     |     |
|     |     |     | uˆ 1,2 | 1,2  | 2,1 | ,   |     |     |     |
g − g
|     |     |     | dˆ  | =aid  | +biu  | ,   |     |     |     |
| --- | --- | --- | --- | ----- | ----- | --- | --- | --- | --- |
|     |     |     | 1,2 | g 1,2 | g 2,1 |     |     |     |     |
where
|     | (1  | g/2)(1+g/2)+m2/4 |     |     |     |     |     |     |     |
| --- | --- | ---------------- | --- | --- | --- | --- | --- | --- | --- |
im
|     | ai = | −   |            |     | , bi = |          |      | .   |     |
| --- | ---- | --- | ---------- | --- | ------ | -------- | ---- | --- | --- |
|     | g    | (1  | g/2)2 m2/4 |     | g      | (1 g/2)2 | m2/4 |     |     |
|     |      |     | − −        |     |        | −        | −    |     |     |
|     |      | e   | e          | e   |        |          | e    |     |     |
Notethat, gisevaluatedattimeτ,i.e. thereisnoiterationoverthenonlinearity.
bie2=1,heencethecollisionmatreixisnoteunitary.
| Weobservethat | ai  | 2+   |      |     |     |     |     | Thisimplies |     |
| ------------- | --- | ---- | ---- | --- | --- | --- | --- | ----------- | --- |
|               | |   | g| | | g| 6 |     |     |     |     |             |     |
that the model does not automatically verify the normalization condition. This is usual
for modelswhich computethegroundstatesolutionby solving dynamical equationsin
fictitioustime,suchasEq.(7.2). Hence,thenormalization conditionmustbeimposedat
eachtimestepbydirectlyre-normalizingthewavefunction[33,35,38].
Inanalogywithreal-timeqLB,weintroducethewavefunctions
1
|     |     |     | φ =  | eωcτ(u | +id | ).  |     |     | (7.5) |
| --- | --- | --- | ---- | ------ | --- | --- | --- | --- | ----- |
|     |     |     | 1±,2 |        | 1,2 | 2,1 |     |     |       |
√2
In[21],theequationgoverningφ 1−,2 isderivedanditisshown,byinspectingitsdispersion
relation,thatφ 1−,2 fulfillsEq.(7.2). Inparticular,itisshownthatthedispersionrelationof
theequationgoverningφ 1−,2 coincideswiththedispersionrelationofEq.(7.2)apartfrom
anadditionalmodewhoseeffectistouniformlyamplify φ 1−,2 ,butthenormalization step
compensates this anomalous behavior. For the details of this procedure, we refer to the
original reference [21]. Here we will revise the simpler free-particle case (V =0 and
ext
NU =0).
1
Sinceu andd fulfillEq.(7.3),φ 1±,2 satisfythefollowingequations:
| 1,2 | 1,2 |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
φ+
|     |     |     | ∂      | ic∂ φ  | =0, |     |     |     | (7.6) |
| --- | --- | --- | ------ | ------ | --- | --- | --- | --- | ----- |
|     |     |     | τ 1,2− | z 1−,2 |     |     |     |     |       |
+ =2ω
|     |     |     | ∂ τ φ 1−,2− | ic∂ z φ | c φ | 1−,2 . |     |     | (7.7) |
| --- | --- | --- | ----------- | ------- | --- | ------ | --- | --- | ----- |
1 ,2
By multiplying (7.6) by c and taking the z derivative, multiplying Eq. (7.7) by i and de-
riving it with respect to τ and then subtracting the resulting equations, we obtain the
followinggoverningequationforφ
1−,2
|     |     |     |              | h¯2  | h¯2 |      |     |     |       |
| --- | --- | --- | ------------ | ---- | --- | ---- | --- | --- | ----- |
|     |     |     |              | ∂2φ  |     | ∂2φ  |     |     |       |
|     |     |     | h¯∂ φ 1−,2 = | 1−,2 | +   | 1−,2 | .   |     | (7.8) |
|     |     |     | τ            | 2m z | 2ω  | τ    |     |     |       |
c
The second order time derivative term drives an instability which tends to amplify φ 1−,2
while preserving its spatial profile. However, the normalization step tames the effect of

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007 993
this term. This point is clarified by the study of the dispersion relation we mentioned
above. Here,wejustwanttoobservethat,forthefree-particlecase(V =0andNU =0),
ext d
φ obeyadiffusionequationwiththecorrectdiffusioncoefficient(seeEq.(7.2)).
1−,2
Notethat, in thisimaginary-time extensionofthemodel, noadiabatic assumptionis
needed and it is not required to have “small” potential interaction. As we shall see, in
thiscasethenonlinearitycoefficientcanbesettomuchlargervaluesthaninthereal-time
case.
8 Numerical results
In this section, we will present some numerical results for the two-dimensional linear
Schro¨dinger,one-andtwo-dimensionaltime-dependentGPEandtwo-dimensionaltime-
independentGPE.
8.1 LinearSchro¨dinger equation
ForthelinearSchro¨dingerequation,wereviewsomeofthetwo-dimensionalresultspre-
sentedin[20]. One-dimensionalnumericalresultscanbefoundin[19].
Letusconsider,asinitialcondition,aminimumuncertaintywavepacket
(z z )2 (y y )2
ψ 0 (z,y)= 2π∆ 0z ∆ 0y − 1/2 exp (cid:18) − 4 − ∆2 0 0 z (cid:19) exp − 4 − ∆2 0 0 y ! exp − im(v z z+v y y) .
(cid:0) (cid:1) (cid:0) ((cid:1)8.1)
This is a wave packet centered about (z ,y ) with initial spreads ∆ , ∆ along z and y
0 0 0z 0y
respectively and propagating at speed (v ,v ). With this initial condition, the analytical
z y
solutionoftheSchro¨dingerequationforafreepropagatingparticleisgivenby[40]:
ψ (z,y,t)= 2π ∆ + it ∆ + it − 1/2 exp (z − z 0 − v z t)2
an 0z 2m∆ 0y 2m∆ −4∆2 +2it/m
(cid:20) (cid:18) 0z(cid:19)(cid:18) 0y(cid:19)(cid:21) (cid:18) 0z (cid:19)
(y y v t)2 im(v2+v2)t
0 y z y
exp − − exp im(v z+v y) exp . (8.2)
× − 4∆2 +2it/m
!
z y − 2
!
0y
(cid:0) (cid:1)
Based on this solution, the mean position (Z(t),Y(t)) and the mean spreads ∆ (t) and
z
∆ (t)evolveaccordingtotheequations
y
Z(t)=z +v t, Y(t)=y +v t, (8.3)
0 z 0 y
and
t2 1/2 t2 1/2
∆ (t)= ∆2 + , ∆ (t)= ∆2 + . (8.4)
z
(cid:20)
0z 4m2∆2
0z(cid:21)
y
"
0y 4m2∆2
0y#

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007
994
| y=y    | andt=100 |      |         | z=z andt=100 |      |
| ------ | -------- | ---- | ------- | ------------ | ---- |
|        | 0        |      |         | 0            |      |
| x 10−3 |          |      | x 10−3  |              |      |
| 1      |          |      | 1       |              |      |
|        |          | 128  |         |              | 128  |
|        |          | 256  |         |              | 256  |
|        |          | 512  |         |              | 512  |
| 0.8    |          |      | 0.8     |              |      |
|        |          | 1024 |         |              | 1024 |
| )001,  |          |      | )001,y, |              |      |
| 0.6    |          |      | 0.6     |              |      |
y,z(e 0
0
| 0.4    |          |         | z(e 0.4 |              |         |
| ------ | -------- | ------- | ------- | ------------ | ------- |
| 0.2    |          |         | 0.2     |              |         |
| 0      |          |         | 0       |              |         |
| 0 100  | 200 300  | 400 500 | 0 100   | 200 300      | 400 500 |
|        | z        |         |         | y            |         |
| y=y    | andt=300 |         |         | z=z andt=300 |         |
|        | 0        |         |         | 0            |         |
| x 10−3 |          |         | x 10−3  |              |         |
| 2      |          |         | 2       |              |         |
|        |          | 128     |         |              | 128     |
|        |          | 256     |         |              | 256     |
|        |          | 512     |         |              | 512     |
| 1.5    |          | 1024    | 1.5     |              | 1024    |
| )003,  |          |         | )003,y, |              |         |
| 1      |          |         | 1       |              |         |
y,z(e 0
0
z(e
| 0.5    |          |         | 0.5     |              |         |
| ------ | -------- | ------- | ------- | ------------ | ------- |
| 0      |          |         | 0       |              |         |
| 0 100  | 200 300  | 400 500 | 0 100   | 200 300      | 400 500 |
|        | z        |         |         | y            |         |
| y=y    | andt=500 |         |         | z=z andt=500 |         |
|        | 0        |         |         | 0            |         |
| x 10−3 |          |         | x 10−3  |              |         |
| 2.5    |          |         | 2.5     |              |         |
|        |          | 128     |         |              | 128     |
|        |          | 256     |         |              | 256     |
| 2      |          | 512     | 2       |              | 512     |
|        |          | 1024    |         |              | 1024    |
| )005,  |          |         | )005,y, |              |         |
| 1.5    |          |         | 1.5     |              |         |
y,z(e 0
0
| 1     |         |         | z(e 1 |         |         |
| ----- | ------- | ------- | ----- | ------- | ------- |
| 0.5   |         |         | 0.5   |         |         |
| 0     |         |         | 0     |         |         |
| 0 100 | 200 300 | 400 500 | 0 100 | 200 300 | 400 500 |
|       | z       |         |       | y       |         |
Figure1: Differencebetweentherealpartoftheanalyticalsolutionandthemodelresultfory=y 0 (leftcolumn)
andz=z 0 (rightcolumn)(seeEq.(8.5))andfordifferentvaluesof Nz =Ny attimest=100,t=300andt=500.
Parameters are set as follows: vz =0.02, vy =0.04, ∆ =∆ =40 and m=(1/8)h. Solid line: Nz =Ny =128;
0z 0y
Dashed line: Nz =Ny =256; Dotted line: Nz =Ny =512; Dash-dotted line: Nz =Ny =1024. Time and space
| areexpressed in lattice | units. |     |     |     |     |
| ----------------------- | ------ | --- | --- | --- | --- |
e

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007
995
Wecompare theqLBnumerical resultswith theanalytical solutionwhile decreasingthe
|     |     | ∆t=∆z=∆y |     |     | c=1 |     |     |     |
| --- | --- | -------- | --- | --- | --- | --- | --- | --- |
discretization step h (recall that in lattice units). For this numerical
≡
test,wesetthecomputationaldomainas[0,512] [0,512] andy =z =256inlatticeunits
|     |     |     |     |     |     |     | 0 0 |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
×
and the domain is discretized with N =N =128, 256, 512 and 1024 lattice points. The
|                                     |     |     |     | z y  |       |         |           |     |
| ----------------------------------- | --- | --- | --- | ---- | ----- | ------- | --------- | --- |
|                                     |     |     |     | ∆ =∆ | =40,v | =0.02,v | =0.04andm | ∆t= |
| remainingparametersaresetasfollows: |     |     |     | 0z   | 0y    | z       | y         | ω c |
≡
(1/8)h. Periodicboundaryconditionsareimposedinallthesimulations. Theerrorwith
respecttotheanalyticalsolutionEq.(8.2)iscomputedinL2normandisfoundteodecrease
from0.08 to0.009 as thegrid resolutionwas increasedfrom128 to1024 points,butwith
no clear evidence of a specific convergence rate. This could be the effect of concurrence
errorsources: thetimediscretizationerror( (h2)),thesplittingerror( (h))andthelack
O O
| ofadiabaticityform=ω |     | ∆t  | 0. InFig.1thefunction |     |     |     |     |     |
| -------------------- | --- | --- | --------------------- | --- | --- | --- | --- | --- |
|                      |     | c → |                       |     |     |     |     |     |
(φ+(z,y,t))
|     |     | e(z,y,t)= |     | (ψ (z,y,t)) |     |     |     | (8.5) |
| --- | --- | --------- | --- | ----------- | --- | --- | --- | ----- |
|     | e   |           | |ℜ  | an          | −ℜ  |     | |   |       |
takenatthecrosssectionsy=y 0 andz=z 0 forthedifferentvaluesofN z andN y isplotted
attimes100,300and500.
InTable1thepropagationvelocityandthemeanspreadofthepacketwithincreasing
resolution are also shown. For the present setting the expected velocity is v =0.02 and
z
v =0.04and thespreadat time t=500 is computedby Eq.(8.4) and is 64.03 since ∆ =
| y   |     |     |     |     |     |     |     | 0z  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
∆ =40.
0y
Table 1: Propagation velocity and spread of the packet at time t=500 for different values of Nz =Ny . The
expected values are: vz =0.02, vy =0.04and ∆ (500)=∆ (500)=64.03. Here m=(1/8)h and ∆ =∆ =40.
|     |     |      |     | z   | y       |     |         | 0z 0y |
| --- | --- | ---- | --- | --- | ------- | --- | ------- | ----- |
|     |     | N =N | v   | v   | ∆ (500) |     | ∆ (500) |       |
|     |     | z y  | z   | y   | z       |     | y       |       |
e
|     |     | 128  | 0.0175 | 0.0355 | 60.20 |     | 60.19 |     |
| --- | --- | ---- | ------ | ------ | ----- | --- | ----- | --- |
|     |     | 256  | 0.0189 | 0.0379 | 62.41 |     | 62.40 |     |
|     |     | 512  | 0.0191 | 0.0384 | 62.97 |     | 62.95 |     |
|     |     | 1024 | 0.0193 | 0.0386 | 63.11 |     | 63.09 |     |
Now that the convergence of the scheme has been shown, we set h=1 as it is usual
inlatticeBoltzmannschemes. Wewanttoshowtheabilityofthemodeltoreproducethe
meanquantitiesZ(t),Y(t),∆ (t)and∆ (t)asfunctionsoftimefollowingthetheoretical
|     |     | z   |     | y   |     |     |              |        |
| --- | --- | --- | --- | --- | --- | --- | ------------ | ------ |
|     |     |     |     |     |     | =N  | =1024 =0.05, | =0.02, |
relations of Eqs. (8.3), (8.4). In particular, we set N z y and v z v y
∆ =50, ∆ =32 and m=0.2. In Fig. 2, we compare the numerical curve of Z(t), Y(t),
| 0z 0y |     |     |     |     |     |     |     |     |
| ----- | --- | --- | --- | --- | --- | --- | --- | --- |
∆ (t)and∆ (t)withtheanalyticalfunctionsgivenbyEqs.(8.3)and(8.4).
| z y |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
As a second exampele, we consider the introduction of an harmonic oscillator. The
harmonicpotentialisgivenby
1
|     |     | V(z,y)= |     | mω2[(z | z )2+(y | y   | )2], |     |
| --- | --- | ------- | --- | ------ | ------- | --- | ---- | --- |
|     |     |         |     | 0      | 0       |     | 0    |     |
|     |     |         | 2   |        | −       | −   |      |     |
whereweareassumingforsimplicityω =ω =ω ,notisotropicharmonicpotentialsare
|     |     |     |     | z   | y 0 |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
considered in [20]. It is known that, in this case, the mean position satisfy the classical

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007
996
|     | 80  |     |     |     |     |     |     | 100 |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
∆
|     | 70  |     |     |     | (t) |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
z
|     |     |     |     |     |     |     |     | 80  |     |     | ∆ (t) |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ----- | --- |
60
y
50
60
Z(t)
40
|     | 30  |     |     |     |     |     |     | 40  |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Y(t)
20
20
10
|     | 0   |     |     |     |          |      |     | 0   |         |     |          |      |
| --- | --- | --- | --- | --- | -------- | ---- | --- | --- | ------- | --- | -------- | ---- |
|     | 0   | 200 | 400 | 600 | 800 1000 | 1200 |     | 0   | 200 400 | 600 | 800 1000 | 1200 |
|     |     |     |     | t   |          |      |     |     |         | t   |          |      |
Figure 2: Comparison between Z(t), Y(t), ∆ z (t) and ∆ y (t) and the expected curves given by Eqs. (8.3) and
(8.4) forthefollowing setting: Nz =Ny =1024, vz =0.05, vy =0.02, ∆ =50, ∆ =32and m=0.2. Solid lines
|     |     |     |     |     |     |     |     |     | 0z  | 0y  |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
representnumericalresults;dashedlinesaretheexpectedcurves. Timeandspaceareexpressedinlatticeunits.
e
equationofmotionoftheharmonicoscillatorandthisyields
|                     |     |     |        |     | v         |        |            |     | v         |     |     |       |
| ------------------- | --- | --- | ------ | --- | --------- | ------ | ---------- | --- | --------- | --- | --- | ----- |
|                     |     |     | Z(t)=z |     | + z sin(ω |        | t), Y(t)=y |     | + y sin(ω | t). |     |       |
|                     |     |     |        |     | 0         |        | 0          |     | 0         | 0   |     | (8.6) |
|                     |     |     |        |     | ω         |        |            |     | ω         |     |     |       |
|                     |     |     |        |     | 0         |        |            |     | 0         |     |     |       |
| Moreover,bysetting∆ |     |     |        | ∆   | =∆        | sothat |            |     |           |     |     |       |
|                     |     |     |        | 0   | 0z 0y     |        |            |     |           |     |     |       |
≡
1
|     |     |     |     |     |     | ω   | 0 = | ,   |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
2m∆2
0
the initial spreading is preserved all along the evolution. We want to check the ability
∆
of the model to preserve for different parameter settings. In Table 2, the results are
0
|     |     | ∆   | ∆   |     |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
shown, here z and y are the packet spreads averaged over two periods. In all of the
| simulationswesetv |     |     |     | =0.02andv |     | =0.04. |     |     |     |     |     |     |
| ----------------- | --- | --- | --- | --------- | --- | ------ | --- | --- | --- | --- | --- | --- |
|                   |     |     | z   |           | y   |        |     |     |     |     |     |     |
Table2: Averagedvariancesofthepacketalongzandyfordifferentsettingoftheparameters Nz =Ny ,mand
| ω   | . Here | vz =0.02 and | vy  | =0.04. |     |     |     |     |     |     |     |     |
| --- | ------ | ------------ | --- | ------ | --- | --- | --- | --- | --- | --- | --- | --- |
0
|     |     | N =N |     | ω     | m    |       | ∆    |       | ∆    | Expected∆ |     | e   |
| --- | --- | ---- | --- | ----- | ---- | ----- | ---- | ----- | ---- | --------- | --- | --- |
|     |     | z    | y   | 0     |      |       | z    |       | y    |           |     |     |
|     |     | 1024 |     | 8/642 | 1/16 | 64.35 | 1.33 | 64.35 | 1.33 |           | 64  |     |
|     |     |      |     | 2/322 |      |       | ±    |       | ±    |           |     |     |
|     |     | 1024 |     |       | 1/e4 | 32.25 | 0.70 | 32.27 | 0.75 |           | 32  |     |
|     |     |      |     |       |      |       | ±    |       | ±    |           |     |     |
|     |     | 512  |     | 4/322 | 1/8  | 32.16 | 0.70 | 32.17 | 0.69 |           | 32  |     |
|     |     |      |     |       |      |       | ±    |       | ±    |           |     |     |
|     |     | 512  |     | 2/322 | 1/4  | 31.87 | 0.27 | 31.87 | 0.29 |           | 32  |     |
|     |     |      |     |       |      |       | ±    |       | ±    |           |     |     |
|     |     | 512  |     | 1/162 | 1/2  | 16.01 | 0.69 | 16.02 | 0.70 |           | 16  |     |
|     |     |      |     |       |      |       | ±    |       | ±    |           |     |     |
2/162
|     |     | 256 |     |       | 1/4 | 16.05 | 0.37 | 16.05 | 0.38 |     | 16  |     |
| --- | --- | --- | --- | ----- | --- | ----- | ---- | ----- | ---- | --- | --- | --- |
|     |     |     |     |       |     |       | ±    |       | ±    |     |     |     |
|     |     | 256 |     | 1/162 | 1/2 | 15.74 | 0.32 | 15.74 | 0.32 |     | 16  |     |
|     |     |     |     |       |     |       | ±    |       | ±    |     |     |     |
InFig.3,Z(t),Y(t),∆ (t)and∆ (t)areshownfor =N =512,ω =2/322,m=1/4.
|     |     |     |     | z   |     | y   |     |     | N z y |     | 0   |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ----- | --- | --- | --- |
e

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007
997
40 40
30
30
20
| 20  |     |     | Z(t) |     | Y(t) |
| --- | --- | --- | ---- | --- | ---- |
|     |     |     | ∆(t) |     | ∆(t) |
10
|     |     |     | z   |     | y   |
| --- | --- | --- | --- | --- | --- |
10
0
0
−10
−10 −20
−20 −30
0 1000 2000 3000 4000 5000 6000 0 1000 2000 3000 4000 5000 6000
|     |     | t   |     | t   |     |
| --- | --- | --- | --- | --- | --- |
Figure 3: Z(t), Y(t), ∆ (t) and ∆ (t) for the harmonic oscillator with parameters Nz =Ny =512, ω =2/322,
|     | z   | y   |     |     | 0   |
| --- | --- | --- | --- | --- | --- |
m=1/4,vz =0.02andvy =0.04. Thesolid linesare Z(t)andY(t),whilethedottedonesare∆ (t)and∆ (t).
z y
| Time andspace | are expressed | in lattice | units. |     |     |
| ------------- | ------------- | ---------- | ------ | --- | --- |
e
8.2 Time-dependentGPE
Letusconsiderthefollowingone-dimensionalnonlinearpotential
1
|     |     | V(z,t)= | mω2(z z )2+V ψ(z,t) | 2,  |     |
| --- | --- | ------- | ------------------- | --- | --- |
0 0 nl
|     |     |     | 2 − | | |   |     |
| --- | --- | --- | ----- | --- | --- |
whereV NU . Thispotentialrepresentsaself-interactingparticleconfinedbyanhar-
| nl  | ≡ 1 |     |     |     |     |
| --- | --- | --- | --- | --- | --- |
monic trap. As an initial condition, we use again the one-dimensional Gaussian packet
of minimum uncertainty. We compare qLB numerical results with the ones given by a
classical Crank-Nicolson (CN) scheme. Forbothschemesthe discretization stepsare set
∆z=∆t=0.5
| ash | andDirichletboundaryconditionsareimposed. |     |     |     |     |
| --- | ----------------------------------------- | --- | --- | --- | --- |
≡
It is well known that if the number of particles exceeds a given threshold (i.e. the
nonlinearitycoefficientislargerthanafixedthreshold),theBECbecomesunstable. This
isduetothefactthat,whenrepulsiveforcesbetweentheparticlesofthecondensatepre-
vail over the confining effect of the harmonic trap, the BEC droplet breaks up. From a
mathematical point of view, this correspondsto thesituation where the potentialdevel-
opsa doubly humpedstructure, so that an initially Gaussian wave packet, representing
theBECdroplet,wouldbreakupintotwoormoreseparatedropletsandfinallydissolve
intoapurelychaoticconfiguration[27].
Inonespatialdimension,thecritical value, V,above whichtheBECbecomesunsta-
c
ble(expressedinlatticeunits)isgivenby
=√2π∆3mω2.
|     |     |     | V c |     | (8.7) |
| --- | --- | --- | --- | --- | ----- |
0 0
Letusconsiderthefollowingsetofparameters: Nez =512,m=1/4,∆ =32andω =1/512.
|     |     |     |     | 0   | 0   |
| --- | --- | --- | --- | --- | --- |
=0.19625.
| FromEq.(8.7),weobtainV |     | c   |     |     |     |
| ---------------------- | --- | --- | --- | --- | --- |
e

998 S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007
Theeffectofthenonlineartermistocyclically contractandexpandthewavepacket.
The amplitude of these contractions/expansions is proportional to the strength of the
nonlinear term (i.e. to the value of V ). However, the condensate remains confined by
nl
theharmonictrap.
>
For V V, after a large number of time steps, this cyclic behavior is broken and a
nl c
transition occurs which corresponds to strong variation of the wave function and leads
toachaoticbehavior,endingwiththetotaldispersionofthewavepacket.
40
38
36
34
32
30
0 0.5 1 1.5 2 2.5 3
t x 104
)t(
∆
z
qLB CN
Figure 4: Mean spreading ∆ z (t) for qLB and CN. Parameters are set as follows: Nz =512, z 0 =256, m=1/4,
∆
0
=32, ω
0
=1/512, vz =0.0, V
nl
=0.1.
e
We initially set V =0.1 < V. In Fig. 4, the mean spreading produced by the two
nl c
schemes up to time t=32200 ( 10T, where T=2π/ω ) is reported. From this figure,
0
∼
wenotethatthepacketspreadingisoscillatingfrom32(whichcorrespondstotheinitial
value) to about 39.2, with an oscillation period which is given by T/2. However, the
twoschemescomputeaslightlydifferentvalueforω . Inparticular, ω =1.926 10 3,
0 qLB −
×
while ω =1.951 10 3, leading to two different periods, T 3260 and T 3220
CN − qLB CN
× ∼ ∼
whiletheexpectedvalueis T 3217. Weconcludethat,qLBisslightlylessaccuratethan
∼
CN in computing the oscillation frequency. This is in line with the general observation
that qLB is very efficient in achieving a reasonable accuracy (within a few percent), but
cannot easily be pushedto high-accuracy because of the lack of adiabaticity in the limit
ω ∆t 0. In Fig. 5, the wave function densities computed by qLB and CN at times 0,
c
→
T /4, T /2, 3T /4, T for qLB and 0, T /4, T /2, 3T /4, T for CN are
qLB qLB qLB qLB CN CN CN CN
shown. We note that the expansion/contraction behavior is well visible and the results
are in good agreement. In Fig. 6, kinetic, potential and total energies computed by qLB
andCNschemesarealsoshown. Asatisfactoryenergyconservationisachievedforboth
models and the mean values of the total energy computed by qLB and CN are E =
qLB
1.413 10 3 andE =1.417 10 3.
− CN −
× ×
>
Let us briefly investigate the transition that occurs for V V after a large number
nl c
oftime steps. In particular, we setV =0.21 > V , leaving unchangedthe remaining pa-
nl c

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007
999
|         |       |         | qLB |     |        |     |               |         | CN  |     |         |
| ------- | ----- | ------- | --- | --- | ------ | --- | ------------- | ------- | --- | --- | ------- |
|         | 0.014 |         |     |     |        |     | 0.014         |         |     |     |         |
|         |       |         |     |     | t=0    |     |               |         |     |     | t=0     |
|         |       |         |     |     | t=T /4 |     |               |         |     |     | t=T /4  |
|         | 0.012 |         |     |     | qLB    |     | 0.012         |         |     |     | CN      |
|         |       |         |     |     | t=T /2 |     |               |         |     |     | t=T /2  |
|         |       |         |     |     | qLB    |     |               |         |     |     | CN      |
|         | 0.01  |         |     |     |        |     | 0.01          |         |     |     |         |
|         |       |         |     |     | t=3T   | /4  |               |         |     |     | t=3T /4 |
|         |       |         |     |     | qLB    |     |               |         |     |     | CN      |
|         |       |         |     |     | t=T    |     |               |         |     |     | t=T     |
| 2|)z(ψ| | 0.008 |         |     |     | qLB    |     | 2|)z(ψ| 0.008 |         |     |     | CN      |
|         | 0.006 |         |     |     |        |     | 0.006         |         |     |     |         |
|         | 0.004 |         |     |     |        |     | 0.004         |         |     |     |         |
|         | 0.002 |         |     |     |        |     | 0.002         |         |     |     |         |
|         | 0     |         |     |     |        |     | 0             |         |     |     |         |
|         | 100   | 150 200 | 250 | 300 | 350    | 400 | 100           | 150 200 | 250 | 300 | 350 400 |
|         |       |         | z   |     |        |     |               |         | t   |     |         |
Figure5: WavefunctiondensitiescomputedbyqLBandCN.Parametersaresetasfollows: Nz =512,z =256,
0
m=1/4, ∆ =32, ω =1/512, vz =0.0, V =0.1. For qLB ψ(z,t) 2 is computed at times 0, T /4, T /2,
|     |     | 0 0 |     |     | nl  |     | | | |     |     |     | qLB qLB |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ------- |
ψ(z,t)2
3T qLB /4 and T qLB , while for CN is computed at times 0, T CN /4, T CN /2, 3T CN /4 and T CN .
|     |     |     |     | |   | |   |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
e
|     |        |       | qLB |       |         |     |        |       | CN  |     |       |
| --- | ------ | ----- | --- | ----- | ------- | --- | ------ | ----- | --- | --- | ----- |
|     | x 10−4 |       |     |       |         |     | x 10−4 |       |     |     |       |
|     | 16     |       |     |       |         |     | 16     |       |     |     |       |
|     | 14     |       |     |       |         |     | 14     |       |     |     |       |
|     | 12     |       |     |       |         |     | 12     |       |     |     |       |
|     | 10     |       |     |       |         |     | 10     |       |     |     |       |
|     | 8      |       |     |       |         |     | 8      |       |     |     |       |
|     |        | E     | E   | E = E |  + E    |     |        | E     | E   | E=E | +E    |
|     |        | kin   | pot |       | kin pot |     |        | kin   | pot | kin | pot   |
|     | 6      |       |     |       |         |     | 6      |       |     |     |       |
|     | 4      |       |     |       |         |     | 4      |       |     |     |       |
|     | 2      |       |     |       |         |     | 2      |       |     |     |       |
|     | 0      | 0.5 1 | 1.5 | 2     | 2.5     | 3   | 0      | 0.5 1 | 1.5 | 2   | 2.5 3 |
|     |        |       | t   |       | x 104   |     |        |       | t   |     | x 104 |
Figure6: Kinetic,potentialandtotalenergycomputedbyqLBandCN.Parametersaresetasfollows: Nz =512,
| z   | =256, m=1/4, | ∆ =32, | ω   | =1/512, | vz =0.0, | V =0.1. |     |     |     |     |     |
| --- | ------------ | ------ | --- | ------- | -------- | ------- | --- | --- | --- | --- | --- |
| 0   |              | 0      |     | 0       |          | nl      |     |     |     |     |     |
e
rameters. Afterabout98000timesteps( 30T qLB ),thetransitionfromthecyclicbehavior
∼
describedabove to a chaotic statetakes place. In Fig. 7, thewave function densitycom-
| putedbyqLBattimes30T |     |     |     | ,32T |     | and50T | isshown. |     |     |     |     |
| -------------------- | --- | --- | --- | ---- | --- | ------ | -------- | --- | --- | --- | --- |
|                      |     |     |     | qLB  | qLB |        | qLB      |     |     |     |     |
A similar experiment can be performed in two spatial dimensions. In this case the
non-linearpotentialisgivenby:
1
|                     |     | V(z,y,t)= |     |     | mω2[(z | z )2+(y | y )2]+V | ψ(z,y,t) |     | 2,  |     |
| ------------------- | --- | --------- | --- | --- | ------ | ------- | ------- | -------- | --- | --- | --- |
|                     |     |           |     |     | 0      | 0       | 0       | nl       |     |     |     |
|                     |     |           |     | 2   |        | −       | −       | |        |     | |   |     |
| whereweareassumingω |     |           |     | ω   | =ω .   |         |         |          |     |     |     |
|                     |     |           |     | 0 ≡ | z y    |         |         |          |     |     |     |
Using,asaninitialconditionthewavepacketofEq.(8.1),onecancalculatethecritical

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007
1000
0.014
t=30 T
qLB
0.012
|     |     | 0.01 | t=50 T |     |     |     |
| --- | --- | ---- | ------ | --- | --- | --- |
qLB
t=33 T
qLB
2 0.008
|)z(ψ|
0.006
0.004
0.002
0
|     |     | 0   | 100 | 200 300 | 400 | 500 |
| --- | --- | --- | --- | ------- | --- | --- |
z
>Vc
Figure 7: Wave function density computed by qLB when the BEC becomes unstable (V nl ). Parameters
∆
are set as follows: Nz =512, z 0 =256, m=1/4, 0 =32, ω 0 =1/512, vz =0.0, V nl =0.21. The wave function
| density is reported | at times | 30T | , 32T and | 50T . |     |     |
| ------------------- | -------- | --- | --------- | ----- | --- | --- |
|                     |          | qLB | qLB       | qLB   |     |     |
e
valueabovewhichtheBECbecomesunstable[27].
Intwodimensions,thisvalueisgivenby:
=2π∆4mω2,
V (8.8)
|     |     |     |     | c 0 0 |     |     |
| --- | --- | --- | --- | ----- | --- | --- |
whereallthequantitiesareexpressedinlatticeunits.
e
Asabove, wecompare qLBwith CNbysettingthesamediscretization steph ∆z=
≡
∆y=∆t=0.5andimposingDirichletboundaryconditionsforbothschemes. Moreover,for
|     |     |     |     | =256,m=1/2,∆ |     | =∆ ∆ |
| --- | --- | --- | --- | ------------ | --- | ---- |
thissimulation,parametersaresetasfollows: N z =N y 0z 0y 0 =16,
≡
| ω =1/256,hence,fromEq.(8.8),weobtainV |     |     |     | =π. |     |     |
| ------------------------------------- | --- | --- | --- | --- | --- | --- |
| 0                                     |     |     |     | c   |     |     |
>
Inthiscase,wesetV =2π V ,andweexpecttoobserveeacyclicbehaviorforalarge
|     |     | nl  | c   |     |     |     |
| --- | --- | --- | --- | --- | --- | --- |
number of time steps and then a rapid transition from this stable, oscillating regime to
theunstabledispersionoftheBEC.
Indeed, for about 100T, the wave packet contracts and expands itself cyclically, as
shown in Fig. 8. In particular, the oscillation frequencies computed by qLB and CN are
| =3.841 | 3   |     | =3.898 | 3   |     |     |
| ------ | --- | --- | ------ | --- | --- | --- |
ω qLB 10 − and ω CN 10 − respectively, leading to the following values
|     | ×   |     | ×   |     |     |     |
| --- | --- | --- | --- | --- | --- | --- |
fortheperiods: T 1635,T 1611,whiletheexpectedvalueisT 1608. InFig.8,the
|     | qLB |     | CN  |     |     |     |
| --- | --- | --- | --- | --- | --- | --- |
|     | ∼   |     | ∼   |     |     | ∼   |
wave function densities taken at y=y and computed by qLB and CN at time intervals
0
ofaquarteroftheirperiodsareshownandtheywitnessasatisfactoryagreement.
Afterabout100T timesteps,arapidtransitionoccursandthewavefunctionstarts
qLB
to break up. In Fig. 9, the unstable behavior of the wave packet is shown, in particular,
the wave function density computed by qLB at the cross section y=y at times 100T
0 qLB
| and110T qLB | isreported. |     |     |     |     |     |
| ----------- | ----------- | --- | --- | --- | --- | --- |

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007
1001
|          |     | qLB |     |      |     |        |     | CN  |     |         |
| -------- | --- | --- | --- | ---- | --- | ------ | --- | --- | --- | ------- |
| x 10−4   |     |     |     |      |     | x 10−4 |     |     |     |         |
| 8        |     |     |     |      |     | 8      |     |     |     |         |
|          |     |     |     | t=0  |     |        |     |     |     | t=0     |
| 7        |     |     |     | t=T  | /4  | 7      |     |     |     | t=T /4  |
|          |     |     |     |      | qLB |        |     |     |     | CN      |
|          |     |     |     | t=T  | /2  |        |     |     |     | t=T /2  |
| 6        |     |     |     |      | qLB | 6      |     |     |     | CN      |
|          |     |     |     | t=3T | /4  |        |     |     |     | t=3T /4 |
| 2|)) 5   |     |     |     |      | qLB | 5      |     |     |     | CN      |
|          |     |     |     | t=T  |     | 2|)    |     |     |     | t=T     |
|          |     |     |     |      | qLB | 0      |     |     |     | CN      |
| y,z(ψ| 0 |     |     |     |      |     | y,z(ψ| |     |     |     |         |
| 4        |     |     |     |      |     | 4      |     |     |     |         |
| 3        |     |     |     |      |     | 3      |     |     |     |         |
| 2        |     |     |     |      |     | 2      |     |     |     |         |
| 1        |     |     |     |      |     | 1      |     |     |     |         |
| 0        |     |     |     |      |     | 0      |     |     |     |         |
| 50       | 100 |     | 150 |      | 200 | 50     | 100 |     | 150 | 200     |
|          |     |     | z   |      |     |        |     | z   |     |         |
Figure 8: Wavefunction densities computedbyqLBandCN atthecross section y=y . Parameters areset as
0
follows: Nz =Ny =256,z =y =128,m=1/2,∆ =16,ω =1/256,vz =vy =0.0,V =2π. ForqLB ψ(z,y,t) 2
|     |     | 0   | 0   |     | 0   | 0   |     | nl  |     | | | |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
ψ(z,y,t)2
is computed at times 0, T qLB /4, T qLB /2, 3T qLB /4 and T qLB , while for CN is computed at times 0,
|         |        |        |     |     |     |     | |   | |   |     |     |
| ------- | ------ | ------ | --- | --- | --- | --- | --- | --- | --- | --- |
| T /4, T | /2, 3T | /4 and | T . | e   |     |     |     |     |     |     |
| CN CN   | CN     |        | CN  |     |     |     |     |     |     |     |
−4
x 10
8
t=100 T
qLB
7
t=110 T
qLB
6
2 5
|)
0
y,z(ψ| 4
3
2
1
0
|     |     |     | 0   | 50  | 100 | 150 | 200 | 250 |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
z
Figure 9: Wave function density computed by qLB when the BEC becomes unstable (V >Vc ). Parameters
nl
areset asfollows: NzNy =256, z =y =128, m=1/2, ∆ =16, ω =1/256, vz =vy =0.0,V =2π. Thewave
|                    |          |     | 0        | 0        |      | 0   | 0   |     | nl  |     |
| ------------------ | -------- | --- | -------- | -------- | ---- | --- | --- | --- | --- | --- |
| function densityis | reported |     | at times |          | and  | .   |     |     |     |     |
|                    |          |     |          | 100T qLB | 110T | qLB |     |     |     |     |
e
| 8.3 Ground-state |     | computation |     |     | of theGPE |     |     |     |     |     |
| ---------------- | --- | ----------- | --- | --- | --------- | --- | --- | --- | --- | --- |
Inthissection,weshowsomenumericalresultsobtainedbyapplyingtheimaginary-time
qLB model to the computation of the ground state of the GPE. Since, qualitatively, one-
and two- dimensional results are very similar, we will focus only on two-dimensional

1002 S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007
Table 3: Ground state chemical potential µ and maximum value reached by the ground state wave function
φ0
g≡
φg (x
0
,y
0
) for qLB, CN and BEFD models. Numerical results are also compared with the Thomas-Fermi
chemicalpotential(seeEq.(8.10)). TheresultsarecomputedfordifferentvaluesofV . Parametersaresetas
nl
follows: ωz =ωy =1/128, m=1/8, ∆
0z
=∆
0y
=16, Nz =Ny =512.
V µqLB µCN µBEFD µTF φ0 qLB φ0 CN φ0 BEFD
nl g g g
e
0 0.007816 0.007812 0.007812 – 0.01723 0.01763 0.01763
10 0.009219 0.009250 0.009250 0.004928 0.01627 0.01656 0.01656
100 0.017489 0.017597 0.017597 0.015584 0.01218 0.01226 0.01225
500 0.035802 0.035949 0.035949 0.034846 0.00835 0.00837 0.00836
1000 0.049964 0.050125 0.050125 0.049280 0.00702 0.00704 0.00703
2000 0.070161 0.070338 0.070338 0.069692 0.00590 0.00591 0.00591
3000 0.085721 0.085905 0.085905 0.085355 0.00533 0.00534 0.00534
4000 0.098860 0.099050 0.099050 0.098560 0.00496 0.00497 0.00496
5000 0.110447 0.110642 0.110642 0.110193 0.00469 0.00470 0.00470
10000 0.155967 0.156176 0.156176 0.155837 0.00395 0.00395 0.00395
simulations. Fordetailsonone-dimensionalnumericalresults,wereferthereadertothe
originalwork[21].
In order to validate our numerical results, we compare qLB with classical Crank-
NicholsonschemeandwithabackwardEulerfinitedifference(BEFD)scheme[21,33]. We
also use, as a second termof comparison, the analytic solution obtained in the Thomas-
Fermilimit. This approximation is valid in the strong-interaction limit, in which kinetic
energy contributions can be neglected [41]. This limit is reached for large values of the
non-linearitycouplingconstant NU .
d
By considering the time-independentGPE of Eq. (7.1) and neglecting the kinetic en-
ergyterm,weobtain
µ φ(r)=(V (r)+NU φ(r) 2)φ(r),
TF ext d
| |
whereweindicateµwithµ torecallthatthisistheThomas-Fermichemicalpotential.
TF
Inthiscasethesolutionistrivialandthewavefunctiondensityisgivenby
1
φ(r) 2= (µ V (r))Θ(µ V (r)), (8.9)
TF ext TF ext
| | NU − −
d
whereΘistheHeavisidestepfunction. Thechemicalpotentialgivenbythisapproxima-
tion,µ ,canbecomputedbydirectlyimposingthenormalizationconditiontothewave
TF
functiondensityofEq.(8.9). Inthetwo-dimensionalcase,weobtain
mω2 1/2
µ = NU 0 . (8.10)
TF 2
π
(cid:18) (cid:19)
Letusconsiderthefollowingtwo-dimensionalnonlinearpotential
1
V(z,y,τ)= mω2[(z z )2+(y y )2]+V φ(z,t,τ) 2,
2 0 − 0 − 0 nl | |

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007
1003
|     |        |     | qLB-CN |     |     |        |     |        |     | qLB-BEFD |     |          |     |
| --- | ------ | --- | ------ | --- | --- | ------ | --- | ------ | --- | -------- | --- | -------- | --- |
|     | x 10−4 |     |        |     |     |        |     | x 10−4 |     |          |     |          |     |
|     |        |     |        |     | V   | =10    |     |        |     |          |     | V =10    |     |
|     |        |     |        |     |     | nl     |     |        |     |          |     | nl       |     |
|     |        |     |        |     | V   | =100   |     | |)     |     |          |     | V =100   |     |
|     | |)     |     |        |     |     | nl     |     | 0      |     |          |     | nl       |     |
|     | 0      |     |        |     | V   | =1000  |     | y,z(   |     |          |     | V =1000  |     |
|     | y,z(   |     |        |     |     | nl     |     |        |     |          |     | nl       |     |
|     | 2      |     |        |     | V   | =10000 |     | 2      |     |          |     | V =10000 |     |
|     | NC     |     |        |     |     | nl     |     | DFEB   |     |          |     | nl       |     |
φ−)
φ−)
0
|     | y,z( |     |     |     |     |     |     | 0   |     |     |     |     |     |
| --- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
y,z(
|     | 1   |     |     |     |     |     |     | 1   |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
BLq
BLq
φ|
φ|
|     | 0   |     |     |     |     |     |     | 0   |     |     |     |         |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ------- | --- |
|     | 0   | 100 | 200 | 300 | 400 | 500 |     | 0   | 100 | 200 | 300 | 400 500 |     |
|     |     |     |     | z   |     |     |     |     |     |     | z   |         |     |
Figure 10: Deviations of qLB from CN and BEFD in the computation of the ground state wave function
for different values of V . In particular, (φg ) (z,y ) (φg ) (z,y ) and (φg ) (z,y ) (φg ) (z,y )
|     |     |     | nl  |     |     | qLB | 0   | CN  | 0   |           | qLB | 0 BEFD | 0   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --------- | --- | ------ | --- |
|     |     |     |     |     | |   |     | −   |     | |   | m=1/8,∆ | | =∆  | −      | |   |
computedattheqLBnodalareplotted. Simulationparametersaresetas: 0z 0y =16,ω 0 =1/128,
| Nz  | =Ny =512. | Spaceis | expressed |     | in lattice | units. |     |     |     |     |     |     |     |
| --- | --------- | ------- | --------- | --- | ---------- | ------ | --- | --- | --- | --- | --- | --- | --- |
e
0.02
0.015
|     |     |     | )   | 0   |     |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
y,z(
0.01
g
φ
0.005
0
|     |     |     |     |     | 0   | 100 | 200 | 300 | 400 | 500 |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
z
Figure 11: Ground stateprofiles given bytheqLBmodel for different values of V . Simulation parametersare
nl
set as: m=1/8, ∆ =∆ =16, ω =1/128, Nz =Ny =512. V takes the following values: 0, 10, 100, 500,
|       |      |           | 0z 0y |             | 0       |           |     |            | nl     |     |     |     |     |
| ----- | ---- | --------- | ----- | ----------- | ------- | --------- | --- | ---------- | ------ | --- | --- | --- | --- |
| 1000, | 5000 | and 10000 | (top  | to bottom). | Spaceis | expressed |     | in lattice | units. |     |     |     |     |
e
where V =NU and τ is the fictitious time variable introduced in Section 7. As initial
|     | nl  |     | 2   |     |     |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
condition, the wave packet of Eq. (8.1) is again used and Dirichlet boundary conditions
areimposedforallthethreeschemes.
Let [0,512] [0,512] be our numerical domain (expressed in lattice units) and the re-
×
maining parameters are set as follows: ∆ =∆ =16, ω =1/128 and m=1/8. Dis-
|     |     |     |     |     |     |     | 0z  | 0y  |     | 0   |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
cretizationstepsforqLBaresettounity,whileforCNandBEFDweset∆ =∆y=0.5and
z
∆t=0.1.
Thesimulation mustbe stoppedonly whenastationarystateis reeached,hence

1004 S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007
x 10−3
8
7
6
5
4
3
2
1
0
0 100 200 300 400 500
z
)
y,z(
φ
0
g
x 10−3
5
V =1000 qLB
nl TF
4
3
2
1
0
0 100 200 300 400 500
z
)
y,z(
φ
0
g
V =5000 qLB
nl TF
x 10−3
4
3.5
3
2.5
2
1.5
1
0.5
0
0 100 200 300 400 500
z
)
y,z(
φ
0
g
V =10000 qLB
nl TF
Figure12: Groundstateprofileφg (z,y)takenaty=y 0
for different values of V and compared with the nl Thomas-Fermi wave function (see Eq. (8.9)). Simu-
lationparametersaresetas: m=1/8,∆ =∆ =16,
0z 0y
ω
0
=1/128, Nz =Ny =512. Solid lines: qLB model;
dashedlines: Thomas–Fermiaepproximation. Spaceis
expressed in lattice units. The tails associated with
thekineticenergy contribution arewell visible.
thefollowingstopcriterionisimposed
max φn+1 φn < ε,
i,j=0,
···
,Ng | i,j − i,j|
where N isthenumberofnodalpointandε=10 9.
g −
InTable3,thelimit value ofµandthemaximumvalueofφ attheendofthesimula-
tion,(φ (z ,y )),arereportedforthethreeschemes. Moreover,theThomas-Fermichem-
g 0 0
ical potential µ given by Eq. (8.10) is also shown. In Fig. 10, we compare the ground
TF
state wave function φ (z,y) taken at y=y computed by qLB with the profiles given by
g 0
CNandBEFD.Inparticular,inthefigurewereportthedeviationsbetweenqLBandCN
orBEFD: (φ ) (z,y ) (φ ) (z,y ) and (φ ) (z,y ) (φ ) (z,y ) computedat
g qLB 0 g CN 0 g qLB 0 g BEFD 0
| − | | − |
theqLBnodalpoints. InFig.11, φ (z,y) takenat y=y andcomputedbyqLBfordiffer-
g 0
ent values of V is reported to show the qualitative effect of an increasing nonlinearity
nl
coefficient. Finally, in Fig. 12, we compare qLB results with the wave function given by
theThomas-FermiapproximationofEq.(8.9)forsomevaluesofV . Thetailsmoothing,
nl
well visible in the qLB results, is due to the kinetic energy contribution, hence, it is not
reproducedbytheThomas-Fermiapproximationwherethekinetictermisneglected.
Thesedata witness a satisfactory agreement between qLB and the reference CN and

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007 1005
BEFDsolutions, while CN and BEFDare in excellent agreement with each other (this is
duetothehigherresolutionadoptedinthesereferencecases).
9 Conclusions and outlook
In this work, we have reviewed the derivation ofthe quantum lattice Boltzmann model
and itsmostrecentdevelopments. Inparticular, wehave showntheviability oftheqLB
schemeforthenumericalsolutionofthetime-dependentSchro¨dingerequation,aswellas
itsnon-linearextension,suchastheGross-PitaevskiiequationforBose-Einsteinconden-
sates, in multiple spatial dimensions. Moreover, the formulation of an imaginary-time
qLBmodelfortheground-statecomputationoftheGPEhasalsobeenreviewed.
Being based on a unitary, first-order, relativistic formulation, at variance with most
explicitschemesfornon-relativisticquantumwaveequations,theqLBmethodoffersun-
conditional stability with the size of the time-step/mesh-size. In addition, being based
onafirst-order,hyperbolicformulation, stability can bepreservedwithatime-stepscal-
ing linearly with the mesh size, rather than quadratically like most explicit schemes for
quantumwaveequations[34].
However, its accuracy is subject to the condition ω ∆t=∆x/λ 1, λ =c/ω being
c B B c
≤
the De Broglie wavelength of the particle. Since the time-step scales linearly with the
mesh-spacing (a result of the relativistic formulation), qLB can be taken down to very
refinedgridswithoutsufferingthetime-stepcollapsetypicalofnon-relativisticCourant-
Friedrichs-Lewystabilityconditions,∆t < 2m∆x2/h¯,thuscompensatingforitslow-order
accuracy. However, care must be taken to ensure that errors due to lack of adiabaticity
remainundercontrolwhenω ∆tissenttozero.
c
TheqLBmethodisalsoveryinterestingasaprospectivealgorithmforquantumcom-
puters [8,11,13–15]. Indeed, as observed in [42], the stream-and-collide structure of the
quantumlattice Boltzmann equation maps naturally ontothe structureof quantumnet-
works, i.e. quantum computing devices consisting of quantum logic gates, whose com-
putational operation proceeds synchronously in time. In these networks, the output of
some gates are wire-connected to the input of some others(streaming step), and locally
processedbyunitaryoperations(thecollisionstep).
Besidesthis attractive, but still speculative application, qLB is an interestingscheme
that can be easily implemented in classical (electronic) computers. As an explicit nu-
merical scheme,qLB offers an appealing setoffeatures,such as unconditionedstability,
norm-preserving(unitarity)andamenabilitytoparallelprocessing.
In conclusion, the qLB method makes an excellent candidate for implementation on
classicalcomputersaswellasforprospectivequantumcomputingapplications.
References
[1] G. McNamara and G. Zanetti, Use of the Boltzmann equation to simulate lattice-gas au-
tomata,Phys.Rev.Lett.,61(1988),2332–2335.

1006 S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007
[2] F.HigueraandJ.Jimenez,Boltzmannapproachtolatticegassimulations,Europhys.Lett.,9
(1989),663–668.
[3] F.Higuera,S.SucciandR.Benzi,Latticegasdynamicswithenhancedcollisions,Europhys.
Lett.,9(1989),345–349.
[4] R. Benzi, S. Succi and M. Vergassola, The lattice Boltzmann equation: theory and applica-
tions,Phys.Rep.,222(1992),145–197.
[5] J.L.Lebowitz,S.OrszagandY.Qian(eds.),J.Stat.Phys.,SpecialissueonLatticemodels,81
(1995).
[6] I. V. Karlin, S. Ansumali, C. E. Frouzakis and S. S. Chikatamarla, Elements of the lattice
Boltzmann method I: Linear advection equation, Commun. Comput. Phys., 1 (2006), 616-
655.
[7] I.V.Karlin,S.S.ChikatamarlaandS.Ansumali,ElementsofthelatticeBoltzmannmethodII:
Kineticsandhydrodynamicsinonedimension.Commun.Comput.Phys.,2(2007),196-238.
[8] B. M. Boghosian and W. Taylor, Simulating quantum mechanics on a quantum computer,
PhysicaD,120(1998),30–42.
[9] B. M. Boghosian and W. Taylor, Quantum lattice-gas model for the many-particle
Schro¨dingerequationinddimensions,PhysicalReviewE,57(1998),54–66.
[10] I.Bialynicki-Birula, Weyl, Dirac, andMaxwellequations onalatticeasunitarycellularau-
tomata,PhysicalReviewD,49(1994),6920–6927.
[11] D. A. Meyer, From quantum cellular automata to quantum lattice gases, J. Stat. Phys., 85
(1996),551–574.
[12] D.A.Meyer,Quantumlatticegasesandtheirinvariants,Int.J.Mod.Phys.C,8(1997),717–
735.
[13] J. Yepez, Quantum computation for physical modeling, Comput. Phys. Commun., 146
(2002),277–279.
[14] J. Yepez, Quantum lattice-gas model for computational fluid dynamics, Phys. Rev. E, 63
(2001),046702.
[15] G.Vahala,L.VahalaandJ.Yepez,Quantumlatticegasrepresentationofsomeclassicalsoli-
tons,Phys.Lett.A,310(2003),187–196.
[16] S. Succi and R. Benzi, Lattice Boltzmann equation for quantum mechanics, Physica D, 69
(1993),327–332.
[17] R.Feynman,SimulatingPhysicswithcomputers,Int.J.Theoret.Phys.,21(1982),467–488.
[18] D.P.diVincenzo,Quantumcomputation,Science,270(1995),255–261.
[19] S.Succi,NumericalsolutionoftheSchro¨dingerequationusingdiscretekinetictheory,Phys.
Rev.E,53(1996),1969–1975.
[20] S.PalpacelliandS.Succi,NumericalvalidationofthequantumlatticeBoltzmannschemein
twoandthreedimensions,Phys.Rev.E.,75(2007),066704.
[21] S.Palpacelli,S.SucciandR.Spigler,Ground-statecomputationofBose-Einsteincondensates
byanimaginary-timequantumlatticeBoltzmannscheme,Phys.Rev.E,76(2007),036712.
[22] S. Succi, Lattice Boltzmann schemes for quantum applications, Comp. Phys. Comm., 146
(2002),317–323.
[23] E.Madelung,Quantumtheoryinhydrodynamicalform,Z.Phys.,40(1926),332–336.
[24] L.Landau,E.Lifshitz,RelativisticQuantumFieldTheory,(Pergamon,Oxford,1960).
[25] S.Succi,LatticeBoltzmannequationforrelativisticquantummechanics,Phil.Trans.R.Soc.
Lond.A,360(2002),429–436.
[26] L.I.Schiff,QuantumMechanics,3rdeditionMcGraw-Hill,NewYork(1968).
[27] S. Succi, Lattice quantum mechanics: an application to Bose-Einstein condensation, Int. J.

S.PalpacelliandS.Succi/Commun.Comput.Phys.,4(2008),pp.980-1007 1007
Mod.Phys.C,9(1998),1577–1585.
[28] A. Griffin, D. Snoke and S. Stringari (eds.),Bose-Einstein Condensation (Cambridge Univ.
Press,NewYork,1995).
[29] A.D.Jackson,G.M.KavoulakisandC.J.Pethick,SolitarywavesincloudsofBose-Einstein
condensedatoms,Phys.Rev.A,58(1998),2417–2422.
[30] P. Leboeuf and N. Pavloff, Bose-Einstein beams: Coherent propagation through a guide,
Phys.Rev.A,64(2001),033602.
[31] W. Bao and W. Tang, Ground state solution of Bose-Einstein condensate by directly mini-
mizingtheenergyfunctional,J.Comput.Phys.,187(2003),230–254.
[32] S.K. Adhikari, Numericalstudyof the sphericallysymmetric Gross-Pitaevskiiequation in
twospacedimensions,Phys.Rev.E,62(2000),2937–2944.
[33] W.BaoandQ.Du,ComputingthegroundstatesolutionofBose-Einsteincondensatesbya
normalizedgradientflow,SIAMJ.Sci.Comput.,25(2004),1674–1697.
[34] K. C. Kulander (ed.), Time-Dependent Methods for Quantum Dynamics, thematic issue,
Comp.Phys.Commun.63,(North-Holland,Amsterdam,1991).
[35] W.Bao, ChernandF.Y.Lim, Efficientandspectrallyaccuratenumericalmethodsforcom-
puting ground and first excited states in Bose-Einstein condensates, J. Comput. Phys., 219
(2006),836–854.
[36] A. Aftalion and Q. Du, Vortices in a rotating Bose-Einstein condensate: Critical angular
velocitiesandenergydiagramsintheThomas-Fermiregime,Phys.Rev.A,64(2001),063603.
[37] M.M.Cerimele,M.L.Chiofalo,F.Pistella,S.SucciandM.P.Tosi,Numericalsolutionofthe
Gross-Pitaevskiiequationusinganexplicitfinite-differencescheme,Phys.Rev.E,62(2000),
1382–1389.
[38] M. L. Chiofalo, S. Succi and M. P. Tosi, Ground state of trapped interacting Bose-Einstein
condensatesbyanexplicitimaginary-timealgorithm,Phys.Rev.E,62(2000),7438-7444.
[39] A. Minguzzi, S. Succi, F. Toschi, M. P. Tosi and P. Vignolo, Numerical methods for atomic
quantum gases with applications to Bose-Einstein condensates and to ultracold fermions,
Phys.Rep.,395(2004),223–355.
[40] S.Gasiorowicz,QuantumPhysics(Wiley,NewYork,1974).
[41] F.Dalfovo,L.PitaevskiiandS.Stringari,Thecondensatewavefunctionofatrappedatomic
gas,J.Res.Natl.Inst.Stand.Technol.,101(1996),537–544.
[42] P.LoveandB.Boghosian,TypeIIquantumalgorithms,PhysicaA,362(2006),210–214.

