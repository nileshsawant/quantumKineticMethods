PHYSICALREVIEWE83,046706(2011)
Isotropyofthree-dimensionalquantumlatticeBoltzmannschemes
P.J.Dellar*andD.Lapitski
OCIAM,MathematicalInstitute,24–29StGiles’,OxfordOX13LB,UnitedKingdom
S.Palpacelli
NumidiaSocieta`aresponsabilita`limitata,viaBerna31,I-00144Roma,Italy
S.Succi
IstitutoApplicazionidelCalcolo,ConsiglioNazionaledelleRicerche,viadeiTaurini19,I-00185Roma,Italy
(Received8December2010;published7April2011)
NumericalsimulationswithpreviousformulationsofthequantumlatticeBoltzmann(QLB)schemeinthree
spatialdimensionsshowedsignificantlackofisotropy.IntwoormorespatialdimensionstheQLBapproachrelies
uponoperatorsplittingtodecomposethetimeevolutionintoasequenceofapplicationsoftheone-dimensional
QLBschemealongcoordinateaxes.Eachapplicationmustbeaccompaniedbyarotationofthewavefunction
intoabasisofchiraleigenstatesalignedalongtherelevantaxis.Thepreviouslyobservedlackofisotropywasdue
toaninconsistencyintheapplicationoftheserotations.Oncethisinconsistencyisremoved,theQLBschemeis
showntoexhibitisotropicbehaviortowithinanumericalerrorthatscalesapproximatelylinearlywiththelattice
spacing.ThisestablishestheviabilityoftheQLBapproachintwoandthreespatialdimensions.
DOI:10.1103/PhysRevE.83.046706 PACSnumber(s): 02.70.−c,03.65.−w
I. INTRODUCTION discrete Boltzmann equation, Eq. (8) below, is very unusual
|     |     |     |     |     |     |     |     | in having | characteristic |     | curves. | The | most one | may | normally |
| --- | --- | --- | --- | --- | --- | --- | --- | --------- | -------------- | --- | ------- | --- | -------- | --- | -------- |
ThelatticeBoltzmannapproachhasdevelopedintoavery
|            |               |     |             |                |     |                 |          | expect is             | for the | characteristic |           | surfaces,  | the | surfaces | along     |
| ---------- | ------------- | --- | ----------- | -------------- | --- | --------------- | -------- | --------------------- | ------- | -------------- | --------- | ---------- | --- | -------- | --------- |
| attractive | computational |     | tool        | for simulating |     | fluid dynamics, |          |                       |         |                |           |            |     |          |           |
|            |               |     |             |                |     |                 |          | which discontinuities |         | in             | solutions | propagate, |     | to       | be three- |
| especially | in complex    |     | geometries, | and            | on  | modern          | parallel |                       |         |                |           |            |     |          |           |
dimensionalsurfaceswithinthefour-dimensional(x,t)space
| computer | architectures     |     | [1–5].    | These | advantages    | follow | from   |               |           |     |           |     |               |     |          |
| -------- | ----------------- | --- | --------- | ----- | ------------- | ------ | ------ | ------------- | --------- | --- | --------- | --- | ------------- | --- | -------- |
|          |                   |     |           |       |               |        |        | [13]. In more | algebraic |     | language, | the | flux Jacobian |     | matrices |
| solving  | the Navier-Stokes |     | equations |       | not directly, | but    | recov- |               |           |     |           |     |               |     |          |
fortheDiracequationarerealandsymmetricwhentheDirac
| ering them | instead | as  | the equations |     | that | describe | slowly |          |            |             |     |       |       |     |            |
| ---------- | ------- | --- | ------------- | --- | ---- | -------- | ------ | -------- | ---------- | ----------- | --- | ----- | ----- | --- | ---------- |
|            |         |     |               |     |      |          |        | equation | is written | in Majorana |     | form, | as in | Eq. | (6) below. |
varyingsolutionsofadiscreteBoltzmannequation.Thelatter
However,thesematricesarenotdiagonal,andtheycannotbe
| is a linear, | constant |     | coefficient | hyperbolic |     | system | that is |     |     |     |     |     |     |     |     |
| ------------ | -------- | --- | ----------- | ---------- | --- | ------ | ------- | --- | --- | --- | --- | --- | --- | --- | --- |
simultaneouslydiagonalizedbyanychangeofbasis.Therefore
| readily | discretized | by  | integrating | along | characteristics. |     | All |     |     |     |     |     |     |     |     |
| ------- | ----------- | --- | ----------- | ----- | ---------------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
theconstructionofnumericalsolutionstotheDiracequation
| nonlinearity | is  | confined | to algebraic |     | terms | that are | local to |           |         |                |     |       |                |     |         |
| ------------ | --- | -------- | ------------ | --- | ----- | -------- | -------- | --------- | ------- | -------------- | --- | ----- | -------------- | --- | ------- |
|              |     |          |              |     |       |          |          | cannot be | reduced | to integrating |     | along | characteristic |     | curves. |
eachgridpoint,andmodelcollisionsbetweenparticles.
TheultimatesourceofthisdifficultyliesintheDiracequation’s
IndevelopingquantumlatticeBoltzmann(QLB)methods
|     |     |     |     |     |     |     |     | description | of  | particles | with | spin, | and the | noncommuting |     |
| --- | --- | --- | --- | --- | --- | --- | --- | ----------- | --- | --------- | ---- | ----- | ------- | ------------ | --- |
[6–9]weaimtoexploittheseadvantagestodevelopalgorithms
propertiesofthePaulispinmatrices[10].
forsimulatingquantumsystemsdescribedbytheSchro¨dinger
|     |     |     |     |     |     |     |     | A number |     | of quantum | cellular |     | automata | (QCA) | have |
| --- | --- | --- | --- | --- | --- | --- | --- | -------- | --- | ---------- | -------- | --- | -------- | ----- | ---- |
equation.Again,wedonotattempttosimulatetheSchro¨dinger
beenformulated[14–19],followingFeynman’scheckerboard
equationdirectly,butinsteadrecoveritasthelimitingequation
|     |     |     |     |     |     |     |     | representation |     | of solutions | to  | the | one-dimensional |     | Dirac |
| --- | --- | --- | --- | --- | --- | --- | --- | -------------- | --- | ------------ | --- | --- | --------------- | --- | ----- |
thatdescribesslowlyvaryingsolutionsoftheDiracequation.
|                  |     |                           |     |     |     |             |     | equation | in terms | of particles |     | executing | random | walks | in a |
| ---------------- | --- | ------------------------- | --- | --- | --- | ----------- | --- | -------- | -------- | ------------ | --- | --------- | ------ | ----- | ---- |
| TheDiracequation |     | offersaquantum-mechanical |     |     |     | description |     |          |          |              |     |           |        |       |      |
of a particle that is compatible with special relativity [10]. discrete space-time [20]. In particular, quantum lattice gas
|     |     |     |     |     |     |     |     | cellular automata |     | (QLGCA) | [15–19] |     | are a | specific | class of |
| --- | --- | --- | --- | --- | --- | --- | --- | ----------------- | --- | ------- | ------- | --- | ----- | -------- | -------- |
LikethediscreteBoltzmannequation,theDiracequationisa
|                  |     |             |            |     |         |      |           | QCA that | adopt | the two-step |     | “stream-collide” |     | rule | typical |
| ---------------- | --- | ----------- | ---------- | --- | ------- | ---- | --------- | -------- | ----- | ------------ | --- | ---------------- | --- | ---- | ------- |
| linear, constant |     | coefficient | hyperbolic |     | system, | with | algebraic |          |       |              |     |                  |     |      |         |
source terms. Both equations contain only single partial oflatticegascellularautomatafromhydrodynamics[21,22].
|               |      |                  |         |           |        |              |     | All these     | quantum | automata | deal        | with | continuous |     | complex   |
| ------------- | ---- | ---------------- | ------- | --------- | ------ | ------------ | --- | ------------- | ------- | -------- | ----------- | ---- | ---------- | --- | --------- |
| derivatives   | with | respect          | to time | and       | space. | By contrast, | the |               |         |          |             |      |            |     |           |
|               |      |                  |         |           |        |              |     | scalar fields | instead | of       | the Boolean |      | variables  | of  | classical |
| Navier–Stokes |      | and Schro¨dinger |         | equations |        | both contain | one |               |         |          |             |      |            |     |           |
derivativewithrespecttotime,buttwoderivativeswithrespect automata that take just the values 0 or 1. Instead, QLGCA
tospace.Thisasymmetrybetweenspaceandtimerendersthe propagate quantum bits (qbits) that correspond to arbitrary
|     |     |     |     |     |     |     |     | superpositions |     | of the two | basic | quantum |     | states | |0(cid:2) and |
| --- | --- | --- | --- | --- | --- | --- | --- | -------------- | --- | ---------- | ----- | ------- | --- | ------ | ------------- |
Schro¨dingerequationincompatiblewithspecialrelativity,and
|     |     |     |     |     |     |     |     | |1(cid:2). |     |     |     |     |     |     | =   |
| --- | --- | --- | --- | --- | --- | --- | --- | ---------- | --- | --- | --- | --- | --- | --- | --- |
alsoincompatiblewithbeingdiscretizedbyintegrationalong More precisely, a quantum bit may be written as q
|                  |     |     |     |     |     |     |     | |0(cid:2)+c  | |1(cid:2),wherec |        |      |                            |     |     |     |
| ---------------- | --- | --- | --- | --- | --- | --- | --- | ------------ | ---------------- | ------ | ---- | -------------------------- | --- | --- | --- |
| characteristics. |     |     |     |     |     |     |     | c            |                  | andc   |      | arecomplexscalars,suchthat |     |     |     |
|                  |     |     |     |     |     |     |     | 0            | 1                | 0      | 1    |                            |     |     |     |
|                  |     |     |     |     |     |     |     | p =|c |2andp |                  | =|c |2 | =1−p | givetheprobabilityforthe   |     |     |     |
ThediscreteBoltzmannandDiracequationsarebothlinear, 0 0 1 1 0
quantumparticletobeinstate|0(cid:2)(empty)or|1(cid:2)(occupied),
symmetrichyperbolicsystemsinthesenseofFriedrichs[11].
Thesepropertiesguaranteetheexistenceofweaksolutionsto respectively.Qubitsexhibitaone-to-onecorrespondencewith
|               |          |     |                |     |                       |     |     | spin-1/2 | fermions | whose | up-down |     | components |     | propagate |
| ------------- | -------- | --- | -------------- | --- | --------------------- | --- | --- | -------- | -------- | ----- | ------- | --- | ---------- | --- | --------- |
| initial value | problems |     | [12]. However, |     | the three-dimensional |     |     |          |          |       |         |     |            |     |           |
±z,
|     |     |     |     |     |     |     |     | along opposite |     | directions, | say |     | in physical |     | space. As |
| --- | --- | --- | --- | --- | --- | --- | --- | -------------- | --- | ----------- | --- | --- | ----------- | --- | --------- |
QLGCAsupportarbitrarysuperpositionsoftheBooleanstates,
*dellar@maths.ox.ac.uk they are formally equivalent to QLB schemes [16,18]. For
1539-3755/2011/83(4)/046706(9) 046706-1 ©2011AmericanPhysicalSociety

P.J.DELLAR,D.LAPITSKI,S.PALPACELLI,ANDS.SUCCI PHYSICALREVIEWE83,046706(2011)
example,aQCAstudiedbyMeyerisunitarilyequivalenttothe where the off-diagonal blocks are the three Pauli spin
one-dimensional QLB scheme [15]. In this respect, QLGCA matrices[10]
(cid:8) (cid:9) (cid:8) (cid:9) (cid:8) (cid:9)
and QLB can both be regarded as computational “engines” 0 1 0 −i 1 0
thatpropagateandcollidequbitsincompliancewiththerule σx= , σy= , σz= . (4)
1 0 i 0 0 −1
of relativistic quantum mechanics. By contrast, a classical
lattice Boltzmann scheme is obtained by coarse graining the Thesematricessatisfythewell-knowncommutationrelations
microscopicBooleanconfigurationsinaclassicalLGCA,via
σiσj −σjσi =2i(cid:8) σk, (5)
eitherensembleorspace-timeaveraging.Asimilarbottom-up ijk
relation has long been advocated in the quantum context by where(cid:8) isthealternatingLevi-Civitatensor.Theessential
ijk
’tHooft[23,24],whohasproposedthatquantum-mechanical fact that underlies this paper is that the spin matrices do not
wavedynamicsmayemergeasthecontinuumlimitofanun- commute. There is thus no basis in which all three matrices
derlyingdeterministiclocalandreversiblecellularautomaton. aresimultaneouslydiagonal,andthesameholdsforthe4×4
Besides reviving the time-honored “hidden variable” picture matricesαx,αy,αz builtfromthespinmatrices.
in a modern guise, ’t Hooft’s approach offers tantalizing but Forourpurposesitisconvenienttotransformthestandard
stillhighlyspeculativeconnectionswithquantumgravity. formoftheDiracequationintotheMajoranaform[10],
Turning to more immediate applications, nonlinear exten-
[∂ +c(−αx∂ +β∂ −αz∂ )+iω αy −igI]ψ =0, (6)
sionsoftheQLBschememaybeusedtosimulatetheGross- t x y z c
Pitaevskii equation that describes Bose-Einstein condensates bymultiplyingthestandardformontheleftandontheright
[9,25]. The underlying capability of the QLB schemes for bytheinvolutionmatrixU =2−1/2(αy +β).Thetransformed
solving the Dirac equation [26] suggests further applications wave function in Eq. (6) is ψ =U(cid:3). The Majorana trans-
tosimulatingBose-Einsteincondensatesinhoneycomboptical formation interchanges the roles of the αy and β matrices
lattices [27]. The motion of charge carriers in graphene, the between Eqs. (1) and (6). The resulting equation (6) has the
two-dimensional hexagonal lattice form of carbon, is also desirablepropertythatthethreematricesαx,β,αzmultiplying
described by the Dirac equation [28–31], with the Fermi thespatialderivatives∂ ,∂ ,∂ areallreal.Thematrixαywith
x y z
velocityv ≈106ms−1 actingastheeffectivespeedoflight. complexcoefficientshasbeenmovedintothealgebraicterms.
F
Weexpectthatthetwo-dimensionalQLBalgorithmwillfind The quantum lattice Boltzmann (QLB) approach is based
ausefulroleinsimulatingexperimentswithgraphene. on the structural similarity between Eq. (6) and the discrete
Boltzmannequationwithamatrixcollisionoperator,
(cid:10)n (cid:11) (cid:12)
II. DIRACEQUATION
(∂ +ξ ·∇)f = (cid:10) f −f(0) , (7)
t i i ij j j
The Dirac equation offers a quantum-mechanical descrip- j=0
tion of an electron that is compatible with special relativity.
for i =0,...,n. In particular, Eqs. (6) and (7) each contain
Theformwetakeforourstartingpointis
only first derivatives in space and time. The distribution
(∂ t +cα·∇)(cid:3) =−iω c β(cid:3)+ig(cid:3), (1) functionsf i areanalogoustothefourcomponentsofthewave
function ψ, the set of discrete velocities ξ is analogous to
i
wherecisthelightspeed,h¯ isthereducedPlanck’sconstant, the set of streaming matrices αx, β, αz. The collision matrix
and ω c =mc2/h¯ is the Compton frequency for a particle of (cid:10) is analogous to the matrix αy in the Dirac equation.
ij
mass m. The wave function (cid:3) is a column vector with four
ThefirsthydrodynamiclatticeBoltzmannschemescontained
components, β is a 4×4 matrix, and α =(αx,αy,αz) is a
collisionmatricesderivedfromlinearizinglatticegascollision
collection of three 4×4 matrices, so that α·∇ =αx∂ +
x operators [32,33]. These were later replaced by scattering
αy∂ +αz∂ . The last term couples the wave function to an
y z matrices constructed from a known, rotationally symmetric
appliedscalarpotentialV viathecoefficientg =qV/h¯,where
set of eigenvectors [5,34]. More recent work uses collision
q isthemodulusofanelectron’scharge.
matricesdefinedthroughtheiractiononabasisofmomentsof
In the standard representation [10] of the Dirac equation,
thedistributionfunctions[35–37].
theαandβ matricesare
Rewriting Eq. (7) in an explicit matrix form like Eqs. (1)
⎛ ⎞ ⎛ ⎞
0 0 0 1 0 0 0 −i and(6)gives
⎜ ⎟ ⎜ ⎟ ⎡ ⎛ ⎞ ⎛ ⎞
⎜0 0 1 0⎟ ⎜0 0 i 0 ⎟ ξ 0 ··· 0 ξ 0 ··· 0
αx =⎜ ⎟, αy =⎜ ⎟, 0x 0y
⎝0 1 0 0⎠ ⎝0 −i 0 0 ⎠ ⎢ ⎢ ⎜ ⎜ 0 ξ ··· 0 ⎟ ⎟ ⎜ ⎜ 0 ξ ··· 0 ⎟ ⎟
⎛ 1 0 0 0 ⎞ ⎛ i 0 0 0 ⎞(2)
⎢
⎢ ⎣ ∂ t
+⎜
⎜ ⎝ . . .
1
. . .
x
... . . .
⎟
⎟ ⎠ ∂ x
+⎜
⎜ ⎝ . . .
1
. . .
y
... . . .
⎟
⎟ ⎠ ∂ y
0 0 1 0 1 0 0 0
⎜ ⎜0 0 0 −1 ⎟ ⎟ ⎜ ⎜0 1 0 0 ⎟ ⎟ ⎛ 0 0 ··· ξ nx ⎞ ⎤⎛ ⎞ 0 0 ⎛ ··· ξ ny⎞
Th α e z th = re ⎜ ⎝ eα 1 0 m − a 0 tr 1 ices 0 0 may 0 0 b ⎟ ⎠ e , wri β tte = ni ⎜ ⎝ nb 0 0 loc 0 0 kfo − r 0 m 1 as − 0 1 ⎟ ⎠ . + ⎜ ⎜ ⎜ ⎜ ⎝ ξ 0 0 . . . z ξ 0 1 . . . z · · . · · .. · · 0 0 . . . ⎟ ⎟ ⎟ ⎟ ⎠ ∂ z ⎥ ⎥ ⎥ ⎥ ⎦ ⎜ ⎜ ⎜ ⎜ ⎝ f f . . . 0 1 ⎟ ⎟ ⎟ ⎟ ⎠ =(cid:4) ⎜ ⎜ ⎜ ⎜ ⎜ ⎝ f f 0 1 − − . . . f f 0 1 ( ( 0 0 ) ) ⎟ ⎟ ⎟ ⎟ ⎟ ⎠ .
(cid:8) (cid:9)
αi = 0 σi , (3) 0 0 ··· ξ nz f n f n −f n (0)
σi 0 (8)
046706-2

ISOTROPYOFTHREE-DIMENSIONALQUANTUMLATTICE... PHYSICALREVIEWE83,046706(2011)
Each of the matrices multiplying ∂ , ∂ , ∂ is diagonal, so a one-dimensionalDiracequations,
x y z
formal solution to Eq. (8) is readily obtained by integrating
along the characteristics defined by the vectors ξ i , as in ∂ t u 1,2 +c∂ z u 1,2 =ω c d 2,1 +igu 1,2 , (15)
Ref. [38]. He and Luo [39,40] gave another formal solution ∂ d −c∂ d =−ω u +igd ,
t 1,2 z 1,2 c 2,1 1,2
basedonremovingthealgebraictermwithanintegratingfactor
beforeintegratingalongcharacteristics.Neitherapproachcan forthevariables(u ,d )and(u ,d )thatmakeuptherotated
1 2 2 1
beappliedtotheDiracequation,becausethethreestreaming wave function Z−1ψ =(u ,u ,d ,d )T. The components u
1 2 1 2
matrices cannot be simultaneously diagonalized in either the and d propagate up and down the z axis, respectively, and
standard orthe Majorana forms.Instead,we fallback onthe thesubscriptsindicatethespinup(1)andspindown(2)states,
technique of operator splitting to approximate solutions to respectively.Thevariablesu andd arethustheRiemann
1,2 1,2
Eq.(6)byevolvingthesolutioninonespatialdimensionata invariants [13] of the differential operator for streaming
time. alongz.
As observed in Ref. [6], the system (15) may be treated
as a Boltzmann equation for a pair of complex distribution
III. ONE-DIMENSIONALQUANTUMLATTICE
functionsu andd .Equation(15)maythusbediscretized
BOLTZMANNMODEL 1,2 1,2
by following the same approach that leads from the discrete
The y streaming matrix β is diagonal in the Majorana
Boltzmannequation(7)toalatticeBoltzmannscheme[38].
form (6) of the three-dimensional Dirac equation. However,
The left hand sides of Eq. (15) are derivatives along
theotherstreamingmatricesmaybediagonalizedbyapplying the characteristics given by the light cones (z±cs,t +s)
unitarytransformations.Thematrix−αzmaybediagonalized
parametrized by s. Integrating along these characteristics
as therefore gives exactly the difference in variable u or d
⎛ ⎞ 1,2 1,2
1 0 0 0 betweenthetwoendsofthecharacteristic.Therighthandside
⎜ ⎟ of Eq. (15) is discretized using the trapezoidal rule, which
⎜0 1 0 0 ⎟
Z −1(−αz)Z =⎜ ⎟, (9) involvestheaverageofthevaluesattimet andpositionzand
⎝0 0 −1 0 ⎠ timet +(cid:12)t andpositionz±(cid:12)z.Weevaluateu atz+(cid:12)z
1,2
0 0 0 −1 andd atz−(cid:12)z,sincethesearethepointsappearinginthe
1,2
discretization of the left hand side. The resulting system of
usingtheunitarymatrix
⎛ ⎞ algebraicequationsis
0 −1 0 1
1 ⎜ ⎜1 0 −1 0 ⎟ ⎟ (cid:19)u 1,2 −u 1,2 = 1 2 m(cid:20)(d 2,1 +d (cid:19) 2,1 )+ 1 2 i(cid:20)g(u 1,2 +(cid:19)u 1,2 ),
Z = √ 2 ⎜ ⎝0 1 0 1 ⎟ ⎠ . (10) d (cid:19) 1,2 −d 1,2 =−1 2 m(cid:20)(u 2,1 +(cid:19)u 2,1 )+ 1 2 i(cid:20)g(d 1,2 +d (cid:19) 1,2 ), (16)
1 0 1 0 whereasuperscripthat((cid:19))indicatesthatavariableisevaluated
Applying this to Eq. (6), we multiply on the left by Z−1 to somewhereotherthan(z,t),
obtain (cid:19)u =u (z+(cid:12)z,t +(cid:12)t), u =u (z,t),
1,2 1,2 1,2 1,2
Z −1[∂ t +c(−αx∂ x +β∂ y −αz∂ z )+iω c αy −igI]ψ =0, d (cid:19) =d (z−(cid:12)z,t +(cid:12)t), d =d (z,t). (17)
1,2 1,2 1,2 1,2
(11)
ThedimensionlessComptonfrequencyism(cid:20)=ω (cid:12)t,andthe
andinsertafactorofZZ−1 =I,
dimensionlessscalarpotentialis(cid:20)g =g(z,t)(cid:12)t.
c
Z −1[∂ +c(−αx∂ +β∂ −αz∂ ) The pair of equations (16) may be solved algebraically to
t x y z obtainexplicitformulasfor(cid:19)u andd (cid:19) ,
+iω αy −igI]ZZ −1ψ =0. (12) 1,2 1,2
c
Inthiswaythestreamingmatrixalongzisdiagonalized, (cid:19)u 1,2 =au 1,2 +bd 2,1 ,
(18)
[∂ +cZ −1(−αz)Z∂ +Z −1(−cαx∂ +cβ∂ d (cid:19) 1,2 =ad 1,2 −bu 2,1 ,
t z x y
+iω αy −igI)Z]Z −1ψ =0, (13) wherethecoefficientsaandbare
c
whiletherotatedcollisionmatricesare a =
1−(cid:10)/4
, b=
m(cid:20)
, (cid:10)=m(cid:20)2−(cid:20)g2.
−Z −1(iω αy −igI)Z =ω βαx +igI
1+(cid:10)/4−i(cid:20)g 1+(cid:10)/4−i(cid:20)g
c ⎛c ⎞
ig 0 0 ω These coefficients satisfy |a|2+|b|2 =1, so the right hand
c
⎜ ⎟ side of Eq. (18) corresponds to multiplying the rotated
⎜ 0 ig ω 0⎟
=⎜ ⎝ 0 −ω c ig c 0 ⎟ ⎠ . (14) m wa a v tr e ix fu , nctionZ−1ψ =(u 1 ,u 2 ,d 1 ,d 2 )T bytheunitarycollision
−ω 0 0 ig ⎛ ⎞
c a 0 0 b
Weinsertaminussignthroughoutbecausetherotatedcollision ⎜ ⎟
⎜ 0 a b 0⎟
matrices appear on the right hand sides of (15) and (16). Q=⎜ ⎟. (19)
⎝ 0 −b a 0⎠
Neglectinganydependenceofψ onthex andy coordinates,
thecombinationEqs.(13)and(14)maybewrittenasapairof −b 0 0 a
046706-3

P.J.DELLAR,D.LAPITSKI,S.PALPACELLI,ANDS.SUCCI PHYSICALREVIEWE83,046706(2011)
The streaming step propagates u upward, and d down- alongyisalreadydiagonalintheMajoranaform,soY =I is
1,2 1,2
ward, along the light cones given by (cid:12)z=±c(cid:12)t. This is theidentitymatrix.ThematrixXis
also a unitary operation, so the overall QLB scheme evolves ⎛ ⎞
−1 0 1 0
the discrete wave function through a sequence of unitary
⎜ ⎟
operations. 1 ⎜ 0 1 0 −1⎟
X = √ ⎜ ⎟, (22)
2 ⎝ 1 0 1 0 ⎠
0 1 0 1
IV. ATHREE-DIMENSIONALSCHEMEVIA
OPERATORSPLITTING andtheZmatrixisgiveninEq.(10)above.
Wehavechosentosplitthecollisiontermintothreeparts,
We use operator splitting to decompose the three-
eachofwhichiscombinedwithoneofthestreamingsteps.The
dimensional Dirac equation into the sum of three one-
collisionmatrixthuscoincides,uptoaunitarytransformation,
dimensional equations, each involving spatial derivatives
withthecollisionmatrixfortheone-dimensionalQLBscheme
along a single direction. We thus write the Majorana form with a time step of 1(cid:12)t (see Ref. [8]). In particular, Q (cid:19) is
oftheDiracequationas 3
givenby
⎛ ⎞
(cid:21) (cid:22) (cid:23) (cid:19)a 0 0 −(cid:19) b
∂ t
+
+ (cid:22)
cβ
c(
∂
−
y
α
+
x)
1 3
∂
i
x
(ω
+
c α
1 3
y
i(
+
ω c
g
α
I
y
)
(cid:23) +gI)
Q
(cid:19)= ⎜ ⎜ ⎜
⎜ ⎝0
0
−
(cid:19)a
(cid:19) b (cid:19)
(cid:19)
a
b 0
0
⎟ ⎟ ⎟
⎟ ⎠ , (23)
(cid:22) (cid:23)(cid:24)
+ c(−αz)∂ + 1i(ω αy +gI) ψ =0, (20) (cid:19) b 0 0 (cid:19)a
z 3 c
wherethecoefficients
w di h re e c re tio e n ac a h nd te c r o m llis [· io · n ·] . W is e an spl o i p t e th ra e to c r ol f l o is r io s n tre te a r m m in i g nto in th o r n ee e (cid:19)a = 1+ 1 (cid:10) − / (cid:10) 4 3 − /4 i(cid:20)g , (cid:19) b= 1+(cid:10) m(cid:20) /4 3 −i(cid:20)g ,
equalpieces,andincludeonepieceineachbracketedtermto 3 3 3 3
matchtheformoftheone-dimensionalDiracequationofthe are written in terms of the rescaled dimensionless Compton
lastsection. andpotentialfrequencies,
Weuseoperatorsplittingtoapproximatetheevolutionover
(cid:10) =m(cid:20)2−(cid:20)g2, m(cid:20) = 1ω (cid:12)t, (cid:20)g = 1g(cid:12)t.
a time step (cid:12)t by decomposing Eq. (20) into three separate 3 3 3 3 3 c 3 3
equations.Usingthenotationexp{(cid:12)t[···]}ψ 0 forthesolution Thepatternof+and−signsinthe (cid:19) btermsontheoff-diagonal
of the linear evolution equation ∂ t ψ =[···]ψ with initial of Q (cid:19) follows the same pattern as the αy matrix. The rotated
condition ψ 0 over a time interval (cid:12)t, our operator splitting matrices X−1Q (cid:19) X and Z−1Q (cid:19) Z have the same sign pattern as
correspondstocomputing (cid:19)
Q,butQitselfdoesnot.
(cid:21) (cid:22) (cid:23)(cid:24)
ψ (cid:5) =exp (cid:12)t c(−αx)∂ + 1i(ω αy +gI) ψ(t), (21a) V. NUMERICALRESULTS
x 3 c
(cid:21) (cid:22) (cid:23)(cid:24) To verify the isotropic nature of the scheme we choose
ψ (cid:5)(cid:5) =exp (cid:12)t cβ∂ + 1i(ω αy +gI) ψ (cid:5) , (21b)
y 3 c initial conditions in which the positive energy, spin-up com-
(cid:21) (cid:22) (cid:23)(cid:24) ponent φ + is a spherically symmetric Gaussian wave packet
ψ(t +(cid:12)t)=exp (cid:12)t c(−αz)∂ + 1i(ω αy +gI) ψ (cid:5)(cid:5) . 1
z 3 c withspread(cid:12) 0 ,
(21c) (cid:8) (cid:9)
φ + (x,y,z,t)=
(cid:11)
2π(cid:12)2
(cid:12)
−3/4 exp −
x2+y2+z2
. (24)
1 0 4(cid:12)2
The three operators do not commute, so we incur an 0
O((cid:12)t2)splittingerrorbyusingthisproductofexponentialsto The other three components φ + and φ − are initially set to
2 1,2
approximatetheexponentialofthesumofthethreeoperators. zero. As a first measure of isotropy, we study the temporal
EachofthethreestagesinEqs.(21)denotingevolutionbya evolution of the standard deviations of |φ +|2 in each of the
1
timestep(cid:12)t isaccomplishedbyrotatingψ todiagonalizethe threecoordinatedirections.Thesequantitiesaredefinedby
relevantstreamingmatrix,takingonetimestepoftheexisting (cid:8)(cid:25) (cid:9) (cid:26)(cid:8)(cid:25) (cid:9)
1/2 1/2
one-dimensional QLB scheme described above, and rotating (cid:12) = x2|φ +|2dV |φ +|2dV , (25)
ψ back toitsoriginalbasis.Ouralgorithmisthus composed x 1 1
ofthefollowingsteps
(i) Rotateψ withX−1,collidewithX−1Q (cid:19) X,streamalong and similarly for (cid:12) y and (cid:12) z . We calculate discrete approxi-
mationstotheseintegralsusingthetrapezoidalrule,
x,rotatebackwithX.
(ii) Rotateψ withY−1,collidewithY−1Q (cid:19) Y,streamalong ⎛ ⎞ 1/2 (cid:26)⎛ ⎞ 1/2
(cid:10) (cid:10)
y ( , i r ii o ) ta R te o b ta a t c e k ψ w w ith ith Y. Z−1,collidewithZ−1Q (cid:19) Z,streamalong (cid:12) x =⎝ x i 2|φ 1 +|2dV ⎠ ⎝ |φ 1 +|2dV ⎠ . (26)
i,j,k i,j,k
z,rotatebackwithZ.
We have written the algorithm like this to emphasize the Althoughusuallyonlysecond-orderaccurate,thetrapezoidal
symmetrybetweenthethreesteps,butthematrixforstreaming rulebecomesexponentiallyaccuratewhenthesumsaretaken
046706-4

ISOTROPYOFTHREE-DIMENSIONALQUANTUMLATTICE... PHYSICALREVIEWE83,046706(2011)
overgridsofequallyspacedpointsinadomainwithperiodic
|     |     |     |     |     |     |     |     |     |     | 0.03 |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- |
boundaryconditions[41].
128
256
|     |     |     |     |              |     |     |     |     |     | 0.025 | 512 |     |     |     |
| --- | --- | --- | --- | ------------ | --- | --- | --- | --- | --- | ----- | --- | --- | --- | --- |
|     |     |     | A.  | Freeparticle |     |     |     |     |     |       |     |     |     |     |
1024
ThesolutionoftheSchro¨dingerequationforafreeparticle
=0)withinitialconditionsgivenbyEq.(24)is
(V
0.02
|     |     |     |     | (cid:11) |     | (cid:12) |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | -------- | --- | -------- | --- | --- | --- | --- | --- | --- | --- | --- |
− x2+y2+z2
exp
|     | φ(x,y,z,t)= |     | √   |               | 4(cid:12)2 +2ih¯t/m     |     |     |      |          |       |     |     |     |     |
| --- | ----------- | --- | --- | ------------- | ----------------------- | --- | --- | ---- | -------- | ----- | --- | --- | --- | --- |
|     |             |     |     |               | 0                       |     | .   | (27) |          |       |     |     |     |     |
|     |             |     |     |               | +ih¯t/(2m(cid:12))])3/2 |     |     |      | Δ(cid:2) |       |     |     |     |     |
|     |             |     | (   | 2π[(cid:12) 0 |                         |     |     |      |          | 0.015 |     |     |     |     |
Thethreespreadsarethusgivenby
|     |     |          |      | (cid:8)        |       | (cid:9) |     |      |     |      |     |     |     |     |
| --- | --- | -------- | ---- | -------------- | ----- | ------- | --- | ---- | --- | ---- | --- | --- | --- | --- |
|     |     |          |      |                | h¯2t2 | 1/2     |     |      |     | 0.01 |     |     |     |     |
|     |     |          | (t)= | (cid:12)2+     |       |         |     |      |     |      |     |     |     |     |
|     |     | (cid:12) |      |                |       | ,       |     | (28) |     |      |     |     |     |     |
|     |     |          | α    | 0 4m2(cid:12)2 |       |         |     |      |     |      |     |     |     |     |
0
=x,y,z.
| for                               | α   | Figure |     | 1 shows | the measured       |     | spreads | for a |     | 0.005 |     |     |     |     |
| --------------------------------- | --- | ------ | --- | ------- | ------------------ | --- | ------- | ----- | --- | ----- | --- | --- | --- | --- |
| freeparticlewithm=0.35and(cid:12) |     |        |     |         | =14inacubewithside |     |         |       |     |       |     |     |     |     |
0
| length(cid:15)=100discretizedusing1283 |          |     |     |     | points.Weusenatural |         |     |                           |     |     |     |     |     |     |
| -------------------------------------- | -------- | --- | --- | --- | ------------------- | ------- | --- | ------------------------- | --- | --- | --- | --- | --- | --- |
|                                        |          | c=1 |     | =1. |                     |         |     |                           |     | 0   |     |     |     |     |
| units                                  | in which |     | and | h¯  | The three           | spreads |     | (cid:12) x , (cid:12) y , |     |     |     |     |     |     |
|                                        |          |     |     |     |                     |         |     |                           |     | 0   | 50  | 100 | 150 | 200 |
(cid:12) are in good agreement with each other, and follow the t
z
generaltrendoftheSchro¨dingersolutionasgivenbyEq.(28).
However,thecomputedspreadsdeviatefromtheSchro¨dinger (cid:20)
|     |     |     |     |     |     |     |     |     |     | FIG.2. (Color | online) Evolution | of the average | difference | (cid:12) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ------------- | ----------------- | -------------- | ---------- | -------- |
solutionthroughsuperimposedsmallamplitudeoscillationsat between the three spreads for a wave packet with initial spread
theComptonfrequency.Theseoscillationsareintrinsictothe (cid:12) =14, mass m=0.35, and no potential in a cube with side
0
Dirac equation, and are due to relativistic effects associated length (cid:15)=100. The cube was discretized using 1283, 2563, 5123,
and10243points.
withaparticleoffinitemass.Inaforthcomingpublicationwe
showthattheQLBsolutionisaconsistentnumericalsolution
oftheoriginalDiracequation[26].
Figure2showsthedeviationfromisotropyasmeasuredby refinement.Figure3showsthattherescaledspreaddifferences
|     |     |     |     |     |     |     |     |     | n2(cid:12) | (cid:20) |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ---------- | -------- | --- | --- | --- | --- |
themeandifferencebetweenthethreespreads, forthefourdifferentgridsagreeveryclosely,confirming
(cid:20) =O(n−2).
|     |            |             |           |                       |             |           |     |     | that | (cid:12) | We thus obtain | second-order | convergence |     |
| --- | ---------- | ----------- | --------- | --------------------- | ----------- | --------- | --- | --- | ---- | -------- | -------------- | ------------ | ----------- | --- |
|     | (cid:20) = | 1(|(cid:12) | −(cid:12) | |+|(cid:12) −(cid:12) | |+|(cid:12) | −(cid:12) | |), |     |      |          |                |              |             |     |
(cid:12) x y x z y z toward isotropy, at least as measured by invariance under
3
|     |     |     |     |     |     |     |     |     | 90◦ | rotations. |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---------- | --- | --- | --- | --- |
forthreedifferentcomputationalgridswith1283,2563,5123,
|     | 10243 |         |           | (cid:20) |          |          |       |      |     |     |     |     |     |     |
| --- | ----- | ------- | --------- | -------- | -------- | -------- | ----- | ---- | --- | --- | --- | --- | --- | --- |
| and |       | points. | As shown, | (cid:12) | tends to | decrease | under | grid |     |     |     |     |     |     |
|     |       |         |           |          |          |          |       |      |     | 450 |     |     |     |     |
26
128
Δ
|     |     |     | x   |     |     |     |     |     |     | 400 | 256 |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Δ
|     | 24  |     | y           |     |     |     |     |     |     |     | 512  |     |     |     |
| --- | --- | --- | ----------- | --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- |
|     |     |     | Δ           |     |     |     |     |     |     |     | 1024 |     |     |     |
|     |     |     | z           |     |     |     |     |     |     | 350 |      |     |     |     |
|     | 22  |     | Schrödinger |     |     |     |     |     |     |     |      |     |     |     |
300
|     | 20  |     |     |     |     |     |     |     |     | 250 |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Δ(cid:2)2n
Δ
|     | 18  |     |     |     |     |     |     |     |     | 200 |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
150
16
100
14
50
12
|     |     | 0   | 50  | 100 |     | 150 |     | 200 |     | 0   |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
t
|     |     |     |     |     |     |     |     |     |     | 0   | 50  | 100 | 150 | 200 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
t
|     | FIG.1.                                       | (Color | online) | Dispersion | of a | wave packet | with | initial |     |     |     |          |     |     |
| --- | -------------------------------------------- | ------ | ------- | ---------- | ---- | ----------- | ---- | ------- | --- | --- | --- | -------- | --- | --- |
|     | =14,particlemassm=0.35,andnopotentialinacube |        |         |            |      |             |      |         |     |     |     | (cid:20) |     |     |
spread(cid:12) FIG.3. (Coloronline)Evolutionofn2(cid:12) forawavepacketwith
0
withsidelength(cid:15)=100andacomputationalgridwith1283points. initialspread(cid:12) =14,massm=0.35,andnopotentialinacube
|     |     |     |     |     |     |     |     |     |     |     | 0   |     | (cid:20) |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | -------- | --- |
Thethreecomputedspreads(cid:12) ,(cid:12) ,(cid:12) collapseontoasinglecurve, withsidelength(cid:15)=100.Thescaledmeananisotropiesn2(cid:12) collapse
|     |     |     |     | x y | z   |     |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
whichdiffersbyasuperimposedhighfrequencyoscillationfromthe ontoasinglecurveforfourdifferentcomputations usingn3 points
withn=128,256,512,1024.
Schro¨dingersolution(28)(smoothcurve).
046706-5

P.J.DELLAR,D.LAPITSKI,S.PALPACELLI,ANDS.SUCCI PHYSICALREVIEWE83,046706(2011)
|     | 19  |     |     |     |     |     |     |     | 0.45 |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- | --- | --- |
Δ
|     |     |     |     |     |     |     | x   |     |     |     |     | 128 |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     |     |     |     |     |     | Δ   |     | 0.4 |     |     |     |     |     |     |
256
y
|     | 18  |     |     |     |     |     | Δ   |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
512
|     |     |     |     |     |     |     | z   |     | 0.35 |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- | --- | --- |
1024
|     | 17  |     |     |     |     |     |     |     | 0.3 |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
0.25
| Δ   | 16  |     |     |     |     |     |     | Δ(cid:2) |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | -------- | --- | --- | --- | --- | --- | --- | --- |
0.2
0.15
15
0.1
14
0.05
|     |        |                |           |        |          |         |            |     | 0             |                          |           |                          |                     |        |      |
| --- | ------ | -------------- | --------- | ------ | -------- | ------- | ---------- | --- | ------------- | ------------------------ | --------- | ------------------------ | ------------------- | ------ | ---- |
|     | 13     |                |           |        |          |         |            |     | 0             | 100                      |           | 200                      | 300                 |        | 400  |
|     | 0      | 100            |           | 200    |          | 300     | 400        |     |               |                          |           |                          |                     |        |      |
|     |        |                |           | t      |          |         |            |     |               |                          |           | t                        |                     |        |      |
|     |        |                |           |        |          |         |            |     | FIG.5. (Color | online)                  | Evolution | of (cid:12) (cid:20) for | a wave              | packet | with |
|     | FIG.4. | (Color online) | Evolution | of the | computed | spreads | (cid:12) , |     |               |                          |           |                          |                     |        |      |
|     |        |                |           |        |          |         | x          |     |               | =14,massm=0.1,potentialω |           |                          | =1/(2m(cid:12)2),in |        |      |
(cid:12) ,(cid:12) forawavepacketwithinitialspread(cid:12) =14,particlemass initialspread(cid:12) 0 0
y z 0 acubewithsidelength(cid:15)=200.Thecubewasdiscretizedusing1283, 0
| m=0.1,potentialω |     |     | =1/(2m(cid:12)2),inacubewithsidelength(cid:15)= |     |     |     |     |     |     |     |     |     |     |     |     |
| ---------------- | --- | --- | ----------------------------------------------- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
0
|     |     |     |     | 0   |     |     |     | 2563,5123,and10243points. |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | ------------------------- | --- | --- | --- | --- | --- | --- | --- |
200discretizedwith1283points.
|     |     |     | B. Harmonicpotential |     |     |     |     |     |     | C. Convergencetosphericalsymmetry |     |     |     |     |     |
| --- | --- | --- | -------------------- | --- | --- | --- | --- | --- | --- | --------------------------------- | --- | --- | --- | --- | --- |
In this second test case, we consider the Gaussian wave A more demanding test of isotropy is given by testing
|        |     |          |          |      |             |           |     | the                                     | convergence | of the | numerical | solutions | toward | spherical    |     |
| ------ | --- | -------- | -------- | ---- | ----------- | --------- | --- | --------------------------------------- | ----------- | ------ | --------- | --------- | ------ | ------------ | --- |
| packet | of  | Eq. (24) | confined | in a | spherically | symmetric |     |                                         |             |        |           |           |        |              |     |
|        |     |          |          |      |             |           |     | symmetry,ratherthanjustsymmetryunder90◦ |             |        |           |           |        | rotations.We |     |
harmonicpotential,
|     |     |     |               |     |     |     |     | divide | the              | computational | domain  | up into         | spherical |        | shells, |
| --- | --- | --- | ------------- | --- | --- | --- | --- | ------ | ---------------- | ------------- | ------- | --------------- | --------- | ------ | ------- |
|     |     |     |               |     |     |     |     |        | ∈[r +(cid:12)r], |               |         |                 |           | |φ     | +|      |
|     |     |     | V(r)=−1mω2r2, |     |     |     |     | r      | ,r               | and           | perform | a least squares |           | fit of | to      |
|     |     |     |               | 0   |     |     |     |        | i i              |               |         |                 |           |        | 1       |
2
r2 =x2+y2+z2,
| in our | sign      | convention, | where |          |              |     | and ω 0 is |     |     |     |     |     |     |     |     |
| ------ | --------- | ----------- | ----- | -------- | ------------ | --- | ---------- | --- | --- | --- | --- | --- | --- | --- | --- |
|        |           |             |       |          |              |     |            |     | 110 |     |     |     |     |     |     |
| the    | frequency | associated  | with  | harmonic | oscillations |     | in this    |     |     |     |     |     |     |     |     |
128
| potential. |     | Again, the | analytical | solution | of  | the Schro¨dinger |     |     |     |     |     |     |     |     |     |
| ---------- | --- | ---------- | ---------- | -------- | --- | ---------------- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|            |     |            |            |          |     |                  |     |     | 100 | 256 |     |     |     |     |     |
equationforthispotentialiswellknown,andsetting
512
90
1024
1
|     |     |     | ω   | =   |     |     |     |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     |     | 0   |     |     |     |     |     | 80  |     |     |     |     |     |     |
2m(cid:12)2
0
70
| leads     | to  | the initial   | spread | (cid:12) 0 of the | wave | packet           | being |           |     |     |     |     |     |     |     |
| --------- | --- | ------------- | ------ | ----------------- | ---- | ---------------- | ----- | --------- | --- | --- | --- | --- | --- | --- | --- |
| preserved |     | by subsequent |        | evolution under   |      | the Schro¨dinger |       |           | 60  |     |     |     |     |     |     |
| equation. |     |               |        |                   |      |                  |       | Δ(cid:2)n |     |     |     |     |     |     |     |
50
Weinvestigatetheabilityoftheschemetopreserveisotropy
|      |                                |            | m=0.1       |                | =14            |       |          |     |     |     |     |     |     |     |     |
| ---- | ------------------------------ | ---------- | ----------- | -------------- | -------------- | ----- | -------- | --- | --- | --- | --- | --- | --- | --- | --- |
| of   | the solution.                  | We         | set         | and            | (cid:12)       | thus  | yielding |     |     |     |     |     |     |     |     |
|      |                                |            |             |                | 0              |       |          |     | 40  |     |     |     |     |     |     |
| ω    | =0.0255andhenceaperiodofaboutT |            |             |                |                | =2π/ω | ∼250     |     |     |     |     |     |     |     |     |
| 0    |                                |            |             |                |                |       | 0        |     |     |     |     |     |     |     |     |
| time | units                          | based on   | a box       | side of length | (cid:15)=200   | and   | c=1.     |     | 30  |     |     |     |     |     |     |
| As   | before,                        | we perform | simulations | with           | discretization |       | using    |     |     |     |     |     |     |     |     |
20
| 1283,2563,5123,and10243 |                                  |     |     | points.   |              |             |     |     |     |     |     |     |     |     |     |
| ----------------------- | -------------------------------- | --- | --- | --------- | ------------ | ----------- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|                         | InFig.4,weshowthespreads(cid:12) |     |     | ,(cid:12) | ,and(cid:12) | asfunctions |     |     |     |     |     |     |     |     |     |
|                         |                                  |     |     | x         | y            | z           |     |     | 10  |     |     |     |     |     |     |
oftime.Theyareingoodagreementwitheachother,butshow
|                                                          |            |      |     | (cid:12)=14 |     |          |        |     |     |     |     |     |     |     |     |
| -------------------------------------------------------- | ---------- | ---- | --- | ----------- | --- | -------- | ------ | --- | --- | --- | --- | --- | --- | --- | --- |
| large                                                    | excursions | from | the | value       | one | would    | expect |     | 0   |     |     |     |     |     |     |
|                                                          |            |      |     |             |     |          |        |     | 0   | 100 |     | 200 | 300 |     | 400 |
| themtotakeifthewavepacketevolvedundertheSchro¨dinger     |            |      |     |             |     |          |        |     |     |     |     | t   |     |     |     |
| equation.Figure5showsthemeandifference(cid:12)betweenthe |            |      |     |             |     | (cid:20) |        |     |     |     |     |     |     |     |     |
threespreads.Weobservefirst-orderconvergencetoward90◦ (cid:20)
|     |     |     |     |     |     |     |     |     | FIG.6. (Color | online) | Evolution | of n(cid:12) | for a wave | packet | with |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ------------- | ------- | --------- | ------------ | ---------- | ------ | ---- |
rotationalsymmetryfort (cid:2)80,asconfirmedbytherescaled initialspread(cid:12) =14,massm=0.1,potentialω =1/(2m(cid:12)2),ina
|     |     |     |     |     |     |     |     |     |     | 0   |     |     | 0   |     | 0   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
plot of n(cid:12) (cid:20) in Fig. 6, but no convergence is visible at later cubewithsidelength(cid:15)=200.Thecubewasdiscretizedusing1283,
(cid:20)
times. This appears to be due to the sensitive dependence of 2563,5123,and10243points.Thescaledanistropiesn(cid:12) collapseinto
thesolutionsonspatialresolution,asdescribedinRef.[26]. asinglecurvefort (cid:2)60.
046706-6

ISOTROPYOFTHREE-DIMENSIONALQUANTUMLATTICE... PHYSICALREVIEWE83,046706(2011)
|     | 0.014 |     |     |     |     |     |     |     | 16  |     |     |     |     |     |
| --- | ----- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |       |     |     |     |     | 128 |     |     |     | Δ   |     |     |     |     |
x
256
|     | 0.012 |     |     |     |     |     |     |     |     | Δ   |     |     |     |     |
| --- | ----- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |       |     |     |     |     | 512 |     |     |     |     | y   |     |     |     |
15.5
1024
Δ z
0.01
15
0.008
| σ   |     |     |     |     |     |     |     |     | 14.5 |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- | --- |
Δ
0.006
14
0.004
|     | 0.002 |     |     |     |     |     |     |     | 13.5 |     |     |     |     |     |
| --- | ----- | --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- | --- |
0
|     | 0   |     | 50  |     | 100 |     | 150 |     | 13   |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- | --- |
|     |     |     |     | r   |     |     |     |     | 0    | 100 | 200 |     | 300 | 400 |
t
|        |        |         |      |         | =100 |     | n=  |     |     |     |     |     |     |     |
| ------ | ------ | ------- | ---- | ------- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| FIG.7. | (Color | online) | Plot | of σ(r) | at t | for |     |     |     |     |     |     |     |     |
128,256,512,1024inacubeofside(cid:15)=200.Thefeatureatr ≈100is FIG.9. (Color online) Dispersion of a wave packet with initial
|     |     |     |     |     |     |     |     |     |     | =14, | m=0.1, |     | =1/(2m(cid:12)2), |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---- | ------ | --- | ----------------- | --- |
duetocouplingwithnegativeenergystatesthathavepassedthrough spread (cid:12) 0 particle mass potential ω 0
|     |     |     |     |     |     |     |     |     |     |     | (cid:15)=100 |     |     | 0   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ------------ | --- | --- | --- |
the periodic boundaries. The 1283 grid contains too few points to in a cube with side length and a computational grid with
computeσ forr <3. 1283points.Computedspreads(cid:12) ,(cid:12) ,(cid:12) .
|     |     |     |     |     |     |     |     |     |     |     | x   | y z |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
+b
alinearfunctiona r withineachshell.Wethusconstruct showsthefunctionσ(r)forthecomputedsolutionsatt =100
|     |     | i   | i   |     |     |     |     |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
a piecewise-linear and spherically symmetric solution by and a bin size (cid:12)r =1. The numerical solutions appear to
averaging the computed solution over spherical shells. The (cid:2)50.
|                                      |     |     |     |     |                |     |     | converge |         | toward spherical | symmetry | in the      | region | r      |
| ------------------------------------ | --- | --- | --- | --- | -------------- | --- | --- | -------- | ------- | ---------------- | -------- | ----------- | ------ | ------ |
| standarddeviationσ(r)ofthecomputed|φ |     |     |     |     | +|awayfromthis |     |     |          |         | =100             |          |             |        |        |
|                                      |     |     |     |     |                |     |     | The      | feature | near r           | is due   | to negative | energy | states |
1
spherical average gives a measure of the deviation of the being drawn toward the periodic boundaries of the cube by
computed solution away from spherical symmetry. Figure 7 thepotential.(Ascalarpotentialthatisconfiningforpositive
|     |     |     |     |     |     |     |     | energy | states | is deconfining | for negative | energy | states, | and |
| --- | --- | --- | --- | --- | --- | --- | --- | ------ | ------ | -------------- | ------------ | ------ | ------- | --- |
|     | 1.6 |     |     |     |     |     |     |        |        |                |              |        |         |     |
128
|     |     |     |     |     |     |     |     |     | 30  |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
256
|     | 1.4 |     |     |     |     |      |     |     |     | 128 |     |     |     |     |
| --- | --- | --- | --- | --- | --- | ---- | --- | --- | --- | --- | --- | --- | --- | --- |
|     |     |     |     |     |     | 512  |     |     |     | 256 |     |     |     |     |
|     |     |     |     |     |     | 1024 |     |     |     | 512 |     |     |     |     |
|     | 1.2 |     |     |     |     |      |     |     | 25  |     |     |     |     |     |
1024
1
20
σn
0.8
Δ(cid:2)n
15
0.6
10
0.4
0.2
5
0
|     | 0   |     | 50  |     | 100 | 150 |     |     | 0   |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
r
|     |     |     |     |     |     |     |     |     | 0   | 100 | 200 |     | 300 | 400 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
t
|        |        |         |      |          | =100 |     | n=  |     |     |     |     |     |     |     |
| ------ | ------ | ------- | ---- | -------- | ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| FIG.8. | (Color | online) | Plot | of nσ(r) | at t | for |     |     |     |     |     |     |     |     |
(cid:20)
128,256,512,1024 confirming first-order convergence to spherical FIG.10. (Coloronline)Evolutionofn(cid:12)forawavepacketwith
symmetryundergridrefinementintheregionr (cid:2)50.Thenoncon- initialspread(cid:12) =14,massm=0.1,potentialω =1/(2m(cid:12)2),in
|     |     |     |     |     |     |     |     |     |     | 0   |     |     | 0   | 0   |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
vergingfeaturecenteredatr ≈100isduetocouplingwithnegative acubewithsidelength(cid:15)=100.Thecomputationalgridscontained
n3pointsforn=128,256,512,1024.
energystatesthathavepassedthroughtheperiodicboundaries.
046706-7

P.J.DELLAR,D.LAPITSKI,S.PALPACELLI,ANDS.SUCCI PHYSICALREVIEWE83,046706(2011)
|     |     |     |     |     |     |     | (cid:12) (cid:20) =O(n−1). |           | However, | the  | solution  | in this  | smaller | box is |
| --- | --- | --- | --- | --- | --- | --- | -------------------------- | --------- | -------- | ---- | --------- | -------- | ------- | ------ |
|     |     | 0.2 | 0.4 | 0.6 | 0.8 | 1   |                            |           |          |      |           |          |         |        |
|     |     |     |     |     |     |     | noticeably                 | distorted | away     | from | spherical | symmetry |         | by the |
+|
effectsoftheboundaries,asshownbythecontourplotof|φ
1
|     | 50  |     |     |     |     |     | ontheplanez=0inFig.11,andthusshowsnoconvergence |     |     |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | ----------------------------------------------- | --- | --- | --- | --- | --- | --- | --- |
tosphericalsymmetryundergridrefinement.
|     |     |     |     |     |     |     |     |     | VI. | CONCLUSIONS |     |     |     |     |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ----------- | --- | --- | --- | --- |
25
|     |     |     |     |     |     |     | The Dirac    |           | equation | has many | structural |           | similarities | with        |
| --- | --- | --- | --- | --- | --- | --- | ------------ | --------- | -------- | -------- | ---------- | --------- | ------------ | ----------- |
|     |     |     |     |     |     |     | the discrete | Boltzmann |          | equation | that       | underlies |              | the lattice |
Boltzmannapproachtosimulatinghydrodynamics.However,
|     |     |     |     |     |     |     | it is not | possible | to  | discretize | the three-dimensional |     |     | Dirac |
| --- | --- | --- | --- | --- | --- | --- | --------- | -------- | --- | ---------- | --------------------- | --- | --- | ----- |
y 0 equation by integration along characteristic curves, as in the
|     |     |     |     |     |     |     | derivation          | of lattice | Boltzmann |           | schemes  | for | hydrodynamics. |            |
| --- | --- | --- | --- | --- | --- | --- | ------------------- | ---------- | --------- | --------- | -------- | --- | -------------- | ---------- |
|     |     |     |     |     |     |     | The characteristics |            | of        | the Dirac | equation | are | not            | curves but |
lightcones,three-dimensionalsurfacesinspace-time.Instead,
|     |     |     |     |     |     |     | we employ | operator | splitting |     | to separate | the | solution | of the |
| --- | --- | --- | --- | --- | --- | --- | --------- | -------- | --------- | --- | ----------- | --- | -------- | ------ |
−25
three-dimensionalDiracequationintothreeone-dimensional
steps.Eachstepinvolvesderivativeswithrespecttoonlyone
|     |     |     |     |     |     |     | spatial coordinate,                                      |     | so                      | each step | may       | be discretized |             | using the |
| --- | --- | --- | --- | --- | --- | --- | -------------------------------------------------------- | --- | ----------------------- | --------- | --------- | -------------- | ----------- | --------- |
|     |     |     |     |     |     |     | one-dimensional                                          |     | quantum                 | lattice   | Boltzmann |                | scheme.     | How-      |
|     | −50 |     |     |     |     |     | ever,withineachstepthewavefunctionmustbetransformed      |     |                         |           |           |                |             |           |
|     | −50 |     | −25 | 0   | 25  | 50  |                                                          |     |                         |           |           |                |             |           |
|     |     |     |     | x   |     |     | intocharacteristicvariablesbyrotatingitusingoneoftheX,Y, |     |                         |           |           |                |             |           |
|     |     |     |     |     |     |     | Z matrices,evolved                                       |     | usingtheone-dimensional |           |           |                | scheme,then |           |
FIG.11. (Color online) Contours of |φ+(x,y,0)| on the plane rotatedbacktoitsoriginalformreadyforthenextstep.These
1
z=0 from a 5123 simulation in a cube of side length (cid:15)=100 rotations were applied in the incorrect order in the previous
|     |     |     |     |     |     |     | formulation | of  | the quantum |     | lattice Boltzmann |     | scheme | [6]. |
| --- | --- | --- | --- | --- | --- | --- | ----------- | --- | ----------- | --- | ----------------- | --- | ------ | ---- |
(valuesincreasefromoutercirclestoinnercircles).Thecontoursare
noticeablydistortedawayfromcircularsymmetrybytheboundaries Thisledtoaseverelackofisotropyinthenumericalresults,
ofthecube. asshownbydifferencesinthevaluesofthethreespreads(cid:12) ,
x
|             |         |           |           |                  |              |     | (cid:12) ,(cid:12) | thatpersistedundergridrefinement[8].Inthispaper |          |             |                |          |        |          |
| ----------- | ------- | --------- | --------- | ---------------- | ------------ | --- | ------------------ | ----------------------------------------------- | -------- | ----------- | -------------- | -------- | ------ | -------- |
|             |         |           |           |                  |              |     | y z                |                                                 |          |             |                |          |        |          |
|             |         |           |           |                  |              |     | we have            | shown                                           | that the | corrected   | scheme         | recovers |        | isotropy |
| vice        | versa.) | This part | of the    | solution is thus | not expected | to  |                    |                                                 |          |             |                |          |        |          |
|             |         |           |           |                  |              |     | to within          | a discretization                                |          | error       | that typically |          | scales | linearly |
| demonstrate |         | spherical | symmetry, | but first-order  | convergence  |     |                    |                                                 |          |             |                |          |        |          |
|             |         |           |           |                  |              |     | (or better)        | with                                            | the grid | resolution. | In             | other    | words, | we have  |
(cid:2)50
toward spherical symmetry of the solution in r is shown first-order convergence toward isotropy, at least for
confirmedbytheplotofnσ(r)inFig.8.
|     |     |     |     |     |     |     | times too | short | to allow | interaction |     | with the | boundary, | and |
| --- | --- | --- | --- | --- | --- | --- | --------- | ----- | -------- | ----------- | --- | -------- | --------- | --- |
thusestablishedtheviabilityofthequantumlatticeBoltzmann
|     |                                                |     | D. Effectofasmallerbox |     |           |      | approachinthreespatialdimensions. |     |     |     |     |     |     |     |
| --- | ---------------------------------------------- | --- | ---------------------- | --- | --------- | ---- | --------------------------------- | --- | --- | --- | --- | --- | --- | --- |
|     | Muchbetterconvergenceofthethreespreads(cid:12) |     |                        |     | ,(cid:12) | ,and |                                   |     |     |     |     |     |     |     |
x y
(cid:12) isobservedinasmallerboxofside(cid:15)=100.Figure9shows
z
|     |     |     |     |     |     | =0to |     |     | ACKNOWLEDGMENTS |     |     |     |     |     |
| --- | --- | --- | --- | --- | --- | ---- | --- | --- | --------------- | --- | --- | --- | --- | --- |
thatthethreespreadsremainincloseagreementfromt
t =400, and oscillate around the initial value (cid:12) =14. The P.J.D.’s and D.L.’s research was supported by the UK
0
plots of n(cid:12) (cid:20) in Fig. 10 demonstrate first-order convergence Engineering and Physical Sciences Research Council, Grant
90◦ (cid:2)50. No.EP/E054625/1.Someofthecomputationsmadeuseofthe
| toward |     | rotational | symmetry | for t | For later | times |     |     |     |     |     |     |     |     |
| ------ | --- | ---------- | -------- | ----- | --------- | ----- | --- | --- | --- | --- | --- | --- | --- | --- |
(cid:20) =O(n−1/2)
the convergence appears closer to (cid:12) rather than facilitiesoftheOxfordSupercomputingCentre.
[1]S. Chen and G. D. Doolen, Annu. Rev. Fluid Mech. 30, 329 [8]S.PalpacelliandS.Succi,Phys.Rev.E75,066704(2007).
(1998). [9]S. Palpacelli and S. Succi, Comput. Phys. Commun. 4, 980
| [2]S.Succi,TheLatticeBoltzmannEquation:ForFluidDynamics |     |     |     |     |     |     | (2008). |     |     |     |     |     |     |     |
| ------------------------------------------------------- | --- | --- | --- | --- | --- | --- | ------- | --- | --- | --- | --- | --- | --- | --- |
andBeyond(OxfordUniversityPress,Oxford,2001). [10]V. B. Berestetskii, E. M. Lifshitz, and L. P. Pitaevskii, Quan-
[3]D.Yu,R.Mei,L.-S.Luo,andW.Shyy,ProgressAerospaceSci. tumElectrodynamics(Butterworth-Heinemann,Oxford,1982),
|     | 39,329(2003). |     |     |     |     |     | 2nded. |     |     |     |     |     |     |     |
| --- | ------------- | --- | --- | --- | --- | --- | ------ | --- | --- | --- | --- | --- | --- | --- |
[4]C.K.AidunandJ.R.Clausen,Annu.Rev.FluidMech.42,439 [11]K.O.Friedrichs,Commun.PureAppl.Math.7,345(1954).
(2010). [12]S. Benzoni-Gavage and D. Serre, Multi-Dimensional Hyper-
[5]R. Benzi, S. Succi, and M. Vergassola, Phys. Rep. 222, 145 bolic Partial Differential Equations: First Order Systems and
|     | (1992). |     |     |     |     |     | Applications(OxfordUniversityPress,Oxford,2007). |     |     |     |     |     |     |     |
| --- | ------- | --- | --- | --- | --- | --- | ------------------------------------------------ | --- | --- | --- | --- | --- | --- | --- |
[6]S.SucciandR.Benzi,PhysicaD69,327(1993). [13]G. B. Whitham, Linear and Nonlinear Waves (Wiley Inter-
| [7]S.Succi,Phys.Rev.E53,1969(1996). |     |     |     |     |     |     | science,NewYork,1974). |     |     |     |     |     |     |     |
| ----------------------------------- | --- | --- | --- | --- | --- | --- | ---------------------- | --- | --- | --- | --- | --- | --- | --- |
046706-8

ISOTROPYOFTHREE-DIMENSIONALQUANTUMLATTICE... PHYSICALREVIEWE83,046706(2011)
[14]I.Bialynicki-Birula,Phys.Rev.D49,6920(1994). [29]Y.Zhang,Y.-W.Tan,H.L.Stormer,andP.Kim,Nature(London)
[15]D.A.Meyer,J.Stat.Phys.85,551(1996). 438,201(2005).
[16]D.A.Meyer,Phys.Rev.E55,5261(1997). [30]M.I.Katsnelson,K.S.Novoselov,andA.K.Geim,Nat.Phys.
[17]B.M.BoghosianandW.Taylor,PhysicaD120,30(1998). 2,620(2006).
[18]B. M. Boghosian and W. Taylor, Phys. Rev. E 57, 54 [31]A.H.CastroNeto,F.Guinea,N.M.R.Peres,K.S.Novoselov,
(1998). andA.K.Geim,Rev.Mod.Phys.81,109(2009).
[19]L. Vahala, G. Vahala, and J. Yepez, Phys. Lett. A 306, 227 [32]G. R. McNamara and G. Zanetti, Phys. Rev. Lett. 61, 2332
(2003). (1988).
[20]R.P.FeynmanandA.R.Hibbs,QuantumMechanicsandPath [33]F.J.HigueraandJ.Jime´nez,Europhys.Lett.9,663(1989).
Integrals(McGraw-Hill,NewYork,1965). [34]F. J. Higuera, S. Succi, and R. Benzi, Europhys. Lett. 9, 345
[21]U.Frisch,B.Hasslacher,andY.Pomeau,Phys.Rev.Lett.56, (1989).
1505(1986). [35]D. d’Humie`res, in Rarefied Gas Dynamics: Theory and Sim-
[22]S.Wolfram,J.Stat.Phys.45,471(1986). ulations, edited by B. D. Shizgal and D. P. Weaver, Prog.
[23]G.’tHooft,J.Stat.Phys.53,323(1988). Astronaut.Aeronaut,vol.159(AIAA,Washington,D.C.,1994),
[24]G.’tHooft,Int.J.Mod.Phys.A25,4385(2010). pp.450–458.
[25]S.Palpacelli,S.Succi,andR.Spigler,Phys.Rev.E76,036712 [36]P.LallemandandL.-S.Luo,Phys.Rev.E61,6546(2000).
(2007). [37]P.J.Dellar,J.Comput.Phys.190,351(2003).
[26]D. Lapitski and P. J. Dellar, Philos. Trans. R. Soc. London, [38]X.He,S.Chen,andG.D.Doolen,J.Comput.Phys.146,282
Ser.A,doi:10.1098/rsta.2011.0017. (1998).
[27]L.HaddadandL.Carr,PhysicaD238,1413(2009). [39]X.HeandL.-S.Luo,Phys.Rev.E55,R6333(1997).
[28]K.S.Novoselov, A.K.Geim,S.V.Morozov,D.Jiang,M.I. [40]X.HeandL.-S.Luo,Phys.Rev.E56,6811(1997).
Katsnelson,I.V.Grigorieva,S.V.Dubonos,andA.A.Firsov, [41]J. P. Boyd, Chebyshev and Fourier Spectral Methods, 2nd ed.
Nature(London)438,197(2005). (Dover,NewYork,2000).
046706-9

