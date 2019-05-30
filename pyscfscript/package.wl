(* ::Package:: *)

BeginPackage["Bredas`"]


Begin["`Private`"]


Sector[pointlist_]:=Transpose[{Drop[#[[1]],-1],Drop[#[[2]],1]}&@Table[pointlist,2]];
StructureExtractor[fchkfilestring_]:=MapAt[IntegerPart,MapAt[Partition[#,3]&,Flatten/@(Drop[Drop[#,1],-1]&/@Map[Select[#,NumberQ]&,ImportString[#,"Table"]&/@(StringTake[fchkfilestring,#]&/@Sector[Transpose[Join[StringPosition[fchkfilestring,"Nuclear charges"],StringPosition[fchkfilestring,"Current cartesian coordinates"],StringPosition[fchkfilestring,"Number of symbols"]]][[1]]]),{2}]),2],1];
AtomicNumberConvert[integerlist_]:=ElementData[#,"Abbreviation"]&/@integerlist;
StringCoordinateGenerator[fchkfilestring_]:=StringRiffle[StringRiffle[#,"   "]&/@Flatten/@Transpose[{AtomicNumberConvert[StructureExtractor[fchkfilestring][[1]]],Map[ToString[#,CForm]&,StructureExtractor[fchkfilestring][[2]],{2}]}],"\n"];
CoordinateGenerator[fchkfilestring_]:=Transpose[{AtomicNumberConvert[StructureExtractor[fchkfilestring][[1]]],StructureExtractor[fchkfilestring][[2]]}];
CoordinateToString[coordinatelist_]:=StringRiffle[StringRiffle[#,"   "]&/@Flatten/@coordinatelist,"\n"];
PyscfScriptGeneratorRHF[coordinatelist_,basis_]:=
"import pyscf, math
import numpy
from pyscf import gto, scf, mcscf,dft

mol = gto.Mole()
mol.build(
# verbose = 4,
symmetry = False,
atom =
'''
"<>CoordinateToString[coordinatelist]<>
"\n'''
,
basis='"<>basis<>"'

)
mf.kernel()

mf.analyze()
"

PyscfScriptGeneratorDFT[coordinatelist_,basis_,xc_]:=
"import pyscf, math
import numpy
from pyscf import gto, scf, mcscf,dft

mol = gto.Mole()
mol.build(
# verbose = 4,
symmetry = False,
atom =
'''
"<>CoordinateToString[coordinatelist]<>
"\n'''
,
basis='"<>basis<>"'

)

mf = dft.RKS(mol)
mf.xc ='"<>xc<>"'\n
mf.kernel()

mf.analyze()
"
PyscfScriptGenerator[coordinatelist_,basis_,method_]:=
If[method=="RHF",PyscfScriptGeneratorRHF[coordinatelist,basis],
If[method[[1]]=="DFT",PyscfScriptGeneratorDFT[coordinatelist,basis,method[[2]]]
]
];
BredasPyscfScriptGenerator[coordinatelist1_,coordinatelist2_,basis_,method_]:=
If[method=="RHF",
"import pyscf, math
import numpy
from pyscf import gto, scf, mcscf,dft

mol1 = gto.Mole()
mol1.build(
# verbose = 4,
symmetry = False,
atom =
'''
"<>CoordinateToString[coordinatelist1]<>
"\n'''
,
basis='"<>basis<>"'

)
homoNumber1 = mol1.tot_electrons()//2
mf1 = scf.RHF(mol1)
mf1.conv_tol = 1e-9
mf1.kernel()

mol2 = gto.Mole()
mol2.build(
# verbose = 4,
output = None,
symmetry = False,
atom = 
'''
"<>CoordinateToString[coordinatelist2]<>
"\n'''
,
basis='"<>basis<>"'

)
homoNumber2 = mol2.tot_electrons()//2
mf2 = scf.RHF(mol2)
mf2.conv_tol = 1e-9
mf2.kernel()

mol12 = gto.Mole()
mol12.build(
# verbose = 4,
output = None,
symmetry = False,
atom = 
'''
"<>CoordinateToString[Join[coordinatelist1,coordinatelist2]]<>
"\n'''
,
basis='"<>basis<>"'

)
homoNumber12 = mol12.tot_electrons()//2
mf12 = scf.RHF(mol12)
mf12.conv_tol = 1e-9
mf12.kernel()
mf12.analyze()

# Bredas
fock = mf12.get_fock()*27.2114
fockSmall_11 = fock[:len(fock)//2, :len(fock)//2]
fockSmall_12 = fock[:len(fock)//2, len(fock)//2:]
fockSmall_22 = fock[len(fock)//2:, len(fock)//2:]
overlap = mf12.get_ovlp()
overlap = overlap[:len(overlap)//2, len(overlap)//2:]
homo1 = mf1.mo_coeff[..., homoNumber1 - 1]
homo2 = mf2.mo_coeff[..., homoNumber2 - 1]
h11 = numpy.einsum(\"i,ij,j->\",homo1, fockSmall_11, homo1)
h22 = numpy.einsum(\"i,ij,j->\",homo2, fockSmall_22, homo2)
j12 = numpy.einsum(\"i,ij,j->\",homo1, fockSmall_12, homo2)
j21 = j12
s12 = numpy.einsum(\"i,ij,j->\",homo1, overlap, homo2)
s21 = s12
t12 = (j12 - 0.5 * (h11 + h22) * s12) / (1 - s12 * s12)
print(\"t12 from bredas results:  \", t12)
",
If[method[[1]]=="DFT",
"import pyscf, math
import numpy
from pyscf import gto, scf, mcscf,dft

mol1 = gto.Mole()
mol1.build(
# verbose = 4,
symmetry = False,
atom =
'''
"<>CoordinateToString[coordinatelist1]<>
"\n'''
,
basis='"<>basis<>"'

)
homoNumber1 = mol1.tot_electrons()//2
mf1 = dft.RKS(mol1)
mf1.xc = "<>method[[2]]<>"
mf1.conv_tol = 1e-9
mf1.kernel()

mol2 = gto.Mole()
mol2.build(
# verbose = 4,
output = None,
symmetry = False,
atom = 
'''
"<>CoordinateToString[coordinatelist2]<>
"\n'''
,
basis='"<>basis<>"'

)
homoNumber2 = mol2.tot_electrons()//2
mf2 = dft.RKS(mol2)
mf2.xc = "<>method[[2]]<>"
mf2.conv_tol = 1e-9
mf2.kernel()

mol12 = gto.Mole()
mol12.build(
# verbose = 4,
output = None,
symmetry = False,
atom = 
'''
"<>CoordinateToString[Join[coordinatelist1,coordinatelist2]]<>
"\n'''
,
basis='"<>basis<>"'

)
homoNumber12 = mol12.tot_electrons()//2
mf12 = dft.RKS(mol12)
mf1.xc = "<>method[[2]]<>"
mf12.conv_tol = 1e-9
mf12.kernel()
mf12.analyze()

# Bredas
fock = mf12.get_fock()*27.2114
fockSmall_11 = fock[:len(fock)//2, :len(fock)//2]
fockSmall_12 = fock[:len(fock)//2, len(fock)//2:]
fockSmall_22 = fock[len(fock)//2:, len(fock)//2:]
overlap = mf12.get_ovlp()
overlap = overlap[:len(overlap)//2, len(overlap)//2:]
homo1 = mf1.mo_coeff[..., homoNumber1 - 1]
homo2 = mf2.mo_coeff[..., homoNumber2 - 1]
h11 = numpy.einsum(\"i,ij,j->\",homo1, fockSmall_11, homo1)
h22 = numpy.einsum(\"i,ij,j->\",homo2, fockSmall_22, homo2)
j12 = numpy.einsum(\"i,ij,j->\",homo1, fockSmall_12, homo2)
j21 = j12
s12 = numpy.einsum(\"i,ij,j->\",homo1, overlap, homo2)
s21 = s12
t12 = (j12 - 0.5 * (h11 + h22) * s12) / (1 - s12 * s12)
print(\"t12 from bredas results:  \", t12)
"
]
];
Rotation[center1_,pointgroup_,xangle_,yangle_,zangle_]:=
Module[{translated},
translated =( #-center1&/@pointgroup);
translated={{1,0,0},{0,Cos[xangle],-Sin[xangle]},{0,Sin[xangle],Cos[xangle]}}.{{Cos[yangle],0,Sin[yangle]},{0,1,0},{-Sin[yangle],0,Cos[yangle]}}.{{Cos[zangle],-Sin[zangle],0},{Sin[zangle],Cos[zangle],0},{0,0,1}}.#&/@translated;
 #+center1&/@translated
];
Translation[pointgroup_,transformlist_]:=
#+transformlist&/@pointgroup;
pyscfPBSScriptGenerator[filename_]:=
"#PBS-N Walter
#PBS-l nodes=1:ppn=1
#PBS-q ljtest
#PBS-l walltime=99:60:00
#PBS-j oe
nprocs=`cat $PBS_NODEFILE|wc-l`
\[IndentingNewLine]cd $PBS_O_WORKDIR
cat $PBS_NODEFILE>NODEFILE
source/public/software/profile.d/mpi_openmpi-2.0.0-gnu.sh
source/public/software/compiler/intel2017/mkl/bin/mklvars.sh intel64
export OMP_NUM_THREADS=1
\[IndentingNewLine]python "<>filename<>" > "<>filename<>".out";


End[]


EndPackage[]
