#!/usr/bin/env perl
# -s -> SzW parameter (opt)
# -c -> CorrelMin parameter (opt)
# -i -> Inc parameter (opt)
# -r -> RegulBase parameter (opt)
# -m -> maskfile
# -o : name of output xml file (opt)
# arg1 img1
# arg2 img2

use warnings;
use strict;
use Getopt::Std;

# get & test options
my %opts;
exit(2) unless getopts('s:c:i:r:o:m:', \%opts);
exit(3) unless ($ARGV[0]);
exit(3) unless ($ARGV[1]);

my $input1=$ARGV[0];
my $input2=$ARGV[1];

my @content = &init_param();
if (defined($opts{"s"})) {
	@content = &change_size_window($opts{"s"}, @content);
}
if (defined($opts{"c"})) {
	@content = &change_correl_min($opts{"c"}, @content);
}
if (defined($opts{"i"})) {
	@content = &change_incertitude($opts{"i"}, @content);
}
if (defined($opts{"r"})) {
	@content = &change_regulbase($opts{"r"}, @content);
}
if (defined($opts{"m"})) {
	@content = &change_maskfile($opts{"m"}, @content);
}
@content = &change_inputs($input1, $input2, @content);

my $outputname = "ParamRun.xml";
$outputname = $opts{"o"} if (defined($opts{"o"}));

# write the param file
open(OUT, ">$outputname");
foreach my $line (@content) {
	print OUT "$line\n";
}
close(OUT);

sub change_maskfile()
{
	my ($mask, @params) = @_;
	my @result=();
	foreach my $item (@params) {
		$item =~ s/<Symb> Masq=(.*\.stop)/<Symb> Masq=$mask/ if ($item =~ m/<Symb> Masq=/);
		push(@result, $item)			;
	}
	return @result;
}

sub change_inputs()
{
	my ($im1, $im2, @params) = @_;
	my @result=();
	foreach my $item (@params) {
		$item =~ s/<Symb> Im1=(.*\.tif)/<Symb> Im1=$im1/ if ($item =~ m/<Symb> Im1=/);
		$item =~ s/<Symb> Im2=(.*\.tif)/<Symb> Im2=$im2/ if ($item =~ m/<Symb> Im2=/);
		push(@result, $item)			;
	}
	return @result;
}

sub change_size_window()
{
	my ($val, @params) = @_;
	my @result=();
	foreach my $item (@params) {
		if ($item =~ m/<Symb> SzW/) {
			$item =~ s/<Symb> SzW=(\d+)/<Symb> SzW=$val/;
		}
		push(@result, $item)			;
	}
	return @result;
}

sub change_correl_min()
{
	my ($val, @params) = @_;
	my @result=();
	foreach my $item (@params) {
		if ($item =~ m/<Symb> CorrelMin/) {
			$item =~ s/<Symb> CorrelMin=(\d*\.?\d*)/<Symb> CorrelMin=$val/;
		}
		push(@result, $item)			;
	}
	return @result;
}

sub change_incertitude()
{
	my ($val, @params) = @_;
	my @result=();
	foreach my $item (@params) {
		if ($item =~ m/<Symb> Inc=/) {
			$item =~ s/<Symb> Inc=(\d+)/<Symb> Inc=$val/;
		}
		push(@result, $item)			;
	}
	return @result;
}

sub change_regulbase()
{
	my ($val, @params) = @_;
	my @result=();
	foreach my $item (@params) {
		if ($item =~ m/<Symb> RegulBase=/) {
			$item =~ s/<Symb> RegulBase=(\d*\.?\d*)/<Symb> RegulBase=$val/;
		}
		push(@result, $item)			;
	}
	return @result;
}

# DGN : modifier aussi le fichier param pour le mask

sub init_param()
{
	my $text=<< 'EOF';
<ParamMICMAC>
<DicoLoc>

<Symb> Im1=orthoimg_phr1a_p_201208071034251_sen_571784101-001_r1c1_cropped.tif  </Symb>
<Symb> Im2=orthoimg_phr1a_p_201210051030181_sen_579436201-002_r1c1_cropped.tif  </Symb>
<Symb> Masq=mask.stop  </Symb>

<Symb> Dir=MEC/</Symb>
<Symb> Pyr=PyramSat/</Symb>
<Symb> Inc=5</Symb>
<Symb> Pas=0.2 </Symb>
<Symb> Teta0=0 </Symb>
<Symb> UseMasq=true </Symb>
<Symb> UseTeta=false </Symb>
<Symb> RegulBase=0.2 </Symb>
<Symb> UseDequant=true </Symb>

<Symb> SzW=3 </Symb>
<Symb> CorrelMin=0.2 </Symb>
<Symb> GammaCorrel=2 </Symb>

<Symb> PdsF=0.1 </Symb>

<Symb> NbDir=7 </Symb>
<eSymb> P0= / ${PdsF} +   ${PdsF} / 0 ${NbDir} </eSymb>
<eSymb> P1= / ${PdsF} +   ${PdsF} / 1 ${NbDir} </eSymb>
<eSymb> P2= / ${PdsF} +   ${PdsF} / 2 ${NbDir} </eSymb>
<eSymb> P3= / ${PdsF} +   ${PdsF} / 3 ${NbDir} </eSymb>
<eSymb> P4= / ${PdsF} +   ${PdsF} / 4 ${NbDir} </eSymb>
<eSymb> P5= / ${PdsF} +   ${PdsF} / 5 ${NbDir} </eSymb>
<eSymb> P6= / ${PdsF} +   ${PdsF} / 6 ${NbDir} </eSymb>
<eSymb> P7= / ${PdsF} +   ${PdsF} / 7 ${NbDir} </eSymb>

<Symb>  VPds=[${P0},${P1},${P2},${P3},${P4},${P5},${P6},${P7},${P6},${P7},${P4},${P3},${P2},${P1}] </Symb>
<eSymb> NbDirTot=* 2 ${NbDir} </eSymb>


<eSymb> Regul=* ${RegulBase}  ? ${UseTeta} 3 1 </eSymb>
</DicoLoc>
<!-- *************************************************************
Parametres lies au terrain "physique", independamment de la prise de vue
  *************************************************************-->
<Section_Terrain>
<IntervParalaxe>
<!-- Incertitude en parallaxe -->
<!--Px1IncCalc et Px2IncCalc permettent de definir les deux
nappes encadrantes utilisees pour definir la zone de recherche
au premier niveau de la pyramide-->
<Px1IncCalc>  ${Inc}  </Px1IncCalc>
<Px2IncCalc>  ${Inc}   </Px2IncCalc>

<!--Px1Moy et Px2Moy fixent la valeur moyenne de la parallaxe;
leurs valeurs ont donc une influence directe sur la zone de
recherche exploree lors du premier niveau de la pyramide de
resolution. Accessoirement, elles ont une influence sur le
formatage du resultat (exprime en relatif par rapport a cette
valeur moyenne).-->
<Px1Moy >  0.0     </Px1Moy>
<Px2Moy >  0.0   </Px2Moy>

</IntervParalaxe>

<Planimetrie>
<!-- image de Masque utilisee pour designer l'emprise fine de
la correlation; il doit etre superposable au MNT de resolution
1. Si le fichier n'existe pas, il en sera cree un correspondant
aux point du terrain qui sont vus d'au moins deux images (pour
la parallaxe moyenne) -->
<#WHEN VTEST=${UseMasq}>
<MasqueTerrain>
<MT_Image> ${Masq}.tif </MT_Image>
<MT_Xml>   ${Masq}.xml </MT_Xml>
</MasqueTerrain>
</#WHEN>
</Planimetrie>
</Section_Terrain>

<!-- *************************************************************
Parametres lies a la prise de vue, independamment de son exploitation
par le correlateur
************************************************************* -->
<Section_PriseDeVue>
<GeomImages> eGeomImage_Epip </GeomImages>
<Images>
<Im1> ${Im1} </Im1>
<Im2> ${Im2} </Im2>
</Images>
<!--
<FiltreImageIn>
  <TypeFiltrage>  eFiltrageDeriche  </TypeFiltrage>
  <SzFiltrage>   1.0          </SzFiltrage>
</FiltreImageIn>
-->
</Section_PriseDeVue>

<!--  *************************************************************
Parametres reglant le comportement de l'algo de mise en correspondance

La premiere etape doit obligatoirement avoir le champs
resolution a -1. Elle donne les valeurs par defaut et ne
sera pas executee.

L'ordre des resolutions : les plus basses aux plus grandes.
************************************************************* -->
<Section_MEC>
<ChantierFullImage1> true </ChantierFullImage1>

<ClipMecIsProp> false </ClipMecIsProp>

<EtapeMEC><!-- Etape de Mise En Correspondance -->
<DeZoom > -1 </DeZoom> <!-- le seul fils obligatoire a toutes les etapes-->
<SzW> ${SzW}   </SzW> <!-- la taille de la fenetre de correlation [-4,4]x[-4,4]-->


<CorrelMin>  ${CorrelMin} </CorrelMin>
<GammaCorrel>  ${GammaCorrel} </GammaCorrel>
<DynamiqueCorrel> eCoeffGamma </DynamiqueCorrel>


<AlgoRegul> eAlgo2PrgDyn </AlgoRegul>
<ModulationProgDyn Portee="Globale">
<EtapeProgDyn>
<ModeAgreg>    ePrgDAgrSomme   </ModeAgreg>
<#WHEN VTEST=${UseTeta}>
<Px1MultRegul> ${VPds} </Px1MultRegul>
<Px2MultRegul> ${VPds} </Px2MultRegul>
</#WHEN>
<NbDir>      ${NbDirTot}     </NbDir>
<Teta0>    ${Teta0}    </Teta0>
</EtapeProgDyn>
<Px1PenteMax > 2 </Px1PenteMax>
<Px2PenteMax > 2 </Px2PenteMax>
</ModulationProgDyn>
<Px1Regul>  ${Regul}    </Px1Regul>
<Px2Regul>  ${Regul}   </Px2Regul>

<GenImagesCorrel> true </GenImagesCorrel>

<ModeInterpolation> eInterpolSinCard </ModeInterpolation>
<SzSinCard>  5.0 </SzSinCard>
<SzAppodSinCard>  5.0 </SzAppodSinCard>

<Px1DilatAlti>  2    </Px1DilatAlti>
<Px1DilatPlani> 2    </Px1DilatPlani>
<Px2DilatAlti>  2    </Px2DilatAlti>
<Px2DilatPlani> 2    </Px2DilatPlani>

<Px1Pas>        ${Pas}  </Px1Pas>
<Px2Pas>        ${Pas}  </Px2Pas>
<SsResolOptim> 2 </SsResolOptim>

</EtapeMEC>

<EtapeMEC>
<DeZoom > 2</DeZoom>
<Px1Pas>   1  </Px1Pas>
<Px2Pas>   1 </Px2Pas>
</EtapeMEC>

<EtapeMEC>
<DeZoom > 1</DeZoom>
<Px1Pas>   0.8  </Px1Pas>
<Px2Pas>   0.8 </Px2Pas>
</EtapeMEC>
<EtapeMEC>
<DeZoom > 1</DeZoom>
<Px1Pas>   0.4  </Px1Pas>
<Px2Pas>   0.4 </Px2Pas>
</EtapeMEC>

<EtapeMEC>
<DeZoom > 1</DeZoom>
<Px1Pas>   0.2  </Px1Pas>
<Px2Pas>   0.2 </Px2Pas>
</EtapeMEC>

<EtapeMEC>
<DeZoom > 1</DeZoom>
<Px1Pas>   0.1  </Px1Pas>
<Px2Pas>   0.1 </Px2Pas>
</EtapeMEC>


<#WHEN VTEST=${UseDequant}>
<EtapeMEC>
<DeZoom >  1  </DeZoom>
<Px1Pas>   1.0     </Px1Pas>
<Px2Pas>   1.0     </Px2Pas>
<AlgoRegul> eAlgoDequant </AlgoRegul>
</EtapeMEC>
</#WHEN>

</Section_MEC>


<Section_Results >
<GeomMNT> eGeomPxBiDim </GeomMNT>
</Section_Results>

<Section_WorkSpace >
<WorkDir >  ThisDir </WorkDir>
<TmpMEC>    ${Dir} </TmpMEC>
<TmpResult> ${Dir} </TmpResult>
<TmpPyr>  ${Pyr}  </TmpPyr>
<ByProcess>  ${MMNbProc} </ByProcess>

<NbCelluleMax> 8e7 </NbCelluleMax>

<SzRecouvrtDalles> 50 </SzRecouvrtDalles>
<SzDalleMin> 500 </SzDalleMin>

<DefTileFile>100000</DefTileFile>

</Section_WorkSpace>

<Section_Vrac> </Section_Vrac>

</ParamMICMAC>
EOF
return split("\n", $text);
}
