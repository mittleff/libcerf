/* Library libcerf:
 *   Compute complex error functions, based on a new implementation of
 *   Faddeeva's w_of_z. Also provide Dawson and Voigt functions.
 *
 * File im_w_of_x.c:
 *   Compute scaled Dawson integral im_w_of_x(x) = 2*dawson(x)/sqrt(pi),
 *   equivalent to the imaginary part of the Faddeeva function w(x) for real x.
 *
 * Copyright:
 *   (C) 2012 Massachusetts Institute of Technology
 *   (C) 2013, 2024 Forschungszentrum Jülich GmbH
 *
 * Licence:
 *   Permission is hereby granted, free of charge, to any person obtaining
 *   a copy of this software and associated documentation files (the
 *   "Software"), to deal in the Software without restriction, including
 *   without limitation the rights to use, copy, modify, merge, publish,
 *   distribute, sublicense, and/or sell copies of the Software, and to
 *   permit persons to whom the Software is furnished to do so, subject to
 *   the following conditions:
 *
 *   The above copyright notice and this permission notice shall be
 *   included in all copies or substantial portions of the Software.
 *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *   LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *   OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *   WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * Authors:
 *   Steven G. Johnson, Massachusetts Institute of Technology, 2012
 *   Joachim Wuttke, Forschungszentrum Jülich, 2013, 2024
 *
 * Website:
 *   http://apps.jcns.fz-juelich.de/libcerf
 *
 * Revision history:
 *   The numerics in this file has been almost entirely reworked by Joachim Wuttke
 *   in September 2024. For large |x|, continuous fractions were replaced with
 *   asymptotic expansions. The small |x| range where Taylor series is used was
 *   expanded. In the middle range, the subdivision in subranges was redone. To avoid
 *   a division, the transformation x -> 1/(1+x) was dropped. Chebyshev coefficients were
 *   recomputed using a Python/mpmath script included in the source repository.
 *   Thereby, an accuracy better than 4 epsilon was achieved.
 *
 *   See also ../CHANGELOG
 *
 * Manual page:
 *   man 3 im_w_of_x
 */

#include "cerf.h"
#include <math.h>
#include "defs.h" // defines _cerf_cmplx, NaN, C, cexp, ...

// for analysing the algorithm:
IMPORT extern int faddeeva_algorithm;
IMPORT extern int faddeeva_nofterms;

// The follwowing four functions compute Im w(x) by Chebyshev interpolation.
// Each function addresses a certain range aCheb<n> to bCheb<n>, which is further
// divided in equidistant subranges.
// Each subrange has its own Chebyshev interpolant.
// Maximum degree is 7. A few interpolants have lower degree.
// Coefficients are computed, and corresponding code sections are written,
// by script devtool/chebcoef_imwofx.py.

//--- The following code is generated by devtool/pro_imwofx_chebcoeffs.py
// clang-format off
static const double aCheb1 = 0.94;
static const double bCheb1 = 1.8;
static double chebInterpolant1(double x)
{
    static const int nSubranges = 17;
    static const double invSubwidth = nSubranges / (bCheb1 - aCheb1);
    static const double Coeffs[17 * 8] = {
        +6.0949556181098219e-01, -1.2218524364895115e-03, -3.6011755028194135e-04, +6.3829711316041099e-06, +3.7276086790846396e-08, -1.9975675476173714e-09, +8.3016540844675932e-12, +3.0713420354437454e-13,
        +6.0566298357129666e-01, -2.5846922253612947e-03, -3.2108270840081591e-04, +6.6027725335981830e-06, +1.7881221898088104e-08, -1.8735506315618358e-09, +1.2226910246724807e-11, +2.5262764777535990e-13,
        +5.9926231742749836e-01, -3.7893650152175445e-03, -2.8118370854837986e-04, +6.6729726984491224e-06, -5.3559569644387983e-11, -1.7071488777309801e-09, +1.5354258346154200e-11, +1.9365454397568875e-13,
        +5.9061218186564635e-01, -4.8341594306312054e-03, -2.4127991960043132e-04, +6.6068178043788285e-06, -1.6153343066298003e-08, -1.5082314544734909e-09, +1.7642437968554883e-11, +1.3322624796460099e-13,
        +5.8003129229713890e-01, -5.7206314163884475e-03, -2.0214303284227529e-04, +6.4201536463908699e-06, -3.0143550302103701e-08, -1.2869150355184713e-09, +1.9090844192470047e-11, +7.4116681870335256e-14,
        +5.6783229631485355e-01, -6.4532255545921560e-03, -1.6444388255206023e-04, +6.1306194056093523e-06, -4.1850049695314986e-08, -1.0530993330256279e-09, +1.9735312249754871e-11, +1.8740264806840661e-14,
        +5.5431641259571118e-01, -7.0388533061064512e-03, -1.2874407055255406e-04, +5.7568583738857089e-06, -5.1194907210787752e-08, -8.1606259714458408e-10, +1.9642425704043431e-11, -3.0936418461949013e-14,
        +5.3976894058849767e-01, -7.4864470543107454e-03, -9.5492193497805167e-05, +5.3177778948304038e-06, -5.8188456225543424e-08, -5.8413189864467234e-10, +1.8902836406879722e-11, -7.3447380403814247e-14,
        +5.2445567142155702e-01, -7.8065076598710596e-03, -6.5024295448689825e-05, +4.8318848395420977e-06, -6.2918483585496225e-08, -3.6443706860123418e-10, +1.7624091716200739e-11, -1.0784517980283280e-13,
        +5.0862019675488168e-01, -8.0106614357575833e-03, -3.7568030105364236e-05, +4.3167162816684756e-06, -6.5537388387555946e-08, -1.6275134889470281e-10, +1.5923437744172397e-11, -1.3368704983794371e-13,
        +4.9248208268873944e-01, -8.1112401807986994e-03, -1.3249919767625813e-05, +3.7883781351561514e-06, -6.6248173391030673e-08, +1.6583516258276341e-11, +1.3921012572392245e-11, -1.5099395860919260e-13,
        +4.7623585110455552e-01, -8.1208953324435115e-03, +7.8949581258765582e-06, +3.2611977930450297e-06, -6.5290082845709992e-08, +1.7065454954643955e-10, +1.1733767715789527e-11, -1.6018856222050900e-13,
        +4.6005069140455107e-01, -8.0522545757831902e-03, +2.5911543349134501e-05, +2.7474906181464351e-06, -6.2924613928445190e-08, +2.9790751045822780e-10, +9.4703676062867128e-12, -1.6201864065092479e-13,
        +4.4407081167592916e-01, -7.9176265241911563e-03, +4.0912293135642113e-05, +2.2574347687536074e-06, -5.9422512455505466e-08, +3.9802234362997383e-10, +7.2272234404502438e-12, -1.5747268692261296e-13,
        +4.2841632969721349e-01, -7.7287564956177333e-03, +5.3064232623456473e-05, +1.7990444896096672e-06, -5.5052228230561947e-08, +4.7174820686801724e-10, +5.0857266595146390e-12, -1.4769391779130460e-13,
        +4.1318460055937523e-01, -7.4966340518745692e-03, +6.2576108313604849e-05, +1.3782287665230847e-06, -5.0070162172258202e-08, +5.2071245532854056e-10, +3.1106647107095048e-12, -1.3389818611743943e-13,
        +3.9845187844130409e-01, -7.2313509236870835e-03, +6.9686111918225918e-05, +9.9892015517446934e-07, -4.4712896114014378e-08, +5.4721731638755746e-10, +1.3497321913287904e-12, -1.1730021524601828e-13,
    };
    const int s = (int)((x-aCheb1)*invSubwidth); // index of subrange
    faddeeva_nofterms = s;
    const double center = ((nSubranges-0.5)-s)*(aCheb1/nSubranges) + (s+0.5) * (bCheb1/nSubranges);
    const double t = 2 * invSubwidth * (x-center); // coord in subrange, between -1 and +1
    const double *const c = Coeffs + s * 8;
    return (((((( c[7] * t + c[6] ) * t + c[5] ) * t + c[4] ) * t + c[3] ) * t + c[2] ) * t + c[1] ) * t + c[0];
}
// clang-format on
//--- End of autogenerated code

//--- The following code is generated by devtool/pro_imwofx_chebcoeffs.py
// clang-format off
static const double aCheb2 = 1.8;
static const double bCheb2 = 3.4;
static double chebInterpolant2(double x)
{
    static const int nSubranges = 29;
    static const double invSubwidth = nSubranges / (bCheb2 - aCheb2);
    static const double Coeffs[29 * 8] = {
        +3.8364676093775235e-01, -7.5563004172479300e-03, +8.9005006767709527e-05, +8.4202931439963118e-07, -5.5092302135327613e-08, +8.5470378602471896e-10, -3.8388264364975346e-13, -1.8023002806963244e-13,
        +3.6889606219123722e-01, -7.1918707686213935e-03, +9.2803133616510501e-05, +4.3532018673990900e-07, -4.6616406788514640e-08, +8.3595700899242039e-10, -2.6429763544498246e-12, -1.4247732676492429e-13,
        +3.5488629645102732e-01, -6.8168598100600318e-03, +9.4362410620824359e-05, +9.5328019775006805e-08, -3.8452969847742351e-08, +7.9326301870812413e-10, -4.3772821111537334e-12, -1.0556897312192287e-13,
        +3.4163019894111624e-01, -6.4394341515559275e-03, +9.4073849995317247e-05, -1.8132141638659784e-07, -3.0810324615985395e-08, +7.3280144802845962e-10, -5.6107052697245213e-12, -7.1156981144033598e-14,
        +3.2912570558321086e-01, -6.0662430226069743e-03, +9.2303706137222453e-05, -3.9942654075334551e-07, -2.3836876615458587e-08, +6.6033831726912895e-10, -6.3871032721678830e-12, -4.0439905455313053e-14,
        +3.1735887827732739e-01, -5.7025325127839028e-03, +8.9386331185111468e-05, -5.6475002270590238e-07, -1.7626310338029748e-08, +5.8102446029029153e-10, -6.7640204466328171e-12, -1.4179503337893612e-14,
        +3.0630657671363820e-01, -5.3522830523449576e-03, +8.5619650612314615e-05, -6.8360756968616449e-07, -1.2224442210239298e-08, +4.9926701368294133e-10, -6.8068409032574928e-12, +7.2574943155731007e-15,
        +2.9593884030177009e-01, -5.0183602844018228e-03, +8.1262932715733135e-05, -7.6251576873641979e-07, -7.6370207633687624e-09, +4.1866705276703801e-10, -6.5836555149996272e-12, +2.3840441576894036e-14,
        +2.8622096212451681e-01, -4.7026708868647318e-03, +7.6536480053504571e-05, -8.0790403680032814e-07, -3.8378597048053347e-09, +3.4201356948105177e-10, -6.1610292268153836e-12, +3.5806299100209220e-14,
        +2.7711525218774924e-01, -4.4063164319987948e-03, +7.1622854474783152e-05, -8.2589123805278115e-07, -7.7679671786116169e-10, +2.7132237400874370e-10, -5.6007535467047879e-12, +4.3588393842971556e-14,
        +2.6858249951249646e-01, -4.1297399162469883e-03, +6.6669255278085110e-05, -8.2212393575598436e-07, +1.6129117360412512e-09, +2.0790802233051249e-10, -4.9575824771634806e-12, +4.7747967830420406e-14,
        +2.6058315185812042e-01, -3.8728610669172274e-03, +6.1790696635298915e-05, -8.0167060098263017e-07, +3.4080319336337316e-09, +1.5247732247118900e-10, -4.2778821738339047e-12, +4.8912659064043210e-14,
        +2.5307823828627762e-01, -3.6351978718369795e-03, +5.7073670123731178e-05, -7.6896436652917005e-07, +4.6897938256662404e-09, +1.0523411856135510e-10, -3.5990766904023895e-12, +4.7724625903829296e-14,
        +2.4603006368810582e-01, -3.4159729413503468e-03, +5.2580025798917278e-05, -7.2778600135519626e-07, +5.5393992163639113e-09, +6.5986655013606382e-11, -2.9497436887684395e-12, +4.4799697653802068e-14,
        +2.3940270617946699e-01, -3.2142042769025014e-03, +4.8350856156766443e-05, -6.8127856215950795e-07, +6.0345871125851087e-09, +3.4250661068961742e-11, -2.3502028462471913e-12, +4.0697798124048073e-14,
        +2.3316234832596092e-01, -3.0287807813822710e-03, +4.4410217841221410e-05, -6.3198550847402335e-07, +6.2471870157891934e-09, +9.3431785055713616e-12, -1.8134427756019494e-12, +3.5903997422327878e-14,
        +2.2727747189295663e-01, -2.8585234109896928e-03, +4.0768573212214456e-05, -5.8190479631835905e-07, +6.2415501839352164e-09, -9.5360671314315043e-12, -1.3462460061320417e-12, +3.0818917993092634e-14,
        +2.2171894360284244e-01, -2.7022332538654561e-03, +3.7425875935937256e-05, -5.3255245151033310e-07, +6.0737249361522400e-09, -2.3237384665251090e-11, -9.5039232401168222e-13, +2.5756844142965419e-14,
        +2.1646001655767227e-01, -2.5587280504783027e-03, +3.4374260435587558e-05, -4.8503024225335171e-07, +5.7912324504166046e-09, -3.2608265415708042e-11, -6.2384519426376223e-13, +2.0949730007930191e-14,
        +2.1147626883548812e-01, -2.4268687714526054e-03, +3.1600323872236452e-05, -4.4009321841397348e-07, +5.4333013079524274e-09, -3.8454310648260135e-11, -3.6185125680271971e-13, +1.6555326760264394e-14,
        +2.0674549752353047e-01, -2.3055778675199179e-03, +2.9087011386738887e-05, -3.9821398263036415e-07, +5.0314292548038250e-09, -4.1512593598001052e-11, -1.5790681463958670e-13, +1.2667808105028092e-14,
        +2.0224758328812437e-01, -2.1938507298295130e-03, +2.6815131103625277e-05, -3.5964155035310805e-07, +4.6101564234793974e-09, -4.2435688491941038e-11, -4.5662939639306655e-15, +9.3295165293331207e-15,
        +1.9796433762591212e-01, -2.0907617706930142e-03, +2.4764535623641179e-05, -3.2445350801608229e-07, +4.1879530662256266e-09, -4.1784407052042767e-11, +1.0591503408622740e-13, +6.5427388555151188e-15,
        +1.9387934227669795e-01, -1.9954663753865270e-03, +2.2915012307556406e-05, -2.9260087509215270e-07, +3.7781446221359353e-09, -4.0027274293137758e-11, +1.8109856429665950e-13, +4.2807133940519654e-15,
        +1.8997778794971576e-01, -1.9071998016058729e-03, +2.1246926522690617e-05, -2.6394561578462018e-07, +3.3898160862022665e-09, -3.7544901018394480e-11, +2.2802538648323201e-13, +2.4973461076570600e-15,
        +1.8624631753818838e-01, -1.8252739275464163e-03, +1.9741661099878027e-05, -2.3829113928009749e-07, +3.0286550609488313e-09, -3.4637624262267366e-11, +2.5300452379213694e-13, +1.1353522002757444e-15,
        +1.8267287736478499e-01, -1.7490725818287314e-03, +1.8381892365427963e-05, -2.1540639068219380e-07, +2.6977078738622548e-09, -3.1535049003856833e-11, +2.6150351011384334e-13, +1.3273294446669187e-16,
        +1.7924657869041394e-01, -1.6780460349733889e-03, +1.7151738988917384e-05, -1.9504428870734482e-07, +2.3980354231441582e-09, -2.8406402316259676e-11, +2.5811824393262992e-13, -5.7235749687567946e-16,
        +1.7595757069811629e-01, -1.6117050965818618e-03, +1.6036815112682094e-05, -1.7695533301935693e-07, +2.1292649412550595e-09, -2.5370880350890267e-11, +2.4660056368143051e-13, -1.0382018539783940e-15,
    };
    const int s = (int)((x-aCheb2)*invSubwidth); // index of subrange
    faddeeva_nofterms = s;
    const double center = ((nSubranges-0.5)-s)*(aCheb2/nSubranges) + (s+0.5) * (bCheb2/nSubranges);
    const double t = 2 * invSubwidth * (x-center); // coord in subrange, between -1 and +1
    const double *const c = Coeffs + s * 8;
    return (((((( c[7] * t + c[6] ) * t + c[5] ) * t + c[4] ) * t + c[3] ) * t + c[2] ) * t + c[1] ) * t + c[0];
}
// clang-format on
//--- End of autogenerated code

//--- The following code is generated by devtool/pro_imwofx_chebcoeffs.py
// clang-format off
static const double aCheb3 = 3.4;
static const double bCheb3 = 5.84;
static double chebInterpolant3(double x)
{
    static const int nSubranges = 25;
    static const double invSubwidth = nSubranges / (bCheb3 - aCheb3);
    static const double Coeffs[25 * 8] = {
        +1.7161408356739341e-01, -2.7008914230842759e-03, +4.5875272617772727e-05, -8.5924238090686965e-07, +1.7681261530499215e-08, -3.7181908828783327e-10, +6.8217936660564007e-12, -7.4965301984432387e-14,
        +1.6638919930150747e-01, -2.5271639106587210e-03, +4.1116008942132608e-05, -7.3161660350234390e-07, +1.4350710713939903e-08, -2.9651805278816908e-10, +5.7139577926303676e-12, -8.0686321819319615e-14,
        +1.6149370306131297e-01, -2.3710427117884137e-03, +3.7048322299551829e-05, -6.2780230908532815e-07, +1.1705984338587409e-08, -2.3462763454686850e-10, +4.6165137772956853e-12, -7.4815029962936313e-14,
        +1.5689496858209642e-01, -2.2300263757439346e-03, +3.3544740369234361e-05, -5.4284187687301069e-07, +9.6163766752554768e-09, -1.8524392500750502e-10, +3.6424557586026127e-12, -6.3893876800273352e-14,
        +1.5256490021634114e-01, -2.1020679410112967e-03, +3.0504494898553563e-05, -4.7277249004978288e-07, +7.9653496170582967e-09, -1.4658106742024071e-10, +2.8326131604158508e-12, -5.1821080837774606e-14,
        +1.4847912306375408e-01, -1.9854795454519632e-03, +2.7847947848537949e-05, -4.1448765630008114e-07, +6.6557128557508091e-09, -1.1663620462728580e-10, +2.1871261939931000e-12, -4.0633980148674219e-14,
        +1.4461634275690971e-01, -1.8788575518622918e-03, +2.5511926578099282e-05, -3.6557928073455598e-07, +5.6098285463163664e-09, -9.3540337427026655e-11, +1.6866053058783571e-12, -3.1182835165937800e-14,
        +1.4095783759324509e-01, -1.7810244554231607e-03, +2.3445988072700898e-05, -3.2418911294717437e-07, +4.7673982823906774e-09, -7.5707590101978698e-11, +1.3050048079277167e-12, -2.3637362416946568e-14,
        +1.3748705305806966e-01, -1.6909840321055359e-03, +2.1609512200515772e-05, -2.8888207691398819e-07, +4.0823986725031353e-09, -6.1868040802775031e-11, +1.0165806361213557e-12, -1.7828922899599697e-14,
        +1.3418927538747077e-01, -1.6078866935093394e-03, +1.9969480268816651e-05, -2.5854449435837719e-07, +3.5200193188521838e-09, -5.1042012643251773e-11, +7.9903376243532131e-13, -1.3459392370375849e-14,
        +1.3105136630198602e-01, -1.5310026015704868e-03, +1.8498793431581454e-05, -2.3230537993364724e-07, +3.0539946891880243e-09, -4.2491406691777283e-11, +6.3443892448983651e-13, -1.0214819677748064e-14,
        +1.2806154537304323e-01, -1.4597006465359400e-03, +1.7175003336885761e-05, -2.0947704315707907e-07, +2.6644509777033972e-09, -3.5667606399642431e-11, +5.0904496039331594e-13, -7.8179425811725073e-15,
        +1.2520920979842234e-01, -1.3934318543652032e-03, +1.5979351575430126e-05, -1.8951088901763052e-07, +2.3362498704140214e-09, -3.0165103774623529e-11, +4.1261176597157142e-13, -6.0453888948357888e-15,
        +1.2248478384924304e-01, -1.3317161553686233e-03, +1.4896038116946890e-05, -1.7196472711765053e-07, +2.0577528428548233e-09, -2.5683956285724457e-11, +3.3764764798458220e-13, -4.7269358771577927e-15,
        +1.1987959209633006e-01, -1.2741317234944364e-03, +1.3911659056685170e-05, -1.5647858537243931e-07, +1.8199159575657498e-09, -2.2001081876569501e-11, +2.7871799114187010e-13, -3.7373324401993482e-15,
        +1.1738575188887287e-01, -1.2203063011932810e-03, +1.3014769891087520e-05, -1.4275672287057245e-07, +1.6156327102047936e-09, -1.8949088840316932e-11, +2.3188524609490064e-13, -2.9863509058008563e-15,
        +1.1499608157048448e-01, -1.1699100747727159e-03, +1.2195542510484694e-05, -1.3055413721655219e-07, +1.4392577406576226e-09, -1.6401001628781390e-11, +1.9428401772874779e-13, -2.4096933539004175e-15,
        +1.1270402167330533e-01, -1.1226497739659578e-03, +1.1445492842733829e-05, -1.1966633426040040e-07, +1.2862602352107780e-09, -1.4259365831483992e-11, +1.6381152330898109e-13, -1.9616641194741592e-15,
        +1.1050356690017070e-01, -1.0782637484385771e-03, +1.0757262362208980e-05, -1.0992148040657755e-07, +1.1529696407289877e-09, -1.2448501767583121e-11, +1.3890762730046451e-13, -1.6096719772284504e-15,
        +1.0838920713874818e-01, -1.0365178316442173e-03, +1.0124441152553322e-05, -1.0117431162339750e-07, +1.0363870807486262e-09, -1.0908979331361121e-11, +1.1840059172919074e-13, -1.3302820507773223e-15,
        +1.0635587608609252e-01, -9.9720184495799124e-04, +9.5414233968959639e-06, -9.3301352707972824e-08, +9.3404379238580720e-10, -9.5936461394937826e-12, +1.0139937636843333e-13, -1.1064550617647496e-15,
        +1.0439890632299620e-01, -9.6012662671194618e-04, +9.0032884486391536e-06, -8.6197125956787828e-08, +8.4389353333120364e-10, -8.4647409677386757e-12, +8.7218244222295671e-14, -9.2564415970450200e-16,
        +1.0251398988324899e-01, -9.2512149467712749e-04, +8.5057022820420435e-06, -7.9771116177729420e-08, +7.6422982193259643e-10, -7.4917694450977838e-12, +7.5323572380756642e-14, -7.7848917604186556e-16,
        +1.0069714352675627e-01, -8.9203206881147373e-04, +8.0448353245116107e-06, -7.3945320631186361e-08, +6.9362156598052570e-10, -6.6499202412507145e-12, +6.5295908129770498e-14, -6.5791846059646437e-16,
        +9.8944678057109631e-02, -8.6071839521791071e-04, +7.6172935625133325e-06, -6.8652256185459519e-08, +6.3086248510756785e-10, -5.9188691896812988e-12, +5.6802546951707182e-14, -5.5852612967356663e-16,
    };
    const int s = (int)((x-aCheb3)*invSubwidth); // index of subrange
    faddeeva_nofterms = s;
    const double center = ((nSubranges-0.5)-s)*(aCheb3/nSubranges) + (s+0.5) * (bCheb3/nSubranges);
    const double t = 2 * invSubwidth * (x-center); // coord in subrange, between -1 and +1
    const double *const c = Coeffs + s * 8;
    return (((((( c[7] * t + c[6] ) * t + c[5] ) * t + c[4] ) * t + c[3] ) * t + c[2] ) * t + c[1] ) * t + c[0];
}
// clang-format on
//--- End of autogenerated code

//--- The following code is generated by devtool/pro_imwofx_chebcoeffs.py
// clang-format off
static const double aCheb4 = 5.84;
static const double bCheb4 = 10.9;
static double chebInterpolant4(double x)
{
    static const int nSubranges = 28;
    static const double invSubwidth = nSubranges / (bCheb4 - aCheb4);
    static const double Coeffs[28 * 8] = {
        +9.6550658873895523e-02, -1.5162481471508466e-03, +2.4202282965850871e-05, -3.9301328540687725e-07, +6.4993887579417626e-09, -1.0959030764084057e-10, +1.8879626106647567e-12, -3.3248289507691891e-14,
        +9.3611928205163356e-02, -1.4239556135108731e-03, +2.1991852898091836e-05, -3.4511748758385870e-07, +5.5080183579141205e-09, -8.9492115517822901e-11, +1.4828416956077880e-12, -2.5064417971959295e-14,
        +9.0849308806134593e-02, -1.3399602407910014e-03, +2.0046520640123273e-05, -3.0440919545060972e-07, +4.6954581363613192e-09, -7.3632885319758815e-11, +1.1757038228899609e-12, -1.9116931765458313e-14,
        +8.8247211977470597e-02, -1.2632824871677084e-03, +1.8327135563069185e-05, -2.6961299017329128e-07, +4.0246185600026001e-09, -6.1005018388276936e-11, +9.4024359978573407e-13, -1.4735365251589550e-14,
        +8.5791861141525455e-02, -1.1930852193253270e-03, +1.6801384125989369e-05, -2.3971372206000537e-07, +3.4670799022552322e-09, -5.0866728670840542e-11, +7.5790173064224896e-13, -1.1467704365899274e-14,
        +8.3471032422208435e-02, -1.1286492297719218e-03, +1.5442416791276834e-05, -2.1389665923062143e-07, +3.0008432589756243e-09, -4.2665144356035804e-11, +6.1538680753953463e-13, -9.0036170058277808e-15,
        +8.1273839143014809e-02, -1.0693535945361543e-03, +1.4227785678185235e-05, -1.9150290893206423e-07, +2.6087212374286166e-09, -3.5983572775739499e-11, +5.0305192259941166e-13, -7.1266038881476579e-15,
        +7.9190551692747024e-02, -1.0146597928761005e-03, +1.3138614924621800e-05, -1.7199583877100006e-07, +2.2771707371987486e-09, -3.0504745421845618e-11, +4.1381018280047739e-13, -5.6834496796569005e-15,
        +7.7212446084330782e-02, -9.6409877720137437e-04, +1.2158947207883668e-05, -1.5493552306122552e-07, +1.9954360000926044e-09, -2.5984827829823392e-11, +3.4239766130801416e-13, -4.5643131795897959e-15,
        +7.5331675951374633e-02, -9.1726037574459791e-04, +1.1275224935546244e-05, -1.3995911266515999e-07, +1.7549124341542817e-09, -2.2234800246819142e-11, +2.8486473515217477e-13, -3.6895146958327518e-15,
        +7.3541163811577659e-02, -8.7378455386745912e-04, +1.0475875325005284e-05, -1.2676562391639020e-07, +1.5486695048371910e-09, -1.9106959402367977e-11, +2.3822014886799117e-13, -3.0006254016637247e-15,
        +7.1834508262308305e-02, -8.3335416676778233e-04, +9.7509762990201585e-06, -1.1510405674244892e-07, +1.3710895481498637e-09, -1.6485025748472684e-11, +2.0017879092435479e-13, -2.4543699501858687e-15,
        +7.0205904423928475e-02, -7.9568891683424526e-04, +9.0919857388399205e-06, -1.0476404371094838e-07, +1.2175919535142251e-09, -1.4276819975510020e-11, +1.6898088706726977e-13, -2.0183985427738386e-15,
        +6.8650075456038601e-02, -7.6054029005029799e-04, +8.4915207686462462e-06, -9.5568439367607789e-08, +1.0844208284142313e-09, -1.2408790869632695e-11, +1.4326100127336244e-13, -1.6683194677946308e-15,
        +6.7162213374107577e-02, -7.2768729271319025e-04, +7.9431768161630939e-06, -8.7367408759759650e-08, +9.6848029518857600e-10, -1.0821890964424908e-11, +1.2195218231844151e-13, -1.3855928696048079e-15,
        +6.5737927713689845e-02, -6.9693284593088593e-04, +7.4413784975499690e-06, -8.0033682830865548e-08, +8.6720582832910021e-10, -9.4684428028793220e-12, +1.0421517972572249e-13, -1.1560201969323073e-15,
        +6.4373200845182374e-02, -6.6810072352358103e-04, +6.9812561167096256e-06, -7.3458728290039648e-08, +7.7846307612043544e-10, -8.3097397863465342e-12, +8.9385722227591393e-14, -9.6865034188087377e-16,
        +6.3064348947871496e-02, -6.4103294102184566e-04, +6.5585428970840791e-06, -6.7549638795860712e-08, +7.0046779637388140e-10, -7.3141963052676452e-12, +7.6934932739432183e-14, -8.1498033967514177e-16,
        +6.1807987818557009e-02, -6.1558752085310903e-04, +6.1694890838117831e-06, -6.2226598567290446e-08, +6.3172212503010681e-10, -6.4559118083249710e-12, +6.6439381345910442e-14, -6.8836634432720237e-16,
        +6.0601002825527370e-02, -5.9163657261286824e-04, +5.8107898427345094e-06, -5.7420802885831171e-08, +5.7096356173804429e-10, -5.7135491136453304e-12, +5.7558267908792867e-14, -5.8358609490599391e-16,
        +5.9440522429442789e-02, -5.6906463833281324e-04, +5.4795244965878050e-06, -5.3072745262588680e-08, +5.1712391822098210e-10, -5.0694529347377374e-12, +5.0015919716691145e-14, -4.9651146475949494e-16,
        +5.8323894783699794e-02, -5.4776726149654664e-04, +5.1731051194797189e-06, -4.9130800384722743e-08, +4.6929611739803618e-10, -4.5089532384328436e-12, +4.3588280468158365e-14, -4.2386165508049453e-16,
        +5.7248667001968907e-02, -5.2764974568192803e-04, +4.8892328895187610e-06, -4.5550046832888333e-08, +4.2670721329730596e-10, -4.0198117051823824e-12, +3.8092417313811543e-14, -3.6301592032879143e-16,
        +5.6212566742861854e-02, -5.0862607448652944e-04, +4.6258608995017245e-06, -4.2291285076304129e-08, +3.8869636620263629e-10, -3.5917796370722572e-12, +3.3378324990589587e-14, -3.1187055931336622e-16,
        +5.5213485813507816e-02, -4.9061796909721919e-04, +4.3811623645279894e-06, -3.9320215217598320e-08, +3.5469678475212367e-10, -3.2162431477185307e-12, +2.9322489199338614e-14, -2.6872904704244539e-16,
        +5.4249465537121129e-02, -4.7355406371288154e-04, +4.1535033566666556e-06, -3.6606745979335536e-08, +3.2422085900450359e-10, -2.8859370731506427e-12, +2.5822805289807381e-14, -2.3221714046482350e-16,
        +5.3318683665961195e-02, -4.5736918219010068e-04, +3.9414193506210202e-06, -3.4124411949839696e-08, +2.9684787200229643e-10, -2.5947132652504404e-12, +2.2794546952422841e-14, -2.0121692072330869e-16,
        +5.2419442651656357e-02, -4.4200370188857412e-04, +3.7435949886076293e-06, -3.1849880477785496e-08, +2.7221380401109568e-10, -2.3373521301295662e-12, +2.0167152764527101e-14, -1.7481527905222660e-16,
    };
    const int s = (int)((x-aCheb4)*invSubwidth); // index of subrange
    faddeeva_nofterms = s;
    const double center = ((nSubranges-0.5)-s)*(aCheb4/nSubranges) + (s+0.5) * (bCheb4/nSubranges);
    const double t = 2 * invSubwidth * (x-center); // coord in subrange, between -1 and +1
    const double *const c = Coeffs + s * 8;
    return (((((( c[7] * t + c[6] ) * t + c[5] ) * t + c[4] ) * t + c[3] ) * t + c[2] ) * t + c[1] ) * t + c[0];
}
// clang-format on
//--- End of autogenerated code

/******************************************************************************/
/*  Library function im_w_of_z                                                */
/******************************************************************************/

double im_w_of_x(double x)
{
    // Steven G. Johnson, October 2012.
    // Revised for relative accuracy better than 4 epsilon by Joachim Wuttke, Sept 2024.

    // Uses three different methods:
    // - asymptotic expansion for large |x|,
    // - Chebyshev polynomials for medium |x|,
    // - Taylor (Maclaurin) series for small |x|.

    const double ispi = 0.56418958354775628694807945156; // 1 / sqrt(pi)
    const double ax = fabs(x); // very fast

    if (ax > bCheb4) {
        // Use asymptotic expansion up to N = 0, 3, 6, or 10

	faddeeva_algorithm = 550;
	// With N=15 or 20 we could extend the range to 7.73 or 6.72,
	// but we expect Chebyshev to be faster.
	if (ax > 125) {
	    if (ax > 6.6e7) { // 1-term expansion, important to avoid overflow
		faddeeva_nofterms = 1; // N = 0
		return ispi / x;
	    }
	    faddeeva_nofterms = 4; // N = 3
	    const double r = 1/x;
	    const double r2 = r*r;
	    return ispi * r * ((((
				     + 1.875) * r2 // coefficient (2N-1)!!/2^N
				 + 0.75) * r2
				+ 0.5) * r2
			       + 1);
	}
	const double r = 1/x;
	const double r2 = r*r;
	if (ax > 22.7) {
	    faddeeva_nofterms = 7; // N=6
	    return ispi * r * (((((((
					+ 162.421875 ) * r2
				    + 29.53125 ) * r2
				   + 6.5625 ) * r2
				  + 1.875 ) * r2
				 + 0.75 ) * r2
				+ 0.5 ) * r2
			       + 1);
	}
	faddeeva_nofterms = 11; // N=10
	return ispi * r * (((((((((((
					+ 639383.8623046875 ) * r2
				    + 67303.564453125 ) * r2
				   + 7918.06640625 ) * r2
				  + 1055.7421875 ) * r2
				 + 162.421875 ) * r2
				+ 29.53125 ) * r2
			       + 6.5625 ) * r2
			      + 1.875 ) * r2
			     + 0.75 ) * r2
			    + 0.5 ) * r2
			   + 1);
    }

    if (ax < aCheb1) {
        // Use Taylor expansion (2/sqrt(pi)) * (x - 2/3 x^3  + 4/15 x^5  - 8/105 x^7 + ...)

	faddeeva_algorithm = 500;
	const double x2 = x*x;
	if (ax < 0.016) {
	    faddeeva_nofterms = 4;
	    return (((( - 0.085971746064420005629 ) * x2// x^7
		      + 0.30090111122547001970 ) * x2 // x^5
		     - 0.75225277806367504925 ) * x2 // x^3
		    + 1.1283791670955125739 ) * x;
	}
	if (ax < .29) {
	    faddeeva_nofterms = 9;
	    return ((((((((( + 8.38275934019361123956e-6 ) * x2 // x^17
			   - 7.1253454391645686483238e-5 ) * x2 // x^15
			  + 0.00053440090793734269229 ) * x2 // x^13
			 - 0.0034736059015927275001 ) * x2 // x^11
			+ 0.019104832458760001251 ) * x2 // x^9
		       - 0.085971746064420005629 ) * x2 // x^7
		      + 0.30090111122547001970 ) * x2 // x^5
		     - 0.75225277806367504925 ) * x2 // x^3
		    + 1.1283791670955125739 ) * x;
	}
	faddeeva_nofterms = 17;
        return ((((((((((((((((( + 1.16774718055184835728293189e-14 ) * x2 // x^33
			       - 1.92678284791054972871829131e-13 ) * x2 //x^31
			      + 2.98651341426135223029374655e-12 ) * x2 // x^29
			     - 4.33044445067896090883119155e-11 ) * x2 // x^27
			    + 5.8461000084165966602290712e-10 ) * x2 // x^25
			   - 7.30762501052074563638866034e-9 ) * x2 // x^23
			  + 8.40376876209885782941868884e-8 ) * x2 // x^21
			 - 8.82395720020380130481012927e-7 ) * x2 // x^19
			+ 8.38275934019361123956e-6 ) * x2 // x^17
		       - 7.1253454391645686483238e-5 ) * x2 // x^15
		      + 0.00053440090793734269229 ) * x2 // x^13
		     - 0.0034736059015927275001 ) * x2 // x^11
		    + 0.019104832458760001251 ) * x2 // x^9
		   - 0.085971746064420005629 ) * x2 // x^7
		  + 0.30090111122547001970 ) * x2 // x^5
		 - 0.75225277806367504925 ) * x2 // x^3
		+ 1.1283791670955125739 ) * x;
    }

    // Remaining intermediate range:
    // Use Chebyshev interpolants.

    if (ax < bCheb2) {
	if (ax < bCheb1) {
	    faddeeva_algorithm = 510;
	    return copysign(chebInterpolant1(ax), x);
	}
	faddeeva_algorithm = 520;
	return copysign(chebInterpolant2(ax), x);
    }
    if (ax < bCheb4) {
	if (ax < bCheb3) {
	    faddeeva_algorithm = 530;
	    return copysign(chebInterpolant3(ax), x);
	}
	faddeeva_algorithm = 540;
	return copysign(chebInterpolant4(ax), x);
    }

    return NaN;
} // im_w_of_z
