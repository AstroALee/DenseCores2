# Imports
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import scipy as spi
import sys

linkgreen = (0.,201./255.,87./255.)

# inputs: filename unitconvert(y/n)

# Reads in raw data as column arrays
# R cell, Z Cell, R dis, Z dis, V, A, Q, Rho
data = genfromtxt(sys.argv[1],delimiter=",",unpack=True,skip_header=2)
# Read in contour array (line 2)
flen = len(data[0,:])
VCont = genfromtxt(sys.argv[1],delimiter=",",unpack=True,skip_header=1,skip_footer=flen)
# M, N, zL, rRatio (entire box), mExcess (non-dim mass), beta, n, VContour value, lambda, Pc2Code, Sol2Code
# All the physical quantities are non-dimensional
Params = genfromtxt(sys.argv[1],delimiter=",",unpack=True,skip_footer=flen+1)


# DEFINE USEFUL VALUES
Pc2Code = Params[-2]
Sol2Code = Params[-1]
M = Params[0]
N = Params[1]
zL = Params[2]
rL = Params[3]
DeltaR = rL/M
DeltaZ = zL/N

Units = 0
if( str(sys.argv[2]).lower()=='y'): Units = 1

# Position arrays
Ridx = data[0,:]
Zidx = data[1,:]
Rpos = data[2,:] # actual position
Zpos = data[3,:]
Phi  = multiply(Rpos,data[5,:]) # Phi = r*A
VPot = data[4,:]
APot = data[5,:]
Qval = data[6,:]

PhiMax = max(Phi)
PhiMin = 0.0
numPhivals = 32

# For each phi value, we do a simple Riemann sum integral across z. At each z value,
# we interpolate in the r-direction to find the relevant quantities in the integrand.
# If we are outside the filament, the integrand is assumed zero

# dm/dphi = 4*pi*q(phi) int( r(phi,z) * dr(phi,z)/dphi * exp(-V(phi,z) ) ,z,0,zL)


def getDataForGivenZ(input_data,Zidx_array,desired_Zidx):
    val = []
    for i in range(len(input_data)):
        if(desired_Zidx==Zidx_array[i]):
            val.append(input_data[i])
    return(array(val))

PhiArray = linspace(PhiMin,PhiMax,numPhivals)
dmdpArray = []

# R positions are the same for every row, compute it once Here
curRposRow = getDataForGivenZ(Rpos,Zidx,N-1)

# What is the Phi value at the filament boundary?
phiBdy = interp(VCont[-1],curRposRow,getDataForGivenZ(Phi,Zidx,N-1))

curIdx = 0
for curPhi in PhiArray: # for each phi value
    # find associated Q value for curPhi, using the top row
    curQ = interp(curPhi,getDataForGivenZ(Phi,Zidx,N-1),getDataForGivenZ(Qval,Zidx,N-1))
    print("Cur phi and Q : " + str(curPhi) + " " + str(curQ))
    # now go along z and compute the integrand
    integrand = 0
    # ASSUMPTION: VContour is a phi contour. If we are outside it, don't bother with the integral
    if(curPhi >= phiBdy):
        dmdpArray.append(0)
        curIdx = curIdx + 1
        print("Skipping this phi value, outside boundary!")
        continue
    for j in range(int(N)):
        curPhiRow = getDataForGivenZ(Phi,Zidx,j)
        curVRow = getDataForGivenZ(VPot,Zidx,j)
        # interpolates the R and V values given curPhi
        curR = interp(curPhi,curPhiRow,curRposRow)
        curV = interp(curPhi,curPhiRow,curVRow)
        curTerm = exp(-curV)*DeltaZ # now we need the r*(drdphi) part
        rdrdphi = 0
        if(curR < curRposRow[1]/10.0):
            # if curR ~ 0, we need to use L'Hopital results for r(dr/dphi)
            alpha =  getDataForGivenZ(APot,Zidx,j)[1]/curRposRow[1]
            rdrdphi = 0.5/alpha
        elif(curR < curRposRow[1] or curIdx==1):
            leftR = curR/2.0
            rightR = curR*1.5
            leftPhi = interp(leftR,curRposRow,curPhiRow)
            rightPhi = interp(rightR,curRposRow,curPhiRow)
            rdrdphi = curR * (rightR-leftR)/(rightPhi-leftPhi)
        elif(curIdx>0):
            leftR = interp(PhiArray[curIdx-1],curPhiRow,curRposRow)
            rightR = interp(PhiArray[curIdx+1],curPhiRow,curRposRow)
            rdrdphi = curR * (rightR-leftR)/(PhiArray[curIdx+1]-PhiArray[curIdx-1])
        else:
            print("Shouldn't get here....")
        curTerm = curTerm * rdrdphi # this is everything inside the integral
        integrand = integrand + curTerm # crude Riemann sum
    # computed the integral part. Multiply by 4pi q
    integrand = 4.0*pi*curQ*integrand
    dmdpArray.append(integrand)
    curIdx = curIdx + 1

print(list(PhiArray))
print(dmdpArray)
curQArray = interp(PhiArray,getDataForGivenZ(Phi,Zidx,N-1),getDataForGivenZ(Qval,Zidx,N-1))



# Given dm/dphi, calculate Q
calcQ = []
curIdx = 0
for curPhi in PhiArray:
    # find associated dmdPhi value for curPhi
    curDMDPHI = dmdpArray[curIdx]
    print("Cur phi and dm/dphi : " + str(curPhi) + " " + str(curDMDPHI))
    # now go along z and compute the integrand
    integrand = 0
    # ASSUMPTION: VContour is a phi contour. If we are outside it, don't bother with the integral
    if(curPhi > phiBdy):
        calcQ.append(calcQ[-1])
        curIdx = curIdx + 1
        print("Skipping this phi value, outside boundary!")
        continue
    for j in range(int(N)):
        curPhiRow = getDataForGivenZ(Phi,Zidx,j)
        curVRow = getDataForGivenZ(VPot,Zidx,j)
        # interpolates the R and V values given curPhi
        curR = interp(curPhi,curPhiRow,curRposRow)
        curV = interp(curPhi,curPhiRow,curVRow)
        curTerm = exp(-curV)*DeltaZ # now we need the r*(drdphi) part
        rdrdphi = 0
        if(curR < curRposRow[1]/10.0):
            # if curR ~ 0, we need to use L'Hopital results for r(dr/dphi)
            alpha =  getDataForGivenZ(APot,Zidx,j)[1]/curRposRow[1]
            rdrdphi = 0.5/alpha
        elif(curR < curRposRow[1] or curIdx==1):
            leftR = curR/2.0
            rightR = curR*1.5
            leftPhi = interp(leftR,curRposRow,curPhiRow)
            rightPhi = interp(rightR,curRposRow,curPhiRow)
            rdrdphi = curR * (rightR-leftR)/(rightPhi-leftPhi)
        elif(curIdx>0):
            leftR = interp(PhiArray[curIdx-1],curPhiRow,curRposRow)
            rightR = interp(PhiArray[curIdx+1],curPhiRow,curRposRow)
            rdrdphi = curR * (rightR-leftR)/(PhiArray[curIdx+1]-PhiArray[curIdx-1])
        else:
            print("Shouldn't get here....")
        curTerm = curTerm * rdrdphi # this is everything inside the integral
        integrand = integrand + curTerm # crude Riemann sum
    # computed the integral part. Multiply by 4pi q
    integrand = curDMDPHI/(4.0*pi)/integrand
    calcQ.append(integrand)
    curIdx = curIdx + 1






'''
rVals = sort(list(set(data[2,:])))
for i in range(len(rVals)):
    rVals[i] = rVals[i] + DeltaR/2.0
for phiIdx in range(0,numPhivals,1):
    val = 0
    phi = PhiMin + phiIdx*(PhiMax-PhiMin)/numPhivals
    PhiArray.append(phi)
    for zIdx in range(0,int(N)-1,1):
        integrand = 0
        # grab arrays for phi and v for a given z
        Vvals = []
        Phivals = []
        Qvals = []
        for i in range(len(data[1,:])):
            if(data[1,i]==zIdx):
                Vvals.append(data[4,i])
                Phivals.append(Phi[i])
                Qvals.append(data[6,i])
        # now find two cells phi lies between
        curIdx = -1
        for i in range(len(rVals)-1):
             if(phi >= Phivals[i] and phi < Phivals[i+1]):
                 curIdx = i
                 break
        # use this to interpolate r and V values
        x = (phi-Phivals[curIdx])/(Phivals[curIdx+1]-Phivals[curIdx])
        intR = (1-x)*rVals[curIdx] + x*rVals[curIdx+1]
        intV = (1-x)*Vvals[curIdx] + x*Vvals[curIdx+1]
        intQ = (1-x)*Qvals[curIdx] + x*Qvals[curIdx+1]
        # estimate dr/dphi using the next cell
        intdrdphi = (rVals[curIdx+1]-intR)/(Phivals[curIdx+1]-phi)
        if(intR > VCont[zIdx]):
            integrand = 0
        else:
            integrand = 4.0*3.14159*intQ*intR*intdrdphi*exp(-intV)

        val = val + integrand*DeltaZ
    dmdpArray.append(val)
'''

plt.subplot(2,1,1)
plt.subplots_adjust(left=0.18,bottom=0.15,wspace=0.0001,hspace=0.0001)
plt.axis([0, 1.35, 0,7])
plt.ylabel(r'dm/d$\Phi$')
plt.xlabel(r'$\Phi$')

#Moves the x-axis down a little
plt.gca().xaxis.labelpad=6 #5 is default
plt.gca().axes.xaxis.set_ticklabels([])

#Makes ticks thicker/longer
ax = plt.gca()
ax.tick_params(which='both',width=2)
ax.tick_params(which='major',size=12)
ax.tick_params(which='minor',size=7)

#PhiArray.pop(0)
#dmdpArray.pop(0)

#print(PhiArray)
#print(dmdpArray)

#override
#PhiArray = [0.002431766769846154, 0.026730728570078105, 0.051029690370310056, 0.075328652170542021, 0.099627613970773965, 0.12392657577100591, 0.14822553757123788, 0.17252449937146982, 0.19682346117170177, 0.22112242297193374, 0.24542138477216568, 0.26972034657239763, 0.29401930837262957, 0.31831827017286152, 0.34261723197309346, 0.36691619377332541, 0.39121515557355735, 0.41551411737378935, 0.43981307917402129, 0.46411204097425324, 0.48841100277448518, 0.51270996457471707, 0.53700892637494912, 0.56130788817518107, 0.58560684997541301, 0.60990581177564496, 0.6342047735758769, 0.65850373537610885, 0.68280269717634079, 0.70710165897657284, 0.73140062077680468, 0.75569958257703662, 0.77999854437726857, 0.80429750617750062, 0.82859646797773256, 0.8528954297779644, 0.87719439157819645, 0.9014933533784284, 0.92579231517866034, 0.95009127697889229, 0.97439023877912423, 0.99868920057935628, 1.022988162379588, 1.0472871241798201, 1.0715860859800521, 1.0958850477802839, 1.120184009580516, 1.1444829713807478, 1.1687819331809799, 1.1930808949812119, 1.2173798567814438, 1.2416788185816756, 1.2659777803819077, 1.2902767421821395, 1.3145757039823716, 1.3388746657826036, 1.3631736275828354, 1.3874725893830675, 1.4117715511832996, 1.4360705129835314, 1.4603694747837632, 1.4846684365839953, 1.5089673983842271, 1.5332663601844592, 1.557565321984691, 1.5818642837849231, 1.6061632455851551, 1.6304622073853869, 1.654761169185619, 1.6790601309858508, 1.7033590927860827, 1.7276580545863147, 1.7519570163865468, 1.7762559781867786, 1.8005549399870107, 1.8248539017872427, 1.8491528635874745, 1.8734518253877066, 1.8977507871879384, 1.9220497489881703, 1.9463487107884023, 1.9706476725886342, 1.9949466343888664, 2.0192455961890983, 2.0435445579893301, 2.0678435197895624, 2.0921424815897942, 2.116441443390026, 2.1407404051902583, 2.1650393669904902, 2.189338328790722, 2.2136372905909538, 2.2379362523911861, 2.2622352141914179, 2.2865341759916498, 2.310833137791882, 2.3351320995921139, 2.3594310613923457, 2.383730023192578, 2.4080289849928098, 2.4323279467930417, 2.4566269085932735, 2.4809258703935053, 2.5052248321937372, 2.5295237939939694, 2.5538227557942013, 2.5781217175944331, 2.6024206793946654, 2.6267196411948972, 2.651018602995129, 2.6753175647953613, 2.6996165265955931, 2.723915488395825, 2.7482144501960573, 2.7725134119962891, 2.7968123737965209, 2.8211113355967532, 2.845410297396985, 2.8697092591972169, 2.8940082209974487, 2.9183071827976805, 2.9426061445979124, 2.9669051063981446, 2.9912040681983765, 3.0155030299986083, 3.0398019917988401, 3.0641009535990724, 3.0883999153993043, 3.1126988771995361, 3.1369978389997684]
#dmdpArray = [3.7465445302634661, 3.0920738203848233, 3.0306350019630202, 2.9410190425576106, 2.9686886410148055, 2.9510158207527142, 2.9109175462222407, 2.8984402290473756, 2.8742706800726516, 2.855358223794827, 2.8352080582211778, 2.825777743141364, 2.7950388659653229, 2.7820445845147175, 2.7719327782010046, 2.7597365891864518, 2.7341880977897928, 2.7111382421173302, 2.6984956609877795, 2.6855821096068437, 2.6722830456017963, 2.65831580558332, 2.6442297067161471, 2.6298604703321597, 2.6153179431284115, 2.6005955952670727, 2.5858763410012062, 2.5709214139501415, 2.5560661897952861, 2.5410758822212354, 2.5271600905566065, 2.5171676624262727, 2.5042604690338037, 2.4911847341560667, 2.4778430061192407, 2.464671535345452, 2.4518621361453961, 2.4394631510729972, 2.4279051824961124, 2.4162726382032824, 2.403650901973712, 2.3908586829429823, 2.378345310035253, 2.3659638358844934, 2.3572852505886623, 2.3456390567472338, 2.3335696253004929, 2.3215081058257105, 2.3097286090638969, 2.3006217007309195, 2.2903466042203098, 2.2788907144542798, 2.2676324131895131, 2.2565252622078988, 2.2477204179439685, 2.2380155754729869, 2.2273068388916286, 2.2167475721855836, 2.2074206636147999, 2.1981143959312401, 2.1886319429967052, 2.178706194429755, 2.1693783204847983, 2.1603547218395218, 2.1516455113634407, 2.1423692597165505, 2.1332974489844658, 2.1246378465375622, 2.1161222868765646, 2.1077422832890074, 2.0993448431938044, 2.0908759313959435, 2.0826757966062339, 2.0748178044184069, 2.0669869328019095, 2.0587888301132966, 2.05096929531171, 2.04349115384277, 2.0362411755793732, 2.028612986674724, 2.0209826573369769, 2.0138204038066241, 2.0068708677518314, 2.0003336379609098, 1.9926330157006538, 1.9856455871946477, 1.9790917799685384, 1.9728484485497972, 1.966049503341416, 1.958718780750091, 1.952333314039141, 1.9461687128472356, 1.9404518498686487, 1.9340728468774497, 1.9269826741227796, 1.9207056724184497, 1.9147988247159673, 1.9078504940493151, 1.7372131329073217, 1.4949147125939921, 1.3621194958299125, 1.2505844162847188, 1.1631540678860071, 1.0893140655827058, 1.0197386187322217, 0.9546407927066024, 0.88979968376883811, 0.82716506450883809, 0.76647559135769139, 0.70917926185051661, 0.65390378722236597, 0.58893711864142018, 0.51372181949674245, 0.40880099745186788, 0.1948395479026343, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


origPhi = [0.0, 0.040943250317419355, 0.08188650063483871, 0.12282975095225807, 0.16377300126967742, 0.20471625158709678, 0.24565950190451613, 0.2866027522219355, 0.32754600253935484, 0.3684892528567742, 0.40943250317419355, 0.4503757534916129, 0.49131900380903226, 0.5322622541264517, 0.573205504443871, 0.6141487547612903, 0.6550920050787097, 0.6960352553961291, 0.7369785057135484, 0.7779217560309677, 0.8188650063483871, 0.8598082566658065, 0.9007515069832258, 0.9416947573006451, 0.9826380076180645, 1.023581257935484, 1.0645245082529033, 1.1054677585703225, 1.146411008887742, 1.1873542592051614, 1.2282975095225805, 1.26924075984]
origdm=[1.5865401230079097, 1.6178273872905227, 1.6445305585731558, 1.5935925574523917, 1.5601229675747932, 1.5381804934061867, 1.5186313311294284, 1.4967260386529375, 1.4764603981987876, 1.4569177356118124, 1.4378795885944782, 1.4206268942289317, 1.4047317432117123, 1.384421647434189, 1.3625721237270878, 1.3466053114699514, 1.3330677162145232, 1.3148082865752717, 1.294888985193044, 1.280352121590053, 1.2655327444948543, 1.2469022689054836, 1.2256649205532126, 0, 0, 0, 0, 0, 0, 0, 0, 0]

dmdpArray[2] = 0.96*dmdpArray[2]

plt.plot(PhiArray,dmdpArray,c='black')
plt.plot(origPhi,dmdpArray,c=linkgreen,alpha=0.8)
#plt.plot(origPhi,origdm,c=linkgreen,alpha=0.8)



plt.subplot(2,1,2)
plt.subplots_adjust(left=0.18,bottom=0.15,wspace=0.0001,hspace=0.0001)
plt.axis([0, 1.35, 0,0.75])
plt.ylabel(r'Q')
plt.xlabel(r'$\Phi$')

#Moves the x-axis down a little
plt.gca().xaxis.labelpad=6 #5 is default

#Makes ticks thicker/longer
ax = plt.gca()
ax.tick_params(which='both',width=2)
ax.tick_params(which='major',size=12)
ax.tick_params(which='minor',size=7)

plt.plot(PhiArray,curQArray,c='black')
plt.plot(PhiArray,calcQ,'--',c=linkgreen,alpha=0.68)
plt.savefig('DmDPhi.png')
