# Imports
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import scipy as spi
import sys



plt.subplot(1,2,1)
plt.subplots_adjust(left=0.15,bottom=0.15,wspace=0.0001,hspace=0.0001)
plt.axis([0, 3, 0,4])
plt.ylabel(r'$dm/d\Phi$',fontsize=18)
plt.xlabel(r'$\Phi$',fontsize=18)

#Moves the x-axis down a little
plt.gca().xaxis.labelpad=6 #5 is default

#Makes ticks thicker/longer
ax = plt.gca()
ax.tick_params(which='both',width=2)
ax.tick_params(which='major',size=12)
ax.tick_params(which='minor',size=7)

#PhiArray.pop(0)


PhiArray = [0.002431766769846154, 0.026730728570078105, 0.051029690370310056, 0.075328652170542021, 0.099627613970773965, 0.12392657577100591, 0.14822553757123788, 0.17252449937146982, 0.19682346117170177, 0.22112242297193374, 0.24542138477216568, 0.26972034657239763, 0.29401930837262957, 0.31831827017286152, 0.34261723197309346, 0.36691619377332541, 0.39121515557355735, 0.41551411737378935, 0.43981307917402129, 0.46411204097425324, 0.48841100277448518, 0.51270996457471707, 0.53700892637494912, 0.56130788817518107, 0.58560684997541301, 0.60990581177564496, 0.6342047735758769, 0.65850373537610885, 0.68280269717634079, 0.70710165897657284, 0.73140062077680468, 0.75569958257703662, 0.77999854437726857, 0.80429750617750062, 0.82859646797773256, 0.8528954297779644, 0.87719439157819645, 0.9014933533784284, 0.92579231517866034, 0.95009127697889229, 0.97439023877912423, 0.99868920057935628, 1.022988162379588, 1.0472871241798201, 1.0715860859800521, 1.0958850477802839, 1.120184009580516, 1.1444829713807478, 1.1687819331809799, 1.1930808949812119, 1.2173798567814438, 1.2416788185816756, 1.2659777803819077, 1.2902767421821395, 1.3145757039823716, 1.3388746657826036, 1.3631736275828354, 1.3874725893830675, 1.4117715511832996, 1.4360705129835314, 1.4603694747837632, 1.4846684365839953, 1.5089673983842271, 1.5332663601844592, 1.557565321984691, 1.5818642837849231, 1.6061632455851551, 1.6304622073853869, 1.654761169185619, 1.6790601309858508, 1.7033590927860827, 1.7276580545863147, 1.7519570163865468, 1.7762559781867786, 1.8005549399870107, 1.8248539017872427, 1.8491528635874745, 1.8734518253877066, 1.8977507871879384, 1.9220497489881703, 1.9463487107884023, 1.9706476725886342, 1.9949466343888664, 2.0192455961890983, 2.0435445579893301, 2.0678435197895624, 2.0921424815897942, 2.116441443390026, 2.1407404051902583, 2.1650393669904902, 2.189338328790722, 2.2136372905909538, 2.2379362523911861, 2.2622352141914179, 2.2865341759916498, 2.310833137791882, 2.3351320995921139, 2.3594310613923457, 2.383730023192578, 2.4080289849928098, 2.4323279467930417, 2.4566269085932735, 2.4809258703935053, 2.5052248321937372, 2.5295237939939694, 2.5538227557942013, 2.5781217175944331, 2.6024206793946654, 2.6267196411948972, 2.651018602995129, 2.6753175647953613, 2.6996165265955931, 2.723915488395825, 2.7482144501960573, 2.7725134119962891, 2.7968123737965209, 2.8211113355967532, 2.845410297396985, 2.8697092591972169, 2.8940082209974487, 2.9183071827976805, 2.9426061445979124, 2.9669051063981446, 2.9912040681983765, 3.0155030299986083, 3.0398019917988401, 3.0641009535990724, 3.0883999153993043, 3.1126988771995361, 3.1369978389997684]
dmdpArray = [3.7465445302634661, 3.3920738203848233, 3.2306350019630202, 3.1410190425576106, 3.0686886410148055, 3.009510158207527142, 2.942109175462222407, 2.8984402290473756, 2.8742706800726516, 2.855358223794827, 2.8352080582211778, 2.825777743141364, 2.7950388659653229, 2.7820445845147175, 2.7719327782010046, 2.7597365891864518, 2.7341880977897928, 2.7111382421173302, 2.6984956609877795, 2.6855821096068437, 2.6722830456017963, 2.65831580558332, 2.6442297067161471, 2.6298604703321597, 2.6153179431284115, 2.6005955952670727, 2.5858763410012062, 2.5709214139501415, 2.5560661897952861, 2.5410758822212354, 2.5271600905566065, 2.5171676624262727, 2.5042604690338037, 2.4911847341560667, 2.4778430061192407, 2.464671535345452, 2.4518621361453961, 2.4394631510729972, 2.4279051824961124, 2.4162726382032824, 2.403650901973712, 2.3908586829429823, 2.378345310035253, 2.3659638358844934, 2.3572852505886623, 2.3456390567472338, 2.3335696253004929, 2.3215081058257105, 2.3097286090638969, 2.3006217007309195, 2.2903466042203098, 2.2788907144542798, 2.2676324131895131, 2.2565252622078988, 2.2477204179439685, 2.2380155754729869, 2.2273068388916286, 2.2167475721855836, 2.2074206636147999, 2.1981143959312401, 2.1886319429967052, 2.178706194429755, 2.1693783204847983, 2.1603547218395218, 2.1516455113634407, 2.1423692597165505, 2.1332974489844658, 2.1246378465375622, 2.1161222868765646, 2.1077422832890074, 2.0993448431938044, 2.0908759313959435, 2.0826757966062339, 2.0748178044184069, 2.0669869328019095, 2.0587888301132966, 2.05096929531171, 2.04349115384277, 2.0362411755793732, 2.028612986674724, 2.0209826573369769, 2.0138204038066241, 2.0068708677518314, 2.0003336379609098, 1.9926330157006538, 1.9856455871946477, 1.9790917799685384, 1.9728484485497972, 1.966049503341416, 1.958718780750091, 1.952333314039141, 1.9461687128472356, 1.9404518498686487, 1.9340728468774497, 1.9269826741227796, 1.9207056724184497, 1.9147988247159673, 1.9078504940493151, 1.7372131329073217, 1.4949147125939921, 1.3621194958299125, 1.2505844162847188, 1.1631540678860071, 1.0893140655827058, 1.0197386187322217, 0.9546407927066024, 0.88979968376883811, 0.82716506450883809, 0.76647559135769139, 0.70917926185051661, 0.65390378722236597, 0.58893711864142018, 0.51372181949674245, 0.40880099745186788, 0.1948395479026343, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

plt.plot(PhiArray,dmdpArray,c='black')



plt.subplot(1,2,2)
plt.subplots_adjust(left=0.15,bottom=0.15,wspace=0.5001,hspace=0.0001)
plt.axis([0, 3, 0,4])
plt.ylabel(r'$dm/d\Phi$',fontsize=18)
plt.xlabel(r'$\Phi$',fontsize=18)

#Moves the x-axis down a little
plt.gca().xaxis.labelpad=6 #5 is default

#Makes ticks thicker/longer
ax = plt.gca()
ax.tick_params(which='both',width=2)
ax.tick_params(which='major',size=12)
ax.tick_params(which='minor',size=7)

#PhiArray.pop(0)

plt.plot(PhiArray,dmdpArray,c='blue')



plt.savefig('TwoPanel.pdf')