import matplotlib.pyplot as plt
import numpy as np

plt.plot([0.0007540933440450814, 0.0007540920884863134, 0.000754090946215813, 0.0007540899172342463, 0.0007540890015417245, 0.0007540881991384696, 0.0007540875100242594, 0.0007540869341998713, 0.0007540864716654161, 0.0007540861224206719,
          0.0007540858864659716, 0.0007540857638022036, 0.0007540857544283686, 0.0007540858583456878, 0.0007540860755539391, 0.0007540864060533448, 0.0007540868498443487, 0.000754087406927173, 0.0007540880773013736, 0.0007540888609676166, 0.0007540897579263461])
# plt.plot([1.545929322332995e-07, 1.5459282010077402e-07, 1.5459270852336005e-07, 1.5459259672390147e-07, 1.5459248492444289e-07, 1.54592373013962e-07, 1.5459226110348112e-07, 1.5459214908197794e-07, 1.5459203739354166e-07, 1.5459192559408308e-07, 1.545918140166691e-07, 1.5459170210618822e-07, 1.5459158997366274e-07, 1.5459147828522646e-07, 1.5459136637474558e-07, 1.545912544642647e-07, 1.5459114310889532e-07, 1.5459103064330293e-07, 1.5459091939895586e-07, 1.545908078215419e-07, 1.54590695911061e-07, 1.5459058433364703e-07, 1.5459047197907694e-07, 1.5459036051268527e-07, 1.545902489352713e-07,
        #   1.545901369137681e-07, 1.5459002522533183e-07, 1.5458991298178404e-07, 1.5458980195948158e-07, 1.545896902710453e-07, 1.5458957847158672e-07, 1.5458946689417274e-07, 1.5458935520573647e-07, 1.545892434062779e-07, 1.54589131495797e-07, 1.5458902002940533e-07, 1.5458890878505827e-07, 1.5458879665253278e-07, 1.545886845200073e-07, 1.5458857316463792e-07, 1.5458846158722395e-07, 1.5458834989878767e-07, 1.54588238432396e-07, 1.5458812652191511e-07, 1.5458801505552344e-07, 1.5458790358913177e-07, 1.545877920117178e-07, 1.5458768021225922e-07, 1.5458756830177833e-07, 1.5458745705743127e-07, 1.545873454800173e-07])
# plt.plot([1.5559911092921652e-07, 1.554981103879527e-07, 1.5540245035250422e-07, 1.5531213248820563e-07, 1.5522715635096773e-07, 1.551475214967013e-07, 1.5507322803642865e-07, 1.5500427608117207e-07, 1.5494066607502077e-07, 1.5488239712979635e-07, 1.548294700226549e-07, 1.5478188397644033e-07, 1.5473963999035334e-07, 1.5470273684314861e-07, 1.5467117597811608e-07, 1.5464495617401042e-07, 1.5462407809696543e-07, 1.5460854141391422e-07, 1.5459834590281218e-07, 1.5459349189672622e-07, 1.5459398028383475e-07, 1.5459980917675864e-07, 1.5461098046287702e-07, 1.5462749314298918e-07, 1.546493467730059e-07,
#           1.54676542019061e-07, 1.547090792142214e-07, 1.5474695780337555e-07, 1.5479017767550118e-07, 1.548387391636652e-07, 1.5489264204582298e-07, 1.5495188643299684e-07, 1.550164725472314e-07, 1.55086400166482e-07, 1.551616692907487e-07, 1.5524227947594227e-07, 1.5532823205433033e-07, 1.5541952524955605e-07, 1.5551616039388705e-07, 1.5561813726527873e-07, 1.557254554196419e-07, 1.5583811507902112e-07, 1.5595611602137183e-07, 1.5607945913487242e-07, 1.5620814286521068e-07, 1.5634216854465421e-07, 1.5648153617320304e-07, 1.5662624486267873e-07, 1.567762950571705e-07, 1.5693168675667835e-07, 1.5709242040529148e-07])
# plt.plot([0.009708456107103913, 0.008742836760692985, 0.007827363296351919, 0.006962129508943926, 0.0061472240465533234, 0.005382730401316982, 0.004668726900793674, 0.004005286699863553, 0.0033924777731667577, 0.0028303629080732584, 0.002318999698191826, 0.001858440537414796, 0.00144873261450007, 0.0010899179081933497, 0.0007820331828877203, 0.0005251099848189167, 0.00031917463880914987, 0.00016424824553884498, 6.034667937027116e-05, 7.480586702413028e-06, 5.655384871738889e-06, 5.4871261591760145e-05, 0.0001551231749326032, 0.00030640085384348037, 0.000508688799212953,
#           0.0007619662854739806, 0.001066207362742544, 0.0014213808595050503, 0.0018274503858373148, 0.0022843743371713243, 0.002792105898591246, 0.0033505930496776637, 0.003959778569885941, 0.004619600044459049, 0.005329989870887064, 0.006090875265887696, 0.006902178272939485, 0.00776381577033769, 0.008675699479788523, 0.009637735975538186, 0.01064982669402914, 0.011711867944093957, 0.012823750917672738, 0.013985361701064347, 0.015196581286702093, 0.016457285585460335, 0.017767345439482107, 0.01912662663553122, 0.020534989918868263, 0.02199229100765021, 0.023498380607843816])
plt.xlabel('t_mesh')
plt.ylabel('loss')
plt.show()
