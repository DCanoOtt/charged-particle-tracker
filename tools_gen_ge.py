import math, json, xraylib as x
NA=6.02214076e23; A=72.630; Z=32; RHO=5.323; re2=7.9407877e-26; mec2=0.51099895
nist={1.0e-2:37.42,1.5e-2:91.52,2.0e-2:42.22,3.0e-2:13.85,4.0e-2:6.207,5.0e-2:3.335,6.0e-2:2.023,8.0e-2:0.9501,0.1:0.5550,0.15:0.2491,0.2:0.1661,0.3:0.1131,0.4:0.09327,0.5:0.08212,0.6:0.07452,0.8:0.06426,
 1.0:0.05727,1.25:0.05101,1.5:0.04657,2.0:0.04086,3.0:0.03524,4.0:0.03275,5.0:0.03158,6.0:0.03107,8.0:0.03103,10.0:0.03156}   # NIST X-ray mass attenuation coefficients, Ge, mu/rho in cm2/g
def kn(E):
    k=E/mec2; l=math.log(1+2*k)
    return 2*math.pi*re2*((1+k)/k**2*(2*(1+k)/(1+2*k)-l/k)+l/(2*k)-(1+3*k)/(1+2*k)**2)
grid=[10,11.1030,11.1032,15,20,30,40,50,60,80,100,150,200,300,400,500,600,800,1000,1250,1500,2000,3000,4000,5000,6000,8000,10000]
rows=[]; chk=[]
ph800=x.CS_Photo(32,800); ra800=x.CS_Rayl(32,800)
for k in grid:
    E=k/1000
    if k<=800:
        ph,co,ra,tot=x.CS_Photo(32,k),x.CS_Compt(32,k),x.CS_Rayl(32,k),x.CS_Total(32,k); pp=0.0
        if k in [round(e*1000) for e in nist]: chk.append((k,tot/nist[E]-1))
    else:
        co=Z*kn(E)*NA/A; ph=ph800*(800/k)**1.8; ra=ra800*(800/k)**2; tot=nist[E]
        pp=max(0.0,tot-co-ph-ra) if k>1022 else 0.0
        if k<=1022: tot=co+ph+ra   # below threshold: sum of the components (NIST total differs by <0.5 %)
    rows.append(dict(E_keV=k,mu_pe=ph,mu_incoh=co,mu_coh=ra,mu_pair=pp,mu_total=tot if k>800 else tot))
print('xraylib total vs NIST total (E<=800 keV): max |dev| = %.3f %%'%(100*max(abs(c[1]) for c in chk)))
for r in rows: print('%9.4f  pe %.4e  inc %.4e  coh %.4e  pair %.4e  tot %.4e'%(r['E_keV'],r['mu_pe'],r['mu_incoh'],r['mu_coh'],r['mu_pair'],r['mu_total']))
json.dump(dict(element='Ge',Z=Z,A=A,density_g_cm3=RHO,units='cm2/g',rows=rows,source='NIST XrayMassCoef (total, E>800 keV); xraylib 4.x (EPDL97/Elam) photoelectric, incoherent, coherent for E<=800 keV; above: free Klein-Nishina, photoelectric extrapolated E^-1.8, pair = total - others'),open('ge_cross_sections.json','w'),indent=1)
with open('ge_cross_sections.csv','w') as f:
    f.write('# Germanium photon mass attenuation coefficients (cm2/g), density 5.323 g/cm3\n# E<=800 keV: xraylib components; above: NIST total, Klein-Nishina incoherent, photoelectric extrapolated, pair = remainder\nE_keV,photoelectric,incoherent,coherent,pair,total\n')
    for r in rows: f.write('%g,%.5e,%.5e,%.5e,%.5e,%.5e\n'%(r['E_keV'],r['mu_pe'],r['mu_incoh'],r['mu_coh'],r['mu_pair'],r['mu_total']))
