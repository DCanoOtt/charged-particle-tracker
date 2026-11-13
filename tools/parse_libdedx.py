import re,numpy as np
E='/tmp/claude-0/-home-claude/f4455f09-14fc-5dfb-bbf4-c0eee58c6786/scratchpad/ext/libdedx/src/data/embedded/'
def arr(txt,name):
    m=re.search(r'static const (float|int) '+name+r'(\[[^\]]*\])+ = \{(.*?)\};',txt,re.S)
    return np.array([float(x) for x in re.findall(r'-?\d+\.?\d*(?:e[+-]?\d+)?',m.group(3).replace('f',''))])
def load(n):
    t=re.sub(r'/\*.*?\*/','',open(E+'dedx_%s.h'%n).read(),flags=re.S)
    en=arr(t,'dedx_%s_energy'%n); ids=arr(t,'dedx_%s_ion_ids'%n).astype(int); tg=arr(t,'dedx_%s_target_ids'%n).astype(int); s=arr(t,'dedx_%s_stp'%n).reshape(len(ids),len(tg),len(en))
    return en,ids,tg,s
if __name__=='__main__':
    for n in ['astar','icru_astar','pstar','icru73','mstar']:
        en,ids,tg,s=load(n); print(n,en[[0,-1]],ids,tg.tolist(),s.shape)
    en,ids,tg,s=load('astar'); i=list(tg).index(18)
    for e in [0.5,1.0,1.3715,2,5]: print(e,np.interp(e,en,s[0,i]))
