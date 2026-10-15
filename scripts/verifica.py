"""Verificação independente (Python/SciPy) dos números do artigo.

Lê os dados BRUTOS e reimplementa os cálculos sem usar o código Julia
(SciPy para ajustes e Theil–Sen, pymannkendall para MK/TFPW), e compara com as
macros de paper/resultados.tex. Verifica também a consistência interna das macros
(estimativa dentro do IC; probabilidades crescentes com o horizonte).

Uso (a partir da raiz do repositório, depois de `julia --project=. scripts/paper.jl`):
    python3 -m venv .venv-verif
    .venv-verif/bin/pip install -r scripts/requirements-verificacao.txt
    .venv-verif/bin/python scripts/verifica.py
Sai com código 1 se algum valor não conferir.
"""
import re, sys, math, collections, datetime as dt
import numpy as np, pandas as pd
from scipy import stats, optimize
import pymannkendall as mk

RAIZ = sys.argv[1] if len(sys.argv) > 1 else __import__("os").path.normpath(__import__("os").path.join(__import__("os").path.dirname(__file__), ".."))
M = {}
for l in open(f"{RAIZ}/paper/resultados.tex", encoding="utf-8"):
    m = re.match(r"\\newcommand\{\\(\w+)\}\{(.*)\}$", l.strip())
    if m:
        v = m.group(2); n = re.search(r"\\num\{([^}]*)\}", v)
        M[m.group(1)] = n.group(1) if n else v

def mac(k):
    s = M[k].replace("<", "").strip()
    try: return float(s)
    except: return M[k]

res = []
def compara(nome, indep, macro, tol_rel=0.02, tol_abs=None):
    mv = mac(macro)
    if isinstance(mv, str):
        ok = str(indep) == mv
    else:
        tol = tol_abs if tol_abs is not None else max(abs(mv) * tol_rel, 10 ** -(len(M[macro].split(".")[-1]) if "." in M[macro] else 0) * 0.51)
        ok = abs(indep - mv) <= tol
    res.append((ok, nome, indep, M[macro], macro))

# ---------------------------------------------------------------- 1. série diária 82983
def le_bdmep(f):
    linhas = open(f, encoding="utf-8").read().splitlines()
    i = next(k for k, l in enumerate(linhas) if l.startswith("Data Medicao"))
    cab = linhas[i].split(";")
    rows = [l.split(";") for l in linhas[i + 1:] if l.strip()]
    df = pd.DataFrame(rows, columns=cab[:len(rows[0])])
    df = df.rename(columns={cab[0]: "data"})
    df["data"] = pd.to_datetime(df["data"])
    def col(prefix):
        c = next(c for c in cab if c.upper().startswith(prefix))
        return pd.to_numeric(df[c].str.replace(",", "."), errors="coerce")
    out = pd.DataFrame({"data": df["data"], "prec": col("PRECIPITACAO"), "tmax": col("TEMPERATURA MAXIMA"),
                        "tmin": col("TEMPERATURA MINIMA")})
    return out

d = le_bdmep(f"{RAIZ}/data/raw/bdmep/dados_82983_D_1961-01-01_2025-06-30.csv")
n_gt150 = int((d.prec > 150).sum())
print(f"[dados] dias com chuva > 150 mm na 82983 bruta: {n_gt150}; máximo diário = {d.prec.max():.1f} mm")
d["ano"] = d.data.dt.year; d["mes"] = d.data.dt.month
d["ah"] = np.where(d.mes >= 7, d.ano, d.ano - 1)

ref = d[(d.ano >= 1991) & (d.ano <= 2020) & d.prec.notna()]
cd = ref.groupby("mes").prec.mean().to_dict()
dias_mes = {m: (pd.Timestamp(2001, m, 1) + pd.offsets.MonthEnd(0)).day for m in range(1, 13)}
normal = sum(cd[m] * dias_mes[m] for m in range(1, 13))
compara("normal 1991–2020 (mm)", normal, "NormalChuva", tol_abs=0.6)

pos = dict(zip(d.data, d.prec))
hid = []
for a in range(1960, 2025):
    dias = pd.date_range(f"{a}-07-01", f"{a+1}-06-30")
    aus = sum(cd[x.month] for x in dias if pd.isna(pos.get(x, np.nan)))
    tot = np.nansum([pos.get(x, np.nan) for x in dias])
    hid.append((a, tot, aus / normal))
H = pd.DataFrame(hid, columns=["ah", "tot", "fca"])
V = H[H.fca <= 0.05].sort_values("ah")
compara("nº de anos hidrológicos válidos", len(V), "NanosValidos", tol_abs=0)
compara("chuva 2012/13 (mm)", float(V[V.ah == 2012].tot.iloc[0]), "ChuvaDozeTreze", tol_abs=0.06)
compara("chuva mínima (mm)", V.tot.min(), "ChuvaMinima", tol_abs=0.06)
compara("ano mais seco", f"{int(V.loc[V.tot.idxmin(), 'ah'])}/{(int(V.loc[V.tot.idxmin(), 'ah'])+1)%100:02d}", "AnoMinimo")
x = V.tot.values; anos = V.ah.values

# ---------------------------------------------------------------- 2. distribuições
x0 = float(V[V.ah == 2012].tot.iloc[0])
fits = {
    "Gamma": stats.gamma(*stats.gamma.fit(x, floc=0)),
    "Weibull": stats.weibull_min(*stats.weibull_min.fit(x, floc=0)),
    "LogNormal": stats.lognorm(*stats.lognorm.fit(x, floc=0)),
    "Normal": stats.norm(*stats.norm.fit(x)),
}
aic = {k: 2 * 2 - 2 * np.sum(f.logpdf(x)) for k, f in fits.items()}
amin = min(aic.values()); w = {k: math.exp(-(a - amin) / 2) for k, a in aic.items()}; s = sum(w.values()); w = {k: v / s for k, v in w.items()}
for k, f in fits.items():
    compara(f"T {k} (anos)", 1 / f.cdf(x0), "T" + k, tol_rel=0.05)
compara("melhor AIC", min(aic, key=aic.get), "MelhorAIC")
compara("peso de Akaike da Gamma", w["Gamma"], "GammaPesoAkaike", tol_abs=0.006)
pavg = sum(w[k] * fits[k].cdf(x0) for k in fits)
compara("T média por AIC (amostra completa)", 1 / pavg, "TMediaAkaike", tol_rel=0.05)
print(f"[AIC] T da média ponderada na amostra completa = {1/pavg:.1f}; macro TMediaAkaikeB (usada no texto) = {M['TMediaAkaikeB']}")

# perfil da Gamma (independente: otimiza k para cada p numa grade)
def ll_kp(k, p):
    th = x0 / stats.gamma.ppf(p, k)
    return np.sum(stats.gamma.logpdf(x, k, scale=th))
def lp(p):
    r = optimize.minimize_scalar(lambda lk: -ll_kp(math.exp(lk), p), bounds=(math.log(0.05), math.log(200)), method="bounded")
    return -r.fun
g = fits["Gamma"]; lmax = np.sum(g.logpdf(x)); ph = g.cdf(x0)
grid = np.exp(np.linspace(math.log(1e-4), math.log(0.2), 400))
dev = np.array([2 * (lmax - lp(p)) for p in grid])
dentro = grid[dev <= stats.chi2.ppf(0.95, 1)]
compara("T Gamma perfilada (est.)", 1 / ph, "TDozeTreze", tol_rel=0.03)
compara("T Gamma perfilada (IC inf)", 1 / dentro.max(), "TDozeTrezeICinf", tol_rel=0.05)
compara("T Gamma perfilada (IC sup)", 1 / dentro.min(), "TDozeTrezeICsup", tol_rel=0.06)

# ---------------------------------------------------------------- 3. tendência da chuva
r = mk.original_test(x); ts = stats.theilslopes(x, anos, alpha=0.95)
compara("Sen chuva (mm/década)", 10 * ts.slope, "SenChuva", tol_abs=0.6)
compara("Sen chuva IC inf", 10 * ts.low_slope, "SenChuvaICinf", tol_abs=1.6)
compara("Sen chuva IC sup", 10 * ts.high_slope, "SenChuvaICsup", tol_abs=1.6)
for nome_, m_, k_ in (("sem 2011/12–2016/17", ~((anos >= 2011) & (anos <= 2016)), "SenChuvaSemSeca"),
                     ("até 2010/11", anos <= 2010, "SenChuvaAteDezOnze"), ("desde 1971/72", anos >= 1971, "SenChuvaDesdeSetentaUm")):
    t_ = stats.theilslopes(x[m_], anos[m_], alpha=0.95)
    compara(f"Sen chuva {nome_}", 10 * t_.slope, k_, tol_abs=0.6)
print(f"[MK] chuva: p MK = {r.p:.2e}; p TFPW (pymannkendall) = {mk.trend_free_pre_whitening_modification_test(x).p:.4f}; macro PTFPWChuva = {M['PTFPWChuva']}")

# ---------------------------------------------------------------- 4. Gamma não estacionária
t = (anos - anos.mean()) / 10
def nll(th):
    k = math.exp(th[0]); mu = np.exp(th[1] + th[2] * t)
    return -np.sum(stats.gamma.logpdf(x, k, scale=mu / k))
o = optimize.minimize(nll, [math.log(g.args[0]), math.log(x.mean()), 0.0], method="Nelder-Mead", options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 20000})
k1, a1, b1 = math.exp(o.x[0]), o.x[1], o.x[2]
compara("efeito NS (%/década)", 100 * (math.exp(b1) - 1), "NSTempoEfeito", tol_abs=0.15)
LR = 2 * (-o.fun - lmax)
compara("ΔAIC NS − estacionário", (6 + 2 * o.fun) - (4 - 2 * lmax), "NSTempoDAIC", tol_abs=0.15)
tu = (2020 - anos.mean()) / 10
Tult = 1 / stats.gamma.cdf(x0, k1, scale=math.exp(a1 + b1 * tu) / k1)
compara("T NS MV em 2020/21", Tult, "TTempoUlt", tol_rel=0.05)

# média por pesos de Akaike das candidatas não estacionárias: posição log-linear e forma
# fixa ou log-linear no tempo (log s = s0 + c·t)
def ns_fit(nome, var=False):
    def mk_(th, tt):
        s_ = np.exp(th[0] + (th[3] * tt if var else 0.0)); eta = th[1] + th[2] * tt
        return {"Normal": lambda: stats.norm(np.exp(eta), s_), "LogNormal": lambda: stats.lognorm(s_, scale=np.exp(eta)),
                "Gamma": lambda: stats.gamma(s_, scale=np.exp(eta) / s_), "Weibull": lambda: stats.weibull_min(s_, scale=np.exp(eta))}[nome]()
    f0 = fits[nome]
    s0 = {"Normal": x.std(), "LogNormal": np.log(x).std(), "Gamma": f0.args[0], "Weibull": f0.args[0]}[nome]
    a0 = {"Normal": math.log(x.mean()), "LogNormal": np.log(x).mean(), "Gamma": math.log(x.mean()), "Weibull": math.log(f0.kwds.get("scale", f0.args[-1]))}[nome]
    nl = lambda th: (lambda v: v if np.isfinite(v) else 1e12)(-np.sum(mk_(th, t).logpdf(x)))
    th0 = [math.log(s0), a0, 0.0]
    if var:
        th0 = list(ns_fit(nome)[2]) + [0.0]
    o_ = optimize.minimize(nl, th0, method="Nelder-Mead", options={"xatol": 1e-9, "fatol": 1e-10, "maxiter": 40000})
    o_ = optimize.minimize(nl, o_.x, method="Nelder-Mead", options={"xatol": 1e-10, "fatol": 1e-11, "maxiter": 40000})
    return 2 * len(th0) + 2 * o_.fun, float(mk_(o_.x, tu).cdf(x0)), o_.x, -o_.fun
FAM = ("Normal", "Weibull", "Gamma", "LogNormal")
ns = {(k, v): ns_fit(k, v) for k in FAM for v in (False, True)}
def media(chaves):
    amin_ = min(ns[c][0] for c in chaves)
    w_ = {c: math.exp(-(ns[c][0] - amin_) / 2) for c in chaves}; sw_ = sum(w_.values())
    return w_, sw_, sw_ / sum(w_[c] * ns[c][1] for c in chaves)
for k in ("Normal", "Weibull", "Gamma"):
    compara(f"T NS {k} em 2020/21", 1 / ns[(k, False)][1], "TNS" + k, tol_rel=0.05)
    compara(f"T NS {k} forma variável", 1 / ns[(k, True)][1], "TNS" + k + "Var", tol_rel=0.05)
compara("LR forma variável (Gamma) p", stats.chi2.sf(2 * (ns[("Gamma", True)][3] - ns[("Gamma", False)][3]), 1), "PFormaVarGamma", tol_abs=0.003)
_, _, Tfix = media([(k, False) for k in FAM])
compara("T NS média, só forma fixa", Tfix, "TNSMediaFixaUlt", tol_rel=0.05)
w8, sw8, T8 = media(list(ns))
compara("peso de Akaike NS da normal fixa", w8[("Normal", False)] / sw8, "PesoNSNormal", tol_abs=0.006)
compara("T NS média das 8 candidatas", T8, "TNSMediaUlt", tol_rel=0.05)

# validação fora da amostra: Gamma estacionária ajustada até 2005/06, limiar P10 do ajuste
trv = anos <= 2005; xtr, xte = x[trv], x[~trv]
p10 = np.quantile(xtr, 0.10)            # mesma definição de quantil (tipo 7) do Julia
compara("limiar P10 da validação (mm)", p10, "LimiarPdez", tol_abs=0.6)
gtr = stats.gamma(*stats.gamma.fit(xtr, floc=0))
compara("esperado P10, Gamma estacionária", len(xte) * gtr.cdf(p10), "EspPdezEstGamma", tol_abs=0.011)
compara("observado P10 no teste", int((xte <= p10).sum()), "ObsPdez", tol_abs=0)

# ---------------------------------------------------------------- 5. SPI-12 de junho (2012/13, 2016/17)
mens = d.groupby(["ano", "mes"]).agg(prec=("prec", "sum"), aus=("prec", lambda s: s.isna().sum()), n=("prec", "size")).reset_index()
mens["dias"] = [dias_mes[m] + (1 if m == 2 and (a % 4 == 0 and (a % 100 != 0 or a % 400 == 0)) else 0) for a, m in zip(mens.ano, mens.mes)]
mens["aus"] = mens.aus + mens.dias - mens.n
mens["esp"] = [aus * cd[m] for aus, m in zip(mens.aus, mens.mes)]
mens["valido"] = (mens.aus < mens.dias) & (mens.esp <= np.maximum(1.0, 0.05 * np.array([cd[m] * dias_mes[m] for m in mens.mes])))
mens = mens.set_index(["ano", "mes"])
def soma12(a):  # jul/(a-1) .. jun/a
    ks = [(a - 1, m) for m in range(7, 13)] + [(a, m) for m in range(1, 7)]
    if not all(k in mens.index and mens.loc[k, "valido"] for k in ks): return None
    return sum(mens.loc[k, "prec"] for k in ks)
cal = [v for a in range(1991, 2021) if (v := soma12(a)) is not None]
q = 1 - sum(1 for v in cal if v > 0) / len(cal)
gs = stats.gamma(*stats.gamma.fit([v for v in cal if v > 0], floc=0))
for a, mk_ in ((2013, "SPIDozeTreze"), (2017, "SPIDezesseisDezessete")):
    z = stats.norm.ppf(q + (1 - q) * gs.cdf(soma12(a)))
    compara(f"SPI-12 jun/{a}", z, mk_, tol_abs=0.03)

# ---------------------------------------------------------------- 6. Markov (tercil inferior bruto)
corte = np.quantile(x, 1 / 3)
est = {int(a): ("S" if v <= corte else "N") for a, v in zip(anos, x)}
C = collections.Counter((est[a], est[a + 1]) for a in est if a + 1 in est)
for par, mk_ in ((("S", "S"), "NSS"), (("S", "N"), "NSN"), (("N", "S"), "NNS"), (("N", "N"), "NNN")):
    compara(f"transições {par}", C[par], mk_, tol_abs=0)
compara("corte do tercil (mm)", corte, "CorteSeco", tol_abs=0.6)

# ---------------------------------------------------------------- 7. tendências de temperatura (série anual do Julia, teste independente)
an = pd.read_csv(f"{RAIZ}/results/paper/anual_civil_82983.csv")
for col, k in (("tmin", "Tmin"), ("tmax", "Tmax")):
    s = an[["ano", col]].dropna()
    ts = stats.theilslopes(s[col], s.ano, alpha=0.95)
    compara(f"Sen {k} (°C/década)", 10 * ts.slope, "Sen" + k, tol_abs=0.006)
    s75 = s[s.ano >= 1975]; t75 = stats.theilslopes(s75[col], s75.ano, alpha=0.95)
    compara(f"Sen {k} desde 1975", 10 * t75.slope, "Sen" + k + "DesdeSetentaCinco", tol_abs=0.006)
# Tmáx sem o salto frente a Cabrobó (ano e salto do SNHT do artigo; a checagem é do ajuste e do Sen)
s = an[["ano", "tmax"]].dropna()
a0, sal = int(float(mac("QuebraTmaxCabrobo"))), float(mac("SaltoTmaxCabrobo"))
ts = stats.theilslopes(s.tmax - sal * (s.ano >= a0), s.ano, alpha=0.95)
compara("Sen Tmáx ajustada", 10 * ts.slope, "SenTmaxAjustada", tol_abs=0.011)
# Tmín anual recalculada do bruto (média dos dias, anos com ≥ 300 dias, ≥ 1971) — checagem da agregação
dd = d[d.ano >= 1971]
ta = dd.groupby("ano").tmin.agg(["mean", "count"]); ta = ta[ta["count"] >= 300]
ts = stats.theilslopes(ta["mean"], ta.index, alpha=0.95)
print(f"[temp] Tmín anual pela média simples dos dias (≥300 dias): Sen = {10*ts.slope:.3f} °C/década (macro {M['SenTmin']}) — método diferente, deve ser próximo")

# ---------------------------------------------------------------- 8. teleconexão ONI (bruta)
def psl(nome):
    tab = {}
    for l in open(f"{RAIZ}/data/raw/indices/{nome}.data").read().splitlines()[1:]:
        c = l.split()
        if len(c) == 13 and c[0].isdigit():
            for m in range(12):
                v = float(c[m + 1])
                if v > -99: tab[(int(c[0]), m + 1)] = v
    return tab
oni = psl("oni")
def novabr(tab, a):
    ks = [(a, 11), (a, 12), (a + 1, 1), (a + 1, 2), (a + 1, 3), (a + 1, 4)]
    return np.mean([tab[k] for k in ks]) if all(k in tab for k in ks) else np.nan
o_ = np.array([novabr(oni, int(a)) for a in anos]); ok = ~np.isnan(o_)
compara("r(chuva, ONI)", stats.pearsonr(x[ok], o_[ok])[0], "ROniBruta", tol_abs=0.006)


# ---------------------------------------------------------------- 10. referência regional (Alexandersson & Moberg)
import glob, os
def hidro_validos(df):
    df = df.copy(); df["ano"] = df.data.dt.year; df["mes"] = df.data.dt.month
    r = df[(df.ano >= 1991) & (df.ano <= 2020) & df.prec.notna()]
    if r.data.dt.date.nunique() < 3650: r = df[df.prec.notna()]
    c = r.groupby("mes").prec.mean().to_dict()
    if len(c) < 12: return {}
    nm = sum(c[m] * dias_mes[m] for m in range(1, 13)); pp = dict(zip(df.data, df.prec)); out = {}
    for a in range(1960, 2025):
        dias = pd.date_range(f"{a}-07-01", f"{a+1}-06-30")
        aus = sum(c[x.month] for x in dias if pd.isna(pp.get(x, np.nan)))
        if aus / nm <= 0.05: out[a] = np.nansum([pp.get(x, np.nan) for x in dias])
    return out
cand = dict(zip(V.ah.astype(int), V.tot))
viz = {}
for f in glob.glob(f"{RAIZ}/data/raw/bdmep/dados_8*_D_*.csv"):
    cod = os.path.basename(f).split("_")[1]
    if cod == "82983": continue
    viz[cod] = hidro_validos(le_bdmep(f))
pesos = {}
for cod, vz in viz.items():
    com = [a for a in cand if a in vz]
    dif = [(cand[a+1]-cand[a], vz[a+1]-vz[a]) for a in com if a+1 in cand and a+1 in vz]
    if len(dif) >= 8:
        rho = np.corrcoef(*zip(*dif))[0, 1]
        if len(com) >= 15 and rho >= 0.5: pesos[cod] = rho ** 2
esc = {c: np.mean([cand[a] for a in cand if a in viz[c]]) / np.mean([viz[c][a] for a in cand if a in viz[c]]) for c in pesos}
ref = {}
for a in sorted(cand):
    ns = [c for c in pesos if a in viz[c]]
    if ns: ref[a] = sum(pesos[c] * viz[c][a] * esc[c] for c in ns) / sum(pesos[c] for c in ns)
ra = np.array(sorted(ref)); rv = np.array([ref[a] for a in ra])
compara("nº de vizinhas usadas", len(pesos), "NvizinhasUsadas", tol_abs=0)
compara("Sen da referência regional", 10 * stats.theilslopes(rv, ra).slope, "SenReferencia", tol_abs=0.6)
qv = np.array([cand[a] / ref[a] for a in ra])
compara("Sen da razão Q (%/década)", 1000 * stats.theilslopes(qv, ra).slope, "SenQ", tol_abs=0.06)

# ---------------------------------------------------------------- 9. consistência interna das macros
inc = []
for k in M:
    if k.endswith("ICinf") and k[:-5] in M and k[:-5] + "ICsup" in M:
        try:
            e, lo, hi = mac(k[:-5]), mac(k), mac(k[:-5] + "ICsup")
            if isinstance(e, float) and not (lo - 1e-9 <= e <= hi + 1e-9):
                inc.append(f"{k[:-5]}: {e} fora de [{lo}; {hi}]")
        except Exception: pass
for base in ("ProbDozeTreze", "ProbNSUlt", "ProbNSCong", "ProbNSCont", "ProbEst", "ProbReg",
             "ProbMediaAkaike", "ProbNSMediaUlt", "ProbRegMedia"):
    vals = [mac(base + s) for s in ("Ncinco", "Numzero", "Ndoiszero", "Ntreszero") if base + s in M]
    if vals != sorted(vals): inc.append(f"{base}: probabilidades não crescem com n: {vals}")

# ---------------------------------------------------------------- relatório
print("\n=== Comparação independente (Python/SciPy × macros do artigo) ===")
for ok, nome, iv, mv, k in res:
    ivs = f"{iv:.4g}" if isinstance(iv, (float, np.floating)) else str(iv)
    print(f"{'OK ' if ok else 'DIF'} {nome:38s} indep = {ivs:>12s} | artigo = {mv:>10s}  ({k})")
print(f"\n{sum(r[0] for r in res)}/{len(res)} conferem.")
print("\n=== Consistência interna das macros ===")
print("\n".join(inc) if inc else "Nenhuma inconsistência (toda estimativa dentro do IC; probabilidades crescem com n).")
sys.exit(0 if all(r[0] for r in res) and not inc else 1)
