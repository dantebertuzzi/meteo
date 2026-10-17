"""Baixa as Normais Climatológicas do Brasil (INMET) e extrai a linha da 82983.

Períodos 1961–1990, 1981–2010 e 1991–2020; chuva acumulada, Tmáx e Tmín.
Fonte: https://portal.inmet.gov.br/normais (o portal recusa conexões sem
user-agent de navegador). Os arquivos originais ficam em data/raw/inmet_normais/
e a linha da estação vai para data/raw/inmet_normais/normais_82983.csv, que é o
que scripts/paper.jl lê (o Julia não lê .xls).

Uso (a partir da raiz do repositório):
    .venv-verif/bin/pip install openpyxl==3.1.5 xlrd==2.0.2
    .venv-verif/bin/python scripts/extrai_normais_inmet.py
"""
import csv, os, sys, urllib.parse, urllib.request
import openpyxl, xlrd

RAIZ = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))
DIR = os.path.join(RAIZ, "data", "raw", "inmet_normais")
BASE = "https://portal.inmet.gov.br/uploads/normais/"
UA = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/130.0 Safari/537.36"
ARQUIVOS = {
    ("1961-1990", "prec"): "Precipitacao-Acumulada_NCB_1961-1990.xls",
    ("1961-1990", "tmax"): "Temperatura-Maxima_NCB_1961-1990.xls",
    ("1961-1990", "tmin"): "Temperatura-Minima_NCB_1961-1990.xls",
    ("1981-2010", "prec"): "30-Precipitação-Acumulada-NCB_1981-2010.xls",
    ("1981-2010", "tmax"): "02-Temperatura-Máxima-NCB_1981-2010.xls",
    ("1981-2010", "tmin"): "03-Temperatura-Mínima-NCB_1981-2010.xls",
    ("1991-2020", "prec"): "Normal-Climatologica-PREC.xlsx",
    ("1991-2020", "tmax"): "Normal-Climatologica-TMAX.xlsx",
    ("1991-2020", "tmin"): "Normal-Climatologica-TMIN.xlsx",
}

def linhas(f):
    if f.endswith(".xlsx"):
        ws = openpyxl.load_workbook(f, read_only=True, data_only=True).active
        return [list(r) for r in ws.iter_rows(values_only=True)]
    sh = xlrd.open_workbook(f).sheet_by_index(0)
    return [sh.row_values(i) for i in range(sh.nrows)]

def valor(x):
    try: return float(x)
    except (TypeError, ValueError): return None   # "-" = sem normal

os.makedirs(DIR, exist_ok=True)
saida = []
for (periodo, var), nome in ARQUIVOS.items():
    f = os.path.join(DIR, nome)
    if not os.path.exists(f):
        req = urllib.request.Request(BASE + urllib.parse.quote(nome), headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=60) as r, open(f, "wb") as o:
            o.write(r.read())
    achou = [r for r in linhas(f) if r and str(r[0]).split(".")[0] == "82983"]
    if len(achou) != 1:
        sys.exit(f"{nome}: esperava uma linha da 82983, achei {len(achou)}")
    r = achou[0]   # código, nome, UF, jan..dez, ano
    for mes in range(1, 14):
        saida.append((periodo, var, mes if mes <= 12 else 0, valor(r[2 + mes])))

with open(os.path.join(DIR, "normais_82983.csv"), "w", newline="", encoding="utf-8") as o:
    w = csv.writer(o)
    w.writerow(["periodo", "variavel", "mes", "valor"])   # mes = 0: valor anual
    for p, v, m, x in saida:
        w.writerow([p, v, m, "" if x is None else x])
print(f"{len(saida)} valores em {os.path.join(DIR, 'normais_82983.csv')}")
