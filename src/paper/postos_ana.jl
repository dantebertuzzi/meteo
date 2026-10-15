# =============================================================================
# postos_ana.jl — pluviômetros da rede da ANA usados como verificação independente
#
# Candidatos: postos de operadores diferentes do INMET a até cerca de 70 km da 82983, com
# registro longo, do inventário HidroInventario da ANA (consultado em 15/10/2026).
# Coordenadas e operador vêm do inventário. A seleção final (≥ 30 anos
# hidrológicos válidos em comum com a 82983) é feita em scripts/paper.jl. `chave` nomeia
# as macros do artigo.
# =============================================================================

const POSTOS_ANA = Dict(
    "940024" => (chave = "Juazeiro", nome = "Juazeiro",          resp = "ANA",     lat = -9.4053, lon = -40.5036),
    "940016" => (chave = "JuazeiroDnocs", nome = "Juazeiro (DNOCS)",  resp = "DNOCS",   lat = -9.4167, lon = -40.5),
    "940006" => (chave = "PetrolinaApac", nome = "Petrolina (APAC)",  resp = "APAC-PE", lat = -9.3833, lon = -40.5),
    "940018" => (chave = "Cavalos", nome = "Campo dos Cavalos", resp = "ANA",     lat = -9.5597, lon = -40.6467),
    "940028" => (chave = "Junco", nome = "Junco",             resp = "ANA",     lat = -9.6767, lon = -40.6039),
    "940025" => (chave = "LagoaBoi", nome = "Lagoa do Boi",      resp = "ANA",     lat = -9.4908, lon = -40.2044),
    "940020" => (chave = "Juremal", nome = "Juremal",           resp = "DNOCS",   lat = -9.7286, lon = -40.3578),
    "840015" => (chave = "LagoaGrande", nome = "Lagoa Grande",      resp = "ANA",     lat = -8.9989, lon = -40.2731),
    "840017" => (chave = "PauFerro", nome = "Pau Ferro",         resp = "DNOCS",   lat = -8.95,   lon = -40.7333),
    "941011" => (chave = "CasaNova", nome = "Casa Nova",         resp = "DNOCS",   lat = -9.175,  lon = -40.97),
    "940005" => (chave = "MalhadaReal", nome = "Malhada Real",      resp = "DNOCS",   lat = -9.0333, lon = -40.0167),
)

"Distância aproximada (km) entre dois pontos (equirretangular; suficiente para < 100 km)."
dist_km(lat1, lon1, lat2, lon2) = 111.2 * hypot(lat1 - lat2, (lon1 - lon2) * cosd((lat1 + lat2) / 2))
