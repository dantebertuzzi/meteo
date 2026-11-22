import Mathlib
import MeteoLean.Retorno

/-!
# Cadeia de Markov de dois estados (seco / não seco)

Formaliza as quantidades derivadas em `src/paper/markov.jl`, com
`a = P(seco | não seco)` e `b = P(não seco | seco)`, de modo que `P(seco | seco) = 1 - b`:

* `π_S = a / (a + b)` é estacionária (`estacionaria`);
* a duração média de uma sequência seca é `1 / b` (`duracao_media_seca`), pelo mesmo
  argumento geométrico do tempo de espera;
* sem persistência (`P(S|S) = P(S|N)`, isto é, `1 - b = a`), a estacionária é a própria
  `P(S|N)` (`sem_persistencia`), e a persistência `P(S|S) > P(S|N)` equivale a `a + b < 1`
  (`persistencia_iff`).
-/

namespace Meteo

/-- Distribuição estacionária do estado seco. -/
noncomputable def piSeco (a b : ℝ) : ℝ := a / (a + b)

/-- `π P = π` na coordenada seca: `π_S (1 - b) + (1 - π_S) a = π_S`. -/
theorem estacionaria {a b : ℝ} (h : 0 < a + b) :
    piSeco a b * (1 - b) + (1 - piSeco a b) * a = piSeco a b := by
  unfold piSeco
  field_simp
  ring

theorem piSeco_mem {a b : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b) (h : 0 < a + b) :
    0 ≤ piSeco a b ∧ piSeco a b ≤ 1 := by
  unfold piSeco
  refine ⟨div_nonneg ha h.le, (div_le_one h).2 (by linarith)⟩

/-- Duração média de uma sequência seca: com saída `b` por ano, a sequência dura `k + 1` anos com
probabilidade `(1 - b)^k b`, e a média é `1 / b`. -/
theorem duracao_media_seca {b : ℝ} (hb0 : 0 < b) (hb1 : b ≤ 1) :
    HasSum (fun k : ℕ => ((k : ℝ) + 1) * ((1 - b) ^ k * b)) (1 / b) :=
  tempo_medio_espera hb0 hb1

/-- Sem persistência, a fração de anos secos de longo prazo é `P(S|N)`. -/
theorem sem_persistencia {a b : ℝ} (hind : 1 - b = a) : piSeco a b = a := by
  unfold piSeco
  have : a + b = 1 := by linarith
  rw [this, div_one]

/-- Persistência (`P(S|S) > P(S|N)`) equivale a `a + b < 1`. -/
theorem persistencia_iff (a b : ℝ) : a < 1 - b ↔ a + b < 1 := by
  constructor <;> intro h <;> linarith

end Meteo
