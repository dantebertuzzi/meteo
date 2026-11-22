import Mathlib

/-!
# Período de retorno e probabilidade em `n` anos

Formaliza as fórmulas fechadas de `src/paper/retorno.jl` e da Seção 4.2 do artigo:

* `probAoMenosUm p n = 1 - (1 - p)^n` (Julia: `prob_ao_menos_um`, Equação 4);
* a soma das probabilidades de o *primeiro* ano abaixo do limiar ser o ano `k + 1` é exatamente
  essa probabilidade (`soma_primeiro_evento`);
* o tempo médio de espera até o primeiro evento, com anos independentes e `p` fixo, é `1 / p`
  (`tempo_medio_espera`): é isso que justifica chamar `T = 1/p` de período de retorno;
* com `p` variável ano a ano, `1 - ∏ (1 - pⱼ)` (Julia: `prob_n_variavel`) reduz-se ao caso fixo
  quando `p` é constante e é monótona em cada `pⱼ` (`probVar_mono`). Daí: se a tendência faz
  `p` crescer, o cenário "tendência continua" tem risco ≥ ao do "clima congelado"
  (`continua_ge_congelado`), como na Tabela 6.
-/

open Finset

namespace Meteo

/-- Probabilidade de ao menos um ano com chuva ≤ x₀ em `n` anos independentes, com `p = F(x₀)`. -/
noncomputable def probAoMenosUm (p : ℝ) (n : ℕ) : ℝ := 1 - (1 - p) ^ n

/-- A probabilidade de o primeiro evento cair no ano `k + 1` é `(1 - p)^k · p`; somando de
`k = 0` a `n - 1` obtém-se a probabilidade de ao menos um evento em `n` anos. -/
theorem soma_primeiro_evento (p : ℝ) (n : ℕ) :
    ∑ k ∈ range n, (1 - p) ^ k * p = probAoMenosUm p n := by
  unfold probAoMenosUm
  induction n with
  | zero => simp
  | succ n ih => rw [sum_range_succ, ih, pow_succ]; ring

theorem probAoMenosUm_nonneg {p : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) (n : ℕ) :
    0 ≤ probAoMenosUm p n := by
  unfold probAoMenosUm
  have : (1 - p) ^ n ≤ 1 := pow_le_one₀ (by linarith) (by linarith)
  linarith

theorem probAoMenosUm_le_one {p : ℝ} (hp1 : p ≤ 1) (n : ℕ) : probAoMenosUm p n ≤ 1 := by
  unfold probAoMenosUm
  have : 0 ≤ (1 - p) ^ n := pow_nonneg (by linarith) n
  linarith

/-- Mais anos, mais risco. -/
theorem probAoMenosUm_mono_n {p : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) {m n : ℕ} (h : m ≤ n) :
    probAoMenosUm p m ≤ probAoMenosUm p n := by
  unfold probAoMenosUm
  have : (1 - p) ^ n ≤ (1 - p) ^ m := pow_le_pow_of_le_one (by linarith) (by linarith) h
  linarith

/-- Tempo médio de espera: `E[K] = Σ_{k≥0} (k + 1)(1 - p)^k p = 1/p` para `0 < p ≤ 1`. -/
theorem tempo_medio_espera {p : ℝ} (hp0 : 0 < p) (hp1 : p ≤ 1) :
    HasSum (fun k : ℕ => ((k : ℝ) + 1) * ((1 - p) ^ k * p)) (1 / p) := by
  have hr : ‖(1 - p)‖ < 1 := by
    rw [Real.norm_eq_abs, abs_lt]; constructor <;> linarith
  have h1 := hasSum_coe_mul_geometric_of_norm_lt_one hr
  have h2 := hasSum_geometric_of_norm_lt_one hr
  have h3 := (h1.add h2).mul_right p
  have e : (1 : ℝ) - (1 - p) = p := by ring
  rw [e] at h3
  convert h3 using 1
  · funext k; ring
  · field_simp; ring

/-- Risco em `n` anos quando `p` varia ano a ano: `1 - ∏_{j<n} (1 - pⱼ)`. -/
noncomputable def probVar (p : ℕ → ℝ) (n : ℕ) : ℝ := 1 - ∏ j ∈ range n, (1 - p j)

/-- Com `p` constante, a fórmula variável é a da Equação 4. -/
theorem probVar_const (p : ℝ) (n : ℕ) : probVar (fun _ => p) n = probAoMenosUm p n := by
  simp [probVar, probAoMenosUm, prod_const, card_range]

/-- Monotonia: se `pⱼ ≤ qⱼ ≤ 1` em todo ano, o risco com `q` é ao menos o risco com `p`. -/
theorem probVar_mono {p q : ℕ → ℝ} (hq : ∀ j, q j ≤ 1)
    (hpq : ∀ j, p j ≤ q j) (n : ℕ) : probVar p n ≤ probVar q n := by
  unfold probVar
  have : ∏ j ∈ range n, (1 - q j) ≤ ∏ j ∈ range n, (1 - p j) := by
    apply prod_le_prod₀
    · intro j _; linarith [hq j]
    · intro j _; linarith [hpq j]
  linarith

/-- Tabela 6: se, a partir do ano de referência, a probabilidade anual só cresce (`p₀ ≤ pⱼ`), o
cenário em que a tendência continua tem risco ao menos igual ao do clima congelado em `p₀`. -/
theorem continua_ge_congelado {p : ℕ → ℝ} {p₀ : ℝ} (h1 : ∀ j, p j ≤ 1)
    (hcresce : ∀ j, p₀ ≤ p j) (n : ℕ) : probAoMenosUm p₀ n ≤ probVar p n := by
  rw [← probVar_const]
  exact probVar_mono h1 hcresce n

end Meteo
