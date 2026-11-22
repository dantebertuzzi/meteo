import Mathlib
import MeteoLean.Retorno

/-!
# Recursão Poisson-binomial da validação fora da amostra

`prob_ao_menos_k` (em `src/paper/retorno.jl`) calcula a distribuição do número de anos secos
`N = Σ Bernoulli(pᵢ)` independentes pela recursão
`f ← [f; 0]·(1 - p) + [0; f]·p`. Aqui ela é `pb`, e provamos que:

* `pb` é uma distribuição: não negativa (`pb_nonneg`), nula acima do número de anos
  (`pb_suporte`) e com massa total 1 (`pb_soma`);
* sua média é `Σ pᵢ` (`pb_media`), o "número esperado" da Tabela 4;
* `P(N = 0) = ∏ (1 - pᵢ)` (`pb_zero`), de modo que `P(N ≥ 1)` coincide com a fórmula de
  risco com `p` variável (`pb_ao_menos_um`);
* com `p` igual em todos os anos, `pb` é a binomial (`pb_binomial`).
-/

open Finset

namespace Meteo

/-- Um passo da recursão: acrescenta um ano com probabilidade `p` de ser seco. -/
def passo (f : ℕ → ℝ) (p : ℝ) : ℕ → ℝ :=
  fun j => f j * (1 - p) + (if j = 0 then 0 else f (j - 1)) * p

/-- Distribuição de `N` para a lista de probabilidades anuais. -/
def pb : List ℝ → ℕ → ℝ
  | [] => fun j => if j = 0 then 1 else 0
  | p :: ps => passo (pb ps) p

theorem pb_suporte : ∀ (ps : List ℝ) (j : ℕ), ps.length < j → pb ps j = 0
  | [], j, h => by simp [pb]; omega
  | p :: ps, j, h => by
    simp only [pb, passo]
    have h1 : pb ps j = 0 := pb_suporte ps j (by simp at h; omega)
    have h2 : pb ps (j - 1) = 0 := pb_suporte ps (j - 1) (by simp at h; omega)
    rw [h1, h2]; simp

theorem pb_nonneg : ∀ (ps : List ℝ), (∀ p ∈ ps, 0 ≤ p ∧ p ≤ 1) → ∀ j, 0 ≤ pb ps j
  | [], _, j => by simp only [pb]; split_ifs <;> norm_num
  | p :: ps, h, j => by
    have hp := h p (by simp)
    have ih := pb_nonneg ps (fun q hq => h q (by simp [hq]))
    simp only [pb, passo]
    have a := mul_nonneg (ih j) (by linarith [hp.2] : (0 : ℝ) ≤ 1 - p)
    have b : 0 ≤ (if j = 0 then 0 else pb ps (j - 1)) * p := by
      split_ifs
      · simp
      · exact mul_nonneg (ih _) hp.1
    linarith

/-- Soma de `passo f p` sobre `range (m + 2)` quando `f` se anula a partir de `m + 1`. -/
theorem soma_passo (f : ℕ → ℝ) (p : ℝ) (m : ℕ) (hf : f (m + 1) = 0) (g : ℕ → ℝ) :
    ∑ j ∈ range (m + 2), g j * passo f p j =
      (1 - p) * ∑ j ∈ range (m + 1), g j * f j + p * ∑ j ∈ range (m + 1), g (j + 1) * f j := by
  unfold passo
  simp only [mul_add, sum_add_distrib]
  congr 1
  · rw [sum_range_succ, hf, mul_sum]; simp only [mul_zero, zero_mul, add_zero]
    exact sum_congr rfl fun j _ => by ring
  · rw [sum_range_succ', mul_sum]
    simp only [Nat.add_sub_cancel, Nat.succ_ne_zero, ↓reduceIte, zero_mul, mul_zero, add_zero]
    exact sum_congr rfl fun j _ => by ring

/-- Massa total 1. -/
theorem pb_soma : ∀ ps : List ℝ, ∑ j ∈ range (ps.length + 1), pb ps j = 1
  | [] => by simp [pb]
  | p :: ps => by
    have hf : pb ps (ps.length + 1) = 0 := pb_suporte ps _ (by omega)
    have := soma_passo (pb ps) p ps.length hf (fun _ => 1)
    simp only [one_mul] at this
    simp only [List.length_cons, pb]
    rw [this, pb_soma ps]; ring

/-- Média `Σ pᵢ`: o número esperado de anos abaixo do limiar. -/
theorem pb_media : ∀ ps : List ℝ, ∑ j ∈ range (ps.length + 1), (j : ℝ) * pb ps j = ps.sum
  | [] => by simp [pb]
  | p :: ps => by
    have hf : pb ps (ps.length + 1) = 0 := pb_suporte ps _ (by omega)
    have := soma_passo (pb ps) p ps.length hf (fun j => (j : ℝ))
    simp only [List.length_cons, pb, List.sum_cons]
    rw [this]
    have e : ∑ j ∈ range (ps.length + 1), ((j + 1 : ℕ) : ℝ) * pb ps j =
        ∑ j ∈ range (ps.length + 1), (j : ℝ) * pb ps j + ∑ j ∈ range (ps.length + 1), pb ps j := by
      rw [← sum_add_distrib]; exact sum_congr rfl fun j _ => by push_cast; ring
    rw [e, pb_media ps, pb_soma ps]; ring

/-- `P(N = 0) = ∏ (1 - pᵢ)`. -/
theorem pb_zero : ∀ ps : List ℝ, pb ps 0 = (ps.map (fun p => 1 - p)).prod
  | [] => by simp [pb]
  | p :: ps => by simp [pb, passo, pb_zero ps]; ring

/-- `P(N ≥ 1) = 1 - P(N = 0)` coincide com o risco em `n` anos de `Retorno.lean`. -/
theorem pb_ao_menos_um (p : ℕ → ℝ) (n : ℕ) :
    1 - pb ((List.range n).map p) 0 = probVar p n := by
  rw [pb_zero, probVar, List.map_map]
  congr 1

/-- Com `p` igual em todos os anos, a recursão dá a binomial. -/
theorem pb_binomial (p : ℝ) : ∀ (n j : ℕ),
    pb (List.replicate n p) j = (n.choose j : ℝ) * p ^ j * (1 - p) ^ (n - j)
  | 0, j => by
    cases j <;> simp [pb]
  | n + 1, j => by
    simp only [List.replicate_succ, pb, passo, pb_binomial p n]
    rcases j with _ | j
    · simp; ring
    · simp only [Nat.succ_ne_zero, ↓reduceIte, Nat.add_sub_cancel, Nat.choose_succ_succ, Nat.cast_add]
      rcases lt_trichotomy j n with h | rfl | h
      · obtain ⟨k, rfl⟩ := Nat.exists_eq_add_of_lt h
        have e1 : j + k + 1 - (j + 1) = k := by omega
        have e2 : j + k + 1 + 1 - (j + 1) = k + 1 := by omega
        have e3 : j + k + 1 - j = k + 1 := by omega
        rw [e1, e2, e3]; ring
      · simp [Nat.choose_succ_self]; ring
      · have h1 : n.choose (j + 1) = 0 := Nat.choose_eq_zero_of_lt (by omega)
        have h2 : n.choose j = 0 := Nat.choose_eq_zero_of_lt h
        simp [h1, h2]

end Meteo
