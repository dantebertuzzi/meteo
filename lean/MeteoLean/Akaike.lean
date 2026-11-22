import Mathlib

/-!
# Média por pesos de Akaike

O artigo (Seção 4.2) faz a média de `p = F(x₀)` entre as candidatas, com pesos de Akaike
`wᵢ = exp(-Δᵢ/2) / Σⱼ exp(-Δⱼ/2)`, e só depois toma `T = 1 / p̄`. Aqui:

* os pesos são positivos e somam 1 (`pesos_soma_um`, `pesos_pos`);
* `p̄` fica entre o menor e o maior `pᵢ` (`media_entre`);
* `T = 1 / p̄` nunca excede a média ponderada dos `Tᵢ = 1 / pᵢ` (`T_media_le_media_T`,
  desigualdade de Jensen para `x ↦ 1/x`). Escolher a média de `p` em vez da de `T` é, portanto,
  a escolha que dá o período de retorno **menor**: a média de Akaike do artigo não pode ser
  puxada para cima por uma candidata de cauda fina.
-/

open Finset

namespace Meteo

variable {ι : Type*} (s : Finset ι)

/-- Peso de Akaike da candidata `i`, a partir das diferenças de AIC `d`. -/
noncomputable def pesoAkaike (d : ι → ℝ) (i : ι) : ℝ :=
  Real.exp (-d i / 2) / ∑ j ∈ s, Real.exp (-d j / 2)

theorem pesos_pos {d : ι → ℝ} (hs : s.Nonempty) (i : ι) : 0 < pesoAkaike s d i :=
  div_pos (Real.exp_pos _) (sum_pos (fun _ _ => Real.exp_pos _) hs)

theorem pesos_soma_um {d : ι → ℝ} (hs : s.Nonempty) : ∑ i ∈ s, pesoAkaike s d i = 1 := by
  unfold pesoAkaike
  rw [← sum_div]
  exact div_self (sum_pos (fun j _ => Real.exp_pos _) hs).ne'

/-- Uma média ponderada (pesos ≥ 0 somando 1) fica entre o mínimo e o máximo. -/
theorem media_entre {w p : ι → ℝ} (hw : ∀ i ∈ s, 0 ≤ w i) (h1 : ∑ i ∈ s, w i = 1)
    {lo hi : ℝ} (hlo : ∀ i ∈ s, lo ≤ p i) (hhi : ∀ i ∈ s, p i ≤ hi) :
    lo ≤ ∑ i ∈ s, w i * p i ∧ ∑ i ∈ s, w i * p i ≤ hi := by
  constructor
  · calc lo = ∑ i ∈ s, w i * lo := by rw [← sum_mul, h1, one_mul]
      _ ≤ ∑ i ∈ s, w i * p i :=
          sum_le_sum fun i hi' => mul_le_mul_of_nonneg_left (hlo i hi') (hw i hi')
  · calc ∑ i ∈ s, w i * p i ≤ ∑ i ∈ s, w i * hi :=
          sum_le_sum fun i hi' => mul_le_mul_of_nonneg_left (hhi i hi') (hw i hi')
      _ = hi := by rw [← sum_mul, h1, one_mul]

/-- `1 / Σ wᵢ pᵢ ≤ Σ wᵢ / pᵢ`: o `T` da média de `p` é no máximo a média ponderada dos `T`. -/
theorem T_media_le_media_T {w p : ι → ℝ} (hw : ∀ i ∈ s, 0 ≤ w i) (h1 : ∑ i ∈ s, w i = 1)
    (hp : ∀ i ∈ s, 0 < p i) :
    (∑ i ∈ s, w i * p i)⁻¹ ≤ ∑ i ∈ s, w i * (p i)⁻¹ := by
  have hc : ConvexOn ℝ (Set.Ioi 0) (fun z : ℝ => z ^ (-1 : ℤ)) := convexOn_zpow (-1)
  have := hc.map_sum_le hw h1 (fun i hi => Set.mem_Ioi.2 (hp i hi))
  simpa [smul_eq_mul, zpow_neg_one] using this

end Meteo
