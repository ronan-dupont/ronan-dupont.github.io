---
title: "ICMASE 2026 — GPBiCG(m, ℓ, k) and adaptive parameter strategies"
title_fr: "ICMASE 2026 — GPBiCG(m, ℓ, k) et stratégies adaptatives de paramètres"
title_ja: "ICMASE 2026 — GPBiCG(m, ℓ, k)と適応的パラメータ戦略"
collection: talks
type: "Conference presentation"
type_fr: "Présentation en conférence"
type_ja: "学会発表"
permalink: /talks/ICMASE2026/
venue: "ICMASE 2026"
location: "Nagoya, Japan"
location_fr: "Nagoya, Japon"
location_ja: "名古屋, 日本"
date: 2026-09-14
cover: /images/talks/icmase2026-cover.jpg
slides_handout: "/files/presentation/ICMASE2026_ultra_handout.pdf?v=20"
slides_progressive: "/files/presentation/ICMASE2026_ultra_progressive.pdf?v=20"
poster: "/files/poster/Moonshot_GPBiCG_MLK_A0_portrait.pdf?v=19"
poster_label: "Moonshot poster · A0"
poster_label_fr: "Poster Moonshot · A0"
poster_label_ja: "Moonshotポスター · A0"
excerpt: "An extension of GPBiCG(m, ℓ) and preliminary adaptive parameter strategies for nonsymmetric linear systems."
excerpt_fr: "Une extension de GPBiCG(m, ℓ) et des stratégies adaptatives préliminaires pour les systèmes linéaires non symétriques."
excerpt_ja: "非対称線形方程式系に対するGPBiCG(m, ℓ)の拡張と予備的な適応パラメータ戦略。"
---

My presentation at ICMASE 2026 in Nagoya, Japan, on September 14, 2026, focuses on **GPBiCG(m, ℓ, k): An Extension of GPBiCG(m, ℓ) and Preliminary Adaptive Parameter Strategies**. This research is carried out in the Zhang–Sogabe Laboratory at Nagoya University, under the supervision of Tomohiro Sogabe and Shao-Liang Zhang.

The talk focuses on iterative methods for solving nonsymmetric linear systems. GPBiCG(m, ℓ, k) extends GPBiCG(m, ℓ) by adding a third phase based on an orthogonality condition. This gives the solver more ways to reduce the residual, but also introduces another parameter to choose.

Our numerical study compares the method on eight matrices, with and without ILU(0) preconditioning. The additional phase can improve convergence and CPU time, especially without preconditioning, although the best parameter triplet depends on the matrix. The talk also introduces a preliminary controller that adjusts the triplet during a solve. It can rescue some poor initial choices, while monitoring costs and the choice of a starting strategy remain important limitations. A fully automatic method that consistently outperforms the classical baselines is still a research direction.

The appendix includes the complete GPBiCG algorithm and the phase rules for GPBiCG(m, ℓ, k).

You can find the presentation and associated Moonshot poster here:

- [Slides — handout PDF]({{ '/files/presentation/ICMASE2026_ultra_handout.pdf' | relative_url }}?v=20)
- [Slides — progressive presentation PDF]({{ '/files/presentation/ICMASE2026_ultra_progressive.pdf' | relative_url }}?v=20)
- [Associated Moonshot poster — A0 portrait PDF]({{ '/files/poster/Moonshot_GPBiCG_MLK_A0_portrait.pdf' | relative_url }}?v=19)

[![Title slide of Ronan Dupont’s ICMASE 2026 presentation on GPBiCG(m, ℓ, k)]({{ '/images/talks/icmase2026-cover.jpg' | relative_url }})]({{ '/files/presentation/ICMASE2026_ultra_handout.pdf' | relative_url }}?v=20)
