---
title: "PhD Defense"
collection: talks
type: "Presentation"
permalink: /talks/phd_defense
venue: "PhD Defense"
date: 2024-09-30
location: "Montpellier, France"
---

On September 30, 2024, I defended my thesis in Montpellier.

Below is a summary of the work: 

Morphodynamical models in shallow coastal waters is a challenging topic, especially when trying to reproduce physical phenomena such as sandbar creation. Classic models are generally very complex and highly parameterized; they separately solve the physical equations of hydrodynamics and morphodynamics at a very small scale of the order of second in time and meter in  space. 
During this thesis, we developed a numerical model proposing a more global approach to coastal morphodynamics, based on an optimization principle.  

The optimization theory is the study of the evolution of a system while searching systematically for the minimum of a function derived from some of its physical properties. Using mathematical optimization theory, we have designed a model that describes the evolution of the sea bottom elevation while taking into account the coupling between morphodynamic and hydrodynamic processes. Our model is based on the assumption that the sea bottom adapts to minimize a wave energy. The choice of this function determines the driving force behind the morphological evolution of the seabed.  

Models based on the minimization principle rely on the calculation of some derivatives. This can be achieved by heavy methods (automatic differentiation) and lighter ones (analytical solutions); but they all have their drawbacks. Our strategy uses the Hadamard derivative to calculate the gradient of any cost function with respect to shape, allowing us to solve the optimization problem at the heart of the model. This strategy has enabled us to create a generic morphodynamic model that can be used with any hydrodynamic tool. Thus, our model has thus been validated both numerically (convergences, etc.) and experimentally, through flume canal experiments. 

Thanks to these developments, the code is operational in 1D and 2D and is ready to solve optimization problems linked to coastal engineering, aimed at optimizing the positions and  shapes of coastal protection structures.


This thesis was presented to the jury: 

Marissa YATES Chargée de Recherche, HDR École des Ponts, LHSV Rapporteuse
Emma Imen TURKI Maîtresse de conférences, HDR Université de Rouen Normandie Rapporteuse
Patrick MARCHESIELLO Directeur de recherche IRD, LEGOS Examinateur
Mehmet ERSOY Professeur Université de Toulon Examinateur
Frédéric BOUCHETTE Professeur Université de Montpellier Directeur de thèse
Bijan MOHAMMADI Professeur Université de Montpellier Directeur de thèse


You can find a presentation of my work here [[PDF]](http://ronan-dupont.github.io/files/presentation/planches_soutenance.pdf).



[Other information](https://www.gm.umontpellier.fr/soutenance-de-these-de-ronan-dupont-2/). 


![Editing a markdown file for a talk](/images/affiche-Soutenance_RonanDupont.jpg)
