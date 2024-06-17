---
title: Modeling Auxin Signaling in Roots: Auxin Computations
authors: Jaap Rutten, Thea van den Berg, Kirsten ten Tusscher
year: 2022
doi: {[DOI}}
---
# @rutten2022

## Introduction
- Auxin maxima signal the presence of a [[Stem Cell Niche]] in the root tip, which is necessary for the development of root organs ([[@sabatini1999]]).
- Asymmetries in auxin signalling direct roots away from salts (Galvan-Ampudia 2013).

>[!definition]+ Definitions
>- **Gravitropism**: The growth of plant roots and shoots either towards or away from the direction of gravity.
>- **Hydrotropism**: The growth of plant roots with respect to moisture.
>- **Phototropism**: The orientation of a plant with respect to light. 


>[!definition]+ Definitions
>- **Xylem** is vascular tissue transports water and minerals up the root. **Protoxylem** is formed first, followed by **Metaxylem**.
>- **Phloem** transports nutrients produced in photosynthesis to the plant. **Protophloem** is formed first, then **Metaphloem**.
>- **Pericycle** cells encircle the vascular tissue.


- Auxin mediates gravitropism in lateral roots (LR), allowing the root to extract nutrients from a larger area (Roychoudry 2013, Lynch 2001).
- Auxin gives subsets of protoxylem neighbouring pericycle cells the capacity for LR formation and determines LR development.
- Auxin determines protophloem differentiation rates and xylem/phloem polarization (Marhava 2018, Mahonen 2006, Bishopp 2011).

## Modelling Challenges

Auxin is involved in different processes at varying time scales:
- Rapid intracellular relocalization of [[PIN]] proteins.
- Changes in cellular gene-expression.
- Interactions with other hormone gradients.
- Cell growth (which may dilute hormone gradients)
- Cell division (which may shift gene transcription factors)

Additionally, similar auxin signals may be interpreted differently depending on the tissue (Weijers 2005, Rademacher 2011). Unpacking contextual signals is a major area of research in the field.

## Root Anatomy

>[!remark]+ List of Modelling Papers
> 
>Longitudinal Models:
>- 1D ([[@band2012]], [[@muraro2013]], [[@moret2020]])
>- 2D Rectangular ([[@grieneisen2007]], [[@mironova2010]], [[@mironova2012]], [[@mahonen2014]], [[@tian2014]], [[@moore2015]], [[@hong2017]], [[@retzer2019]])
>- 2D Idealized ([[@vandenberg2016]], [[@dimambro2017]], [[@vandenberg2018]], [[@salvi2020]])
>- 2D Anatomical ([[@band2014]], [[@mellor2016]], [[@mellor2020]], [[@muraro2016]], [[@moore2017]], [[@moore2024]], [[@liu2017]])
>
>Transverse Models:
>- 2D (Péret 2013, De Rybel 2014, Muraro 2014, Muraro 2016, el-Showk 2015, Mellor 2017)
>- 1D (Fàbregas 2015, Mellor 2019, Miyashima 2019)

Findings that a wedge-shaped root tip has significant impact on auxin patterning has led most researchers to adopt realistic models of the root.

Anatomical data, while the most realistic, also introduces asymmetries. Therefore, studies of tropisms often opt for 2D idealized models instead.

## Auxin Dynamics

Sophisticated modern models have made the following considerations:
- Auxin production via TAA and YUCCA enzymes (heightened at the QC)
- Auxin degradation mediated by GH3
- Type-specific and zone-specific [[AUX1]] influx carriers
- Prescribed type and zone specific [[PIN]] efflux carriers
- Inclusion of plasmodesmata ([[@mellor2020]])

It remains unclear whether sub-cellular auxin gradients exist and are relevant for auxin patterning.

Other regulatory components of auxin models have included:
- WOX5 and [[PLETHORA]] transcription factors
- Other hormones such as gibberellin, ethylene, and cytokinin

## Cell and Tissue Growth

[[@grieneisen2007]], [[@mironova2010]], [[@band2012]], [[@mahonen2014]], [[@salvi2020]]
- These models focus on how auxin and other hormone gradients are affected by (relatively simplistic) models of tissue growth.

[[@vos2014]], [[@fozard2013]], [[@fozard2016]], [[@weise2019]]
- These are highly accurate models of growth, which focus on physical forces acting on the cells. Reconciling these models with hormone models is an area for further research.

## Results of Auxin Modelling

- Auxin gradients are built from the PIN exporter loop ([[@grieneisen2007]]), and modulated by AUX1, ethylene, cytokinin, and plasmodesmata.
- The antagonistic relationship between auxin and cytokinin has been shown to regulate root meristem size ([[@moore2015]])
- Cytokinin signalling also contributes to the existence of an experimentally observed auxin minimum at the top of the meristem ([[@dimambro2017]])
- A positive feedback loop between auxin and [[PLETHORA]] enables the development of the root meristem after germination ([[@salvi2020]])

## Modelling Direction of Growth

TBD

## Modelling Root Branching

TBD