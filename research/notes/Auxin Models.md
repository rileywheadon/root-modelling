# Auxin Models

As discussed previously, auxin efflux is primarily mediated by PINs. Almost all modelling studies involving auxin incorporate PIN in some way, with varying levels of sophistication. The simplest implementations of PIN prescribe a static gradient which replicates the distributions observed *in vivo* ([[@grieneisen2007]], [[@band2014]], [[@mellor2016]], [[@mellor2020]]). Other models have implemented dynamic PIN distributions based on cytokinin ([[@muraro2013]], [[@moore2015]], [[@dimambro2017]], [[@salvi2020]]). One notable model ([[@mironova2012]]), opted to determine the expression of different PINs based on auxin concentration, as shown below.  

![[mironova-auxin.png|center|400]]



Auxin influx is regulated by the AUX1/LAX family of transporters. Some models opt to ignore these transporters outright, instead opting for a constant auxin influx across all cells ([[@grieneisen2007]]). Other models prescribe AUX1/LAX transporters  ([[@band2014]], [[@mellor2016]], [[@mellor2020]]), which are deposited uniformly around the cell membrane. Some models have explored crosstalk between AUX1/LAX transporters and hormones such as auxin ([[@salvi2020]]) and ethylene ([[@moore2015]]).

## Hormone Crosstalk

>[!remark]+ Remark
> 
>I plan to write about hormone crosstalk models involving PLETHORA, cytokinin, and ethylene, but I haven't got around to it yet. In the meantime, here is a sketch of the existing crosstalk relationships I've found. Everything to the left of PIN has been modelled in at least one study, while the relationships between BR, CLASP, and MTs have not.

![[crosstalk.jpeg]]