# Refocus Tarantula paper on SW cloud

This is the region of the nebula where we can see a very clean progression from the highly ionized lines to the neutral/molecular lines. It is a good analog to the Orion Bar, where we see the PDR at an approximately edge-on orientation.

![CleanShot 2022-06-23 at 13.37.06.jpeg](assets/CleanShot%202022-06-23%20at%2013.37.06.jpeg)

This shows the extremes of ionization, from the highly ionized [Cl IV] line in blue to the [C I] line that comes form the CO dissociation front in red.

## Extinction

We already looked at this in one of the notebooks.  There are two aspects that particularly concern us here:

1. How much of the structure of the high-ionization part is caused by extinction? Is the inner edge a real feature?
   - I *feel that the answers are not much and yes, respectively.*
1. Is there any hidden ionized gas behind the molecular cloud in the bottom right?  If so, does this invalidated the idea of it being an edge-on geometry.
   - *It seems like there is, given that the radio continuum brightness has excess over Paschen alpha in that zone*

> **Dust x-section per H atom in LMC**

> The metallicity is roughly half solar

> Bernard et al 2008 find that the dust/gas ratio is 4.6 to 2.3 times lower than in the MW. This is derived from the FIR surface brightness, so the translation to optical extinction would require finding the extinction curve.  But, if we just take 3 as a typical value and then scale the Orion value, we would get **1.7e-22 cm2/H**

$$
A_V / N_\mathrm{H} = 1.7 \times 10^{-22} \ \mathrm{cm^2/H}
$$

## Geometry and densities

> At a distance of 50 kpc: **1 arcsec = 0.242 parsec**

### To what extent is the IF/PDR seen edge-on?

The CO cloud is quite compact in azimuth (PA = 250 +/- 5) , whereas the ionized arc is much broader (PA = 240 +/- 30), so it doesn’t seem feasible that the cloud is directly responsible for confining the arc. On the other hand, the Raman H I emission *is* more laterally extended, similar to the ionized arc, so maybe the lower density neutral gas is sufficient

The appearance is perhaps consistent with the cloud being seen from behind and obscuring some of the ionized gas behind it. There is also the evidence of the radio and IR emission, which shows that there is an absorption hole in the Br gamma, which is not there in the 21 cm continuum.

### Estimate of column densities

The thickness of the Raman layer is about 5 arcsec, which is about 1.2 parsec.

- Assume that this corresponds to visual τ = 1 and that the dust extinction cross section per H scales with the metallicity, making it be about 2e-22 cm2/H.  *Above I calculated 1.7 instead of 2, but that will not make much difference.*
- So this means a neutral column of 5e21 cm-2
- So average neutral density is 5e21 / 1.2 3.08e18 = **1350 pcc**
   - This could be compared with estimates for the ionized density and the molecular density
   - From the figure in Castro 2018 it looks like the [S II] density is about 500 pcc, which would mean a rather low contrast between the ionized and the neutral density
   - We could also make an estimate from the Ha surface brightness, although that would require estimating the line of sight thickness
   - And we could also do a census of the ionizing sources and estimate the flux incident on the ionized filament, and thus what the emission measure should be

## Radio continuum

In principal could get the absolute extinction from comparing the optical/ir hydrogen lines with the radio continuum.

I got a map at 21 cm from CASDA

### Lazendic et al 2003 radio maps and extinction measurements

Find a supposed supernova remnant at the same place that I had identified anomalous radio emission in the 21 cm maps

> This makes me wonder whether the similar excess in Orion may also be a SNR

![CleanShot 2022-07-04 at 11.26.57.jpeg](assets/CleanShot%202022-07-04%20at%2011.26.57.jpeg)

![CleanShot 2022-07-04 at 11.28.48.jpeg](assets/CleanShot%202022-07-04%20at%2011.28.48.jpeg)

On the left is the 3 cm map from Lazendic. On the right is the 21 cm map that I got.  The green circle is the “source” **MCRX J053831.8–690620** from Lazendic. They think that this is an SNR rather than high extinction but this seems unlikely for the following reasons:

## Kinematics

It looks like there is quite a simple velocity structure

![CleanShot 2022-06-25 at 00.11.54.jpeg](assets/CleanShot%202022-06-25%20at%2000.11.54.jpeg)

![CleanShot 2022-06-25 at 00.13.26.jpeg](assets/CleanShot%202022-06-25%20at%2000.13.26.jpeg)

![CleanShot 2022-06-25 at 00.15.47.jpeg](assets/CleanShot%202022-06-25%20at%2000.15.47.jpeg)

![CleanShot 2022-06-25 at 00.17.24.jpeg](assets/CleanShot%202022-06-25%20at%2000.17.24.jpeg)

Range is 210 to 310 heliocentric in each case.

![CleanShot 2022-06-25 at 00.19.50.jpeg](assets/CleanShot%202022-06-25%20at%2000.19.50.jpeg)

![CleanShot 2022-06-25 at 00.25.36.jpeg](assets/CleanShot%202022-06-25%20at%2000.25.36.jpeg)

## X rays

I am not sure we are going to get much useful out of this, since the resolution is relatively poor

### eROSITA observations

[Sasaki et al 2022](https://ui.adsabs.harvard.edu/abs/2022A&A...661A..37S )

![CleanShot 2022-06-23 at 21.03.42.jpeg](assets/CleanShot%202022-06-23%20at%2021.03.42.jpeg)

![CleanShot 2022-06-23 at 21.07.24.jpeg](assets/CleanShot%202022-06-23%20at%2021.07.24.jpeg)

The region that they highlight as **30 Dor Center** is precisely our region.

![CleanShot 2022-06-23 at 21.09.36.jpeg](assets/CleanShot%202022-06-23%20at%2021.09.36.jpeg)

This is the spectrum, which is quite soft, but with a high foreground extinction.

![CleanShot 2022-06-23 at 21.11.20.jpeg](assets/CleanShot%202022-06-23%20at%2021.11.20.jpeg)

These are the fit parameters.  They have two thermal components with 2 million K and 4 million K, approximately. The foreground extinction of 5.5e21 corresponds to about AV = 3, which seems somewhat high.

### Comparison of eROSITA, Suzaku, and Chandra

First the eROSITA image from Sasaki et al (2022), see above.

![CleanShot 2022-06-24 at 08.44.54.jpeg](assets/CleanShot%202022-06-24%20at%2008.44.54.jpeg)

Next the Chandra image (right panel).  Original paper was Townsley+ (2006) but that does not have very good figures. This version is from Cheng et al (2021).

![CleanShot 2022-06-24 at 08.46.19.jpeg](assets/CleanShot%202022-06-24%20at%2008.46.19.jpeg)

## Infrared observations and associated modeling

### Herschel/Spitzer maps from Chevance et al 2016

![CleanShot 2022-07-03 at 20.54.15.jpeg](assets/CleanShot%202022-07-03%20at%2020.54.15.jpeg)

![CleanShot 2022-07-03 at 20.54.43.jpeg](assets/CleanShot%202022-07-03%20at%2020.54.43.jpeg)

This map shows the ionization stratification from highly ionized [S IV] to PDR [C II].

The top-left part of Region G corresponds to the bottom right corner of our map.

### Spitzer maps from Indebetouw et al 2009

![CleanShot 2022-07-03 at 22.37.27.jpeg](assets/CleanShot%202022-07-03%20at%2022.37.27.jpeg)

This is a similar area to above. I have outlined the approximate MUSE field in black.

![CleanShot 2022-07-03 at 22.47.24.jpeg](assets/CleanShot%202022-07-03%20at%2022.47.24.jpeg)

![CleanShot 2022-07-03 at 22.55.40.jpeg](assets/CleanShot%202022-07-03%20at%2022.55.40.jpeg)

On the left is the mid IR [S IV]/[S III] ratio. On the right is the MUSE [Ar IV]/[O III] ratio, which has a very similar distribution, but is mapped at higher resolution.

