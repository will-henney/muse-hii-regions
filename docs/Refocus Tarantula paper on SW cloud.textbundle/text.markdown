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

> **3-color H line image for 2.16, 0.92, 0.49 micron, which nicely shows the effects of extinction.**

![CleanShot 2022-07-04 at 15.08.45.jpeg](assets/CleanShot%202022-07-04%20at%2015.08.45.jpeg)

> The blue-white bits have the lowest extinction. Yellow is moderate extinction that shows up in the optical, while redder parts are higher extinction that saturate in the optical but show up in the infrared.  These mainly show up as reddish brown filaments. Some of them are absorption filaments in both optical and NIR, such as the ones that cross the center of the star cluster. But some of them are *emission* filaments, such as the one that passes directly under the molecular globule (marked *NIR filament*) in the figure. l

> The black contours are the 21cm radio continuum emission, while the white contours are 13CO.

> In general, the radio free-free continuum follows the H lines, but right on top of the CO emission there is a radio peak, which is completely invisible in the H lines, even at 2 micron.  It is right behind where the NIR filament crosses the CO filament.

> Also, we can see a protrusion of the CO globule, marked *overlap*, which seems to produce extinction of the main ionized shell.  This corresponds to a region where the radio emission extends to a slightly greeter radius than the optical/ir, which is consistent with overlapping extinction. For the main globule though, at the point where our cut crosses it, which is marked *no overlap*, it seems that the ionization front is seen roughly edge-on since the radio emission falls just as fast as the optical emission (although we really need to get our hands on the 3 cm maps to be sure, since the 20 cm map is too low resolution).

> **Dust x-section per H atom in LMC**

> The metallicity is roughly half solar

> Bernard et al 2008 find that the dust/gas ratio is 4.6 to 2.3 times lower than in the MW. This is derived from the FIR surface brightness, so the translation to optical extinction would require finding the extinction curve.  But, if we just take 3 as a typical value and then scale the Orion value, we would get **1.7e-22 cm2/H**

$$
A_V / N_\mathrm{H} = 1.7 \times 10^{-22} \ \mathrm{cm^2/H}
$$

## Geometry and densities

> At a distance to LMC of 50 kpc:

> **1 arcsec = 0.242 parsec**

> **1 arcmin = 14.5 parsec** (size of single MUSE field)

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

> This makes me wonder whether the similar excess in Orion may also be a SNR. **NO, this is very unlikely**

![CleanShot 2022-07-04 at 11.26.57.jpeg](assets/CleanShot%202022-07-04%20at%2011.26.57.jpeg)

![CleanShot 2022-07-04 at 11.28.48.jpeg](assets/CleanShot%202022-07-04%20at%2011.28.48.jpeg)

On the left is the 3 cm map from Lazendic. On the right is the 21 cm map that I got.  The green circle is the “source” **MCRX J053831.8–690620** from Lazendic. They think that this is an SNR rather than high extinction but this seems unlikely for the following reasons:

1. The association with the CO cloud.
2. The spectral index between 3 cm and 21 cm looks very flat. At least, from comparing the above two maps, the contrast between this source and the diffuse ionized emission is very similar at the two wavelengths.
3. It coincides with a feature in the 2 micron to 1 micron reddening, see below.

![CleanShot 2022-07-04 at 12.28.15.jpeg](assets/CleanShot%202022-07-04%20at%2012.28.15.jpeg)

![CleanShot 2022-07-04 at 12.56.18.jpeg](assets/CleanShot%202022-07-04%20at%2012.56.18.jpeg)

![CleanShot 2022-07-04 at 13.30.20.jpeg](assets/CleanShot%202022-07-04%20at%2013.30.20.jpeg)

![CleanShot 2022-07-04 at 12.54.29.jpeg](assets/CleanShot%202022-07-04%20at%2012.54.29.jpeg)

The first panel is from the Lazendic paper and shows the radio-optical extinction. The remaining panels are from one of my notebooks and show the optical-infrared reddening.  The supposed source coincides with high reddening between 2.15 and 0.92 microns, but nothing at shorter wavelengths. This suggests very high extinction so it is not seen at visual wavelengths.

In fact, it turns out that this was already known.  The following paragraph is from Indebetouw et al (2009):

![CleanShot 2022-07-05 at 12.56.55.jpeg](assets/CleanShot%202022-07-05%20at%2012.56.55.jpeg)

Chu et al (2004) look at it in detail and conclude that it cannot be a SNR.

## Kinematics

It looks like there is quite a simple velocity structure

![CleanShot 2022-06-25 at 00.11.54.jpeg](assets/CleanShot%202022-06-25%20at%2000.11.54.jpeg)

![CleanShot 2022-06-25 at 00.13.26.jpeg](assets/CleanShot%202022-06-25%20at%2000.13.26.jpeg)

![CleanShot 2022-06-25 at 00.15.47.jpeg](assets/CleanShot%202022-06-25%20at%2000.15.47.jpeg)

![CleanShot 2022-06-25 at 00.17.24.jpeg](assets/CleanShot%202022-06-25%20at%2000.17.24.jpeg)

Range is 210 to 310 heliocentric in each case.

![CleanShot 2022-06-25 at 00.19.50.jpeg](assets/CleanShot%202022-06-25%20at%2000.19.50.jpeg)

![CleanShot 2022-06-25 at 00.25.36.jpeg](assets/CleanShot%202022-06-25%20at%2000.25.36.jpeg)

## Scattering by dust

I can think of two techniques for looking at this:

1. The visual continuum after subtracting a scaled version of the hydrogen lines to account for the bound-free.  The best thing here would be to get the continuum just redward of the Paschen jump. Maybe around 8500 Å
2. The He II lines. Specifically 4686 looks like it has 3 parts to it:
   1. The emission from the WR stars and early O stars, which is far brighter than the diffuse emission by a factor of 1000
   2. The intrinsic diffuse emission, which is associate with the ArIV and other high ionization lines
   3. The scattered diffuse emission, which seems to trace the PDRs. ***Possibly these last two could be separated better by trying to isolate the narrow and broad components of the line.***

## The compact inner globule

This is a small cometary globule that is best seen in the the Br-gamma, [O I], and 3 cm radio emission.  It is quite heavily extinguished in the optical lines, suggesting it is foreground to the cluster and illuminated from behind.

![CleanShot 2022-07-05 at 08.49.39.jpeg](assets/CleanShot%202022-07-05%20at%2008.49.39.jpeg)

It turns out that the globule is associated with a YSO supposedly in the catalog of Gruendl & Chu (2009).  But it does not make the top-ten list for YSO sources in Walborn et al (2013). It is source **118** in Rubio et al (1998), which SIMBAD has as [`[RBW98]`](https://cds.unistra.fr/cgi-bin/Dic-Simbad?%5bRBW98%5d) `IRSW-118` associates with the O9.5 star `[P93] 702` (see Parker 1993), which is `VFTS 464` in the *VLT-FLAMES Tarantula Survey* (this is the number used in all the Walborn papers, for instance) and has designation `IRSF J05383928-6905527` in the InfraRed Survey Facility Catalog (Kato et al 2007). The *V*=16.8 magnitude is too faint to make Parker’s cut-off (*V* < 16) for the bright blue stars in his Table 11.

![CleanShot 2022-07-05 at 20.53.38.jpeg](assets/CleanShot%202022-07-05%20at%2020.53.38.jpeg)

The above table is from the first VFTS paper (Evans et al 2011).  Note that the offset from the GC09 position is excessively large (1.2 arcsec), which would make me suspect that the star is not the same as the IR source. **However, I have compared the positions myself and find the distance to be only 0.5 arcsec (0.13 pc).** The star is just in front of the globule, with the IR YSO candidate position corresponding to the center of the CO emission. The star position coincides approximately with the peak of the Br-gamma emission, so maybe some of that is from the star (although we would not expect an O9 star to have a H line in emission unless it was a supergiant, so I do not actually believe this).

![CleanShot 2022-07-05 at 22.18.44.jpeg](assets/CleanShot%202022-07-05%20at%2022.18.44.jpeg)

![CleanShot 2022-07-05 at 22.50.19.jpeg](assets/CleanShot%202022-07-05%20at%2022.50.19.jpeg)

This is from Campbell et al (2010) who did NIR adaptive optics imaging, so much higher resolution than we have.  This shows that the stellar source is in fact double, and it is indeed just in front of the globule. It is most likely that the star is foreground to the globule and is not the dominant ionizing source (I should work out the total equivalent Ha luminosity from it). They comment on it being a bow-shock but of course it is not a bow shock, it is clearly an ionization front. In fact, if I had carried on reading the Evans et al (2011) paper, I would have found the figure shown at the right.  In that paper, they already realize that it is not a bow shock: *“Indeed, as noted by Campbell et al. (2010), the bow shock is orientated towards R136 (and Brey 75/BAT99-100), suggesting it might well be related to an ionization front (with associated triggered star-formation) rather than a dynamical shock.”*

More recent observations with adaptive optics of the core of the cluster have a resolution of 12.25 mas/pixel (Khorrami et al 2021). But I am not sure if they go out as far as our source – *I need to check ….*

- Distance from R134 to globule: 8.4 arcsec = 2 pc
- Distance from R136 to globule: 20.5 arcsec = 5 pc
- Distance from R136 to bar i-front: 48 arcsec = 12 pc

The *Stapler Clouds* are named and described in Kalari et al (2018). They have them as being 10 to 20 pc in front of the cluster, but that seems a bit too close given that they have no strong radio continuum emission.

![CleanShot 2022-07-05 at 12.22.53.jpeg](assets/CleanShot%202022-07-05%20at%2012.22.53.jpeg)

### Other globules in the SW region

![CleanShot 2022-07-05 at 22.56.22.jpeg](assets/CleanShot%202022-07-05%20at%2022.56.22.jpeg)

These are bright in low ionization lines, but they have a diffuse peak of high ionization emission associated with them too

![CleanShot 2022-07-05 at 22.59.50.jpeg](assets/CleanShot%202022-07-05%20at%2022.59.50.jpeg)

Velocities are similar to the bar in O I, and also in CO

![CleanShot 2022-07-05 at 23.11.39.jpeg](assets/CleanShot%202022-07-05%20at%2023.11.39.jpeg)

## Stellar census

This Table is from Crowther et al (2016).  In R136 there is a 1e7 Lsun star, plus 8 others that are on average 3e6, so that is getting on for 4e7, or 200 times luminosity of Orion.  The factor in QH will be even larger, since these are hotter, say 300, which would make it 3e51.  In the Javier paper, we have L(Ha) being 200 times higher, so that just about checks out.

![CleanShot 2022-07-05 at 23.35.38.jpeg](assets/CleanShot%202022-07-05%20at%2023.35.38.jpeg)

So the globule being 5 pc away from the cluster would be the same ionizing flux as 0.3 pc away from the Trapezium. This is only slightly farther than the distance to the Orion Bar.  Radius of globule is about 1 arcsec, so 0.24 pc, as opposed to radius of Orion Bar, which is about 10 arcsec, or 0.02 pc, which is 10 times smaller.  The ionized density should go as `sqrt(F/R)`, so that will be 3 times smaller in the Tarantula globule than in Orion Bar.  At the same time, the ionization parameter will  be 3 times higher.

The W filament region is 2.4 times further away, which would be equivalent of 0.7 pc away in Orion, which is about 360 arcsec, which is roughly LL2 distance, I think.  So that is a 10 times smaller ionizing flux than the Orion Bar. The thickness of the ionized filament is about 12 arcsec, or 3 pc, which is 150 times larger than Orion Bar. So density should be `sqrt(1500)=40` times smaller. This would mean about 100 pcc.

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

Their modeling strategy is to use plane parallel static constant pressure PDR models (assuming face-on illumination and viewing angle I guess) and to derive a pressure *P* and a radiation field (*G*_FUV) from line ratios of [O I] and [C II] lines and the dust continuum (all in the FIR). Then comparing predicted/observed absolute fluxes gives an area filling factor Phi_A, which if smaller than one means that the individual clouds are smaller than the beam, or if bigger than one indicate crowding along the line of sight, which is sort-of-like limb brightening.

There is a further parameter, which is the maximum depth of the model AV(max), which they find has a value of 1 to 3 from comparison with CO observations.  Then they explain the higher observed extinctions as being due to the crowding effect:

![CleanShot 2022-07-05 at 13.49.35.jpeg](assets/CleanShot%202022-07-05%20at%2013.49.35.jpeg)

This is possibly reasonable when considering spatial averages over 15 arcsec, which is the resolution of their observations. I am not sure it makes sense at smaller scale though, since the CO is concentrated in small clouds, whereas the neutral gas isn’t . On the other hand, there are only a few places where the extinction is super-high

### Spitzer maps from Indebetouw et al 2009

![CleanShot 2022-07-03 at 22.37.27.jpeg](assets/CleanShot%202022-07-03%20at%2022.37.27.jpeg)

This is a similar area to above. I have outlined the approximate MUSE field in black.

![CleanShot 2022-07-03 at 22.47.24.jpeg](assets/CleanShot%202022-07-03%20at%2022.47.24.jpeg)

![CleanShot 2022-07-03 at 22.55.40.jpeg](assets/CleanShot%202022-07-03%20at%2022.55.40.jpeg)

On the left is the mid IR [S IV]/[S III] ratio. On the right is the MUSE [Ar IV]/[O III] ratio, which has a very similar distribution, but is mapped at higher resolution.

