- [Learning to use MPDAF](#org915755a)
    - [00 MPDAF simple demo](#org988d12c)
  - [N66/NGC 346 in the SMC](#orgcec9bf4)
    - [Velocity maps](#org6ccc1cb)

-   New project <span class="timestamp-wrapper"><span class="timestamp">[2021-05-04 Tue]</span></span>
-   Two principal aims:
    1.  Analyzing velocity statistics, such as structure functions
        -   Collaboration with Javier García Vázquez
    2.  Looking for Raman-scattered Balmer wings
        -   Collaboration with Mabel Valerdi


<a id="org915755a"></a>

# Learning to use MPDAF

I have made a series of notebooks


<a id="org988d12c"></a>

## 00 MPDAF simple demo

-   Notebook: [notebooks/00 MPDAF simple demo.ipynb](notebooks/00 MPDAF simple demo.ipynb)
-   Pure python: [notebooks/00 MPDAF simple demo.py](notebooks/00 MPDAF simple demo.py)
-   Reading in the data cube, summing all waves to get an image, summing all pixels to get a spectrum


<a id="orgcec9bf4"></a>

# N66/NGC 346 in the SMC

-   N66 is the H II region, while NGC 346 is technically the star cluster, but the H II region is also sometimes referred to as this
-   Data are at <file:///Users/will/Work/Muse-Hii-Data/SMC-NGC-346/ADP.2017-10-16T11_04_19.247.fits>
-   Papers on NGC 346
    -   Valerdi:2019a Helium abundance
        -   Measures physical conditions too:
            -   [O III] Te = 13000 K
            -   ne = 20-30 ([O II] and [S II]) up to 100 ([Fe III])


<a id="org6ccc1cb"></a>

## Velocity maps

-   Looks like this will be possible for the brighter lines at least
-   There is a problem that seems to be an overzealous background subtraction, which means that some of the lines come out negative
    -   We can maybe fix this by adding on the spectrum seen in the darkest corner of the image
