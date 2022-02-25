# tesspy


## Introduction 
A package for urban and suburban area tessellations. 

## Installation
You can install ``tesspy`` from PyPI using pip:
```
pip install tesspy
```

and from conda:
```
conda install -c conda-forge tesspy
```



## Creating a new environment for tesspy

`tesspy` depends on `geopandas`, which could make the installation sometimes tricky because of the conflicts with the current packages. Therefore, we recommend creating a new clean environment and installing the dependencies from the conda-forge channel.


Create a new environment:
```shell
conda create -n tesspy_env -c conda-forge
```

Activate this environment:
```shell
conda activate tesspy_env
```

Install tesspy from conda-forge channel:
```shell
conda install -c conda-forge tesspy
```


## Documentation
The official documentation is hosted on **[ReadTheDocs](https://tesspy.readthedocs.io)**

## Examples
Example for a germany city: Frankfurt am Main, an US-city: Key West and the capital of Nairobi. In all examples the local projection was used such that any kind of distortion is minimized. The following picutres shows the boundary polygons. 
![This is an image](Examples/pics/Boundary_polygons.png)

### Squares 
For square based tessellations different resolution were used the get similar looking results. 
![This is an image](Examples/pics/Squares.png)

### Hexagons
For hexagona based tessellations different resolution were used the get similar looking results. 
![This is an image](Examples/pics/Hexagons.png)

### Adaptive Squares
For the adaptive square based tessellations different start resolutions were used. As POI the feature "amenity" was used.
![This is an image](Examples/pics/Adaptive_Squares.png)

### Voronoi Polygons
The Voronoi-Diagram based tessellations used the "building" feature to cluster those data and create the generators.
![This is an image](Examples/pics/Voronoi.png)

### City Blocks
For City Blocks all highway-types were used
![This is an image](Examples/pics/CB.png)

## Contributing to tesspy
Contribution is welcome!

## Acknowledgements
`tesspy` is the result of the research project [ClusterMobil](https://www.frankfurt-university.de/de/hochschule/fachbereich-1-architektur-bauingenieurwesen-geomatik/forschungsinstitut-ffin/fachgruppen-des-ffin/fg-neue-mobilitat/relut/forschungsprojekte-relut/clustermobil/) conducted by the [Research Lab for Urban Transport](https://www.frankfurt-university.de/en/about-us/faculty-1-architecture-civil-engineering-geomatics/research-institute-ffin/specialist-groups-of-the-ffin/specialist-group-new-mobility/relut/). This research project is funded by the state of Hesse and [HOLM](https://frankfurt-holm.de/) funding under the “Innovations in Logistics and Mobility” measure of the Hessian Ministry of Economics, Energy, Transport and Housing. [HA Project No.: 1017/21-19]

