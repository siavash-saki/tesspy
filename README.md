# tesspy


## Introduction 
A package for urban and suburban area tessellations. 

## Installation
```python
pip install tesspy
```
## Documentation
See `tesspy.readthedocs.io`

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
ReLUT, HOLM


