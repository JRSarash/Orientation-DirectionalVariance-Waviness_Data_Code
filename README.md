This is the documentation of code to run 3D orientation, directional variance and waviness analysis on images.
1. Software
The code runs in the MATLAB programming language. To install, download MATLAB from: http://www.mathworks.com/products/matlab/
2. Data format
Theoretically, any format that can be read or loaded by MATLAB is suitable to de tested by the code. We include the example data in the folder “Examples”.
3. Run the analysis
Here are totally 7 MATLAB files, where ‘ParaCalc_main.m’ is the main program, and the others are functions that will be called during the running of the main program. The explanations for the variables within the code, as well as the ideas in organizing each part of the code, have been detailed in the main program. The possible parameters that should be modified accordingly to your data sets have been highlighted as ‘modify x’ (x refers to the numbering).
Generally, the main program can be divided into 4 parts:
Load images to create a 3D stack
Here the 2D images are stacked up and form a 3D stack, which is then used for the 3D orientation determination and 3D directional variance and waviness calculation.
2)	Acquire the voxel-wise 3D orientation
Here the voxel-wise 3D orientation of the 3D stack is acquired. The method is described in our previous papers (Biomed. Opt. Express 6, 2294–2310 (2015); Biomaterials, 116, 34-47 (2017); Biomaterials, 179, 96-108 (2018)).
3)	Calculate the 3D waviness and directional variance
Here the voxel-wise waviness and directional variance of collagen fibers are obtained based on the orientation information. The waviness and directional variance calculation method is described in our previous papers (Laser & Photonics Review 16 (2022), 2100576; Opt. Lett. 47 (2022) 357–360; Biomaterials, 116, 34-47 (2017); Biomaterials, 179, 96-108 (2018)).
4)	Perform post-processing
Here we generate ‘pretty’ images of orientation, waviness and directional variance. To acquire these images, the raw intensity image is used to provide the contrast of fiber features, and the orientation, waviness and directional variance maps are labeled by different colors to show the value information. 
4. Example
Images of collagen fiber as examples to test code are saved in the folder ‘Examples’. 
