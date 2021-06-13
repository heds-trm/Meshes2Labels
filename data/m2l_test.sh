#!/bin/bash

# convert meshes using reference image used for segmentation
m2l -i m2l_config.txt -r cropped_phantom_CT.mha -m -o right_pelvis_femur_label.mha 

# convert meshes using reference image used for segmentation but change spacing
m2l -i m2l_config.txt -r cropped_phantom_CT.mha -m -o right_pelvis_femur_label_1x1x1.mha --spacing 1,1,1

# create an image that fits the meshes, using a small padding of 5 mm and a spacing of 1,1,1
m2l -i m2l_config.txt -m -o right_pelvis_femur_label_fitted.mha --spacing 1,1,1 --padding 5 -f

# create an image that fits the meshes, using a small padding of 5 mm and a spacing of 1,1,1
m2l -i m2l_config.txt -m -o right_pelvis_femur_label_ofitted.mha --spacing 1,1,1 --padding 5 --ofit

# if the labels are < 256, we can create an unisgned char image
m2l -i m2l_config.txt -r cropped_phantom_CT.mha -m -o right_pelvis_femur_label_uchar.mha --uchar

# process meshes one by one, use tight fitting
m2l -i phantom_right_prox_femur.obj -l 200 --spacing 1,1,1 -f -o right_femur_label_fit.mha 

# process meshes one by one, use tight fitting
m2l -i phantom_right_prox_femur.obj -l 200 --spacing 1,1,1 --ofit -o right_femur_label_ofit.mha 

# now process the other mesh but update the label image, which will cause the hipbone to be cropped
m2l -i phantom_right_pelvis.ply -l 255 -r right_femur_label_fit.mha -u -o right_femur_pelvis_label_update.mha 

# convert meshes using reference image used for segmentation without using vtk stencil approach
m2l -i m2l_config.txt -r cropped_phantom_CT.mha -m -o right_pelvis_femur_label_1x1x1_2.mha --spacing 1,1,1 --noVtkStencil


