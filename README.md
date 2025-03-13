# GCD-Rice_codes

## GCD-Rice: A long-term paddy rice distribution dataset in Asia at a 30 m spatial resolution

## Usage
# Instructions for use of TWDTW-SWIR1:
1. 100percent_rice_2_train.jl       -Based on the rice and non-rice samples selected from Google Earth, converted to point file(.csv).
2. download-sample-ori-TWDTW.js     -Upload the CSV file to GEE, extract the SWIR1 sequence values at the points, and download it to calculate the SWIR1 standard curve.
3. download-i-f-TWDTW.js            -Executing the download of remote sensing images in GEE.
4. LU_to_crop.py                    -Using the cropland layer of LUCC, extract the rice field layer.
5. mask_crop.py                     -Input the paddy field layer and the downloaded remote sensing image, output the remote sensing image of the paddy field, preprocessing completed.
6. area.yaml,method_info.yaml       -Modify the attached document.
7. call_distance_new.py             -Identify the pre-processed remote sensing images to obtain the rice location map.
# Instructions for use of TWDTW-SWIR1-VH:
It is generally consistent with the aforementioned TWDTW-SWIR1, but there are differences in steps 6 and 7. The details are as follows:
6. area.yaml,method_info.yaml,recognize_multibands.jl   -Modify the corresponding annex documents.
7. call_distance_multibands.py                          -Identify the pre-processed remote sensing images to obtain the rice location map.
# Instructions for use of RandomForest:
1. download-valid_num-RF.js     -Download the number of observation times for remote sensing images.
2. full-area.jl                 -Determine the extent of the imagery
3. number-distribution.py       -Statistically good observational distribution
4. apply-mask-1.jl              -Mask remote sensing image.
5. area-with-0-observation.jl   -Extract the pixels with 0 remote sensing observations per year.
6. add-crop-mask.py             -Mask the rice planting distribution map generated by the TWDTW method.
7. edge-detection.py            -Extract the edge pixels of the study area in the image
8. pick-sample.jl               -Extract four types of samples.
9. download-sample-ori-RF.js    -Download the time series for obtaining sample points.
10. train-data-x.py             -Training time series.
11. sample_curve_to_gee.py      -Merge curves
12. download-prob-RF.js         -Download the prob file from GEE.
13. apply-mask-2.py             -Mask.
14. recognize.py                -Classification and identification of rice and non-rice, including post-treatment: 1. Calculation of rice planting years, 2. Elimination of rice fields with shorter planting years, 3. Elimination of unreasonable small areas of rice fields.
# Instructions for use of Validation:
1. Based on statistical data:         area0.py,validation_together.py,compute the R², RMAE and Slope.
2. Based on the verification point:   kml-to-parts.py,point_valid.py,point_validation.py, calculate the confusion matrix.

Please refer to the paper for specific methods.
