# ================= USER SETTINGS =================

CONDITION_NAME = "Soft HMW BME"
MAIN_PATH = "PATH/TO/ND2_FILES"

SCALE_FACTORS = {
    # Example values
    # 'image_id_1': pixels_per_micron,
    # 'image_id_2': pixels_per_micron
}

OUTPUT_FILENAME = "analysis_results.csv"

# =================================================

import numpy as np
import pandas as pd
from pathlib import Path

import skimage.filters
import skimage.morphology
from skimage.measure import label, regionprops_table
from skimage.morphology import binary_erosion, disk
from scipy.ndimage import binary_fill_holes

from aicsimageio import AICSImage

# Collect ND2 files for the selected experimental condition
main_path = Path(MAIN_PATH)
file_generator = main_path.glob('*.nd2')

file_list = list(file_generator)
print(file_list)

results=[]
cnt=0

# Iterate through all ND2 files in the directory
for file_path in file_list:
    
    # Load ND2 file using AICSImage
    data=AICSImage(file_path)

    image_dapi=data.get_image_data('ZYX',C=0, T=0)
    image_actin=data.get_image_data('ZYX',C=1, T=0)
    image_CK14=data.get_image_data('ZYX',C=2, T=0)
    image_CK8=data.get_image_data('ZYX',C=3, T=0)

    # Max intensity projection for each channel
    dapi=np.max(image_dapi[:,:,:], axis=0)
    actin=np.max(image_actin[:,:,:], axis=0)
    CK14=np.max(image_CK14[:,:,:], axis=0)
    CK8=np.max(image_CK8[:,:,:], axis=0)

    # Extract image ID from filename to match scale factor
    filename_stem = file_path.stem
    image_id = filename_stem.split('-')[0]

    # Retrieve scale factor (pixels per μm) for spatial calibration
    scale_pixels_per_μm = SCALE_FACTORS.get(image_id)

    # Handle missing scale factor information
    if scale_pixels_per_μm is None:
        print(f"WARNING: No scale factor found for image ID '{image_id}'. Using default value 1.0.")
        scale_pixels_per_μm = 1.0
    else:
        print(f"For file {file_path.name}, using scale factor: {scale_pixels_per_μm}")
    
    # Convert pixel measurements to micrometers 
    μm_per_pixel = 1 / scale_pixels_per_μm 

    # Define 10 μm radius in pixels
    radius_in_pixels_float = 10 / μm_per_pixel
    radius_in_pixels_int = int(round(radius_in_pixels_float))

    # Segmentation based on actin channel
    # Otsu threshold for main organoid body
    # Triangle threshold for capturing peripheral signal
    threshold1 = skimage.filters.threshold_otsu(actin)
    threshold2=skimage.filters.thresholding.threshold_triangle(actin)

    mask1 = actin > threshold1
    mask2 = actin > threshold2

    # Morphological cleanup to obtain a solid organoid mask
    radius=2
    dilated=skimage.morphology.binary_dilation(mask1, skimage.morphology.disk(radius))
    solid_mask = binary_fill_holes(dilated)

    # Define inner regions using erosion
    # inner_mask1: small erosion for shape analysis
    # inner_mask2: erosion corresponding to 10 μm (central region)
    inner_mask1 = binary_erosion(solid_mask, disk(2))
    inner_mask2 = binary_erosion(mask2, disk(radius_in_pixels_int))

    # Ring used for morphological measurements
    clean_ring = solid_mask & ~inner_mask1

    # Define peripheral region as a 10 μm outer ring
    periphery = binary_erosion(mask2, disk(radius_in_pixels_int))
    periphery_ring = mask2 & ~periphery


    # Label regions and extract properties 
    labeled= label(clean_ring)
    
    properties = regionprops_table(label_image=labeled, intensity_image=actin, properties=('label', 'area', 'perimeter', 'mean_intensity'))

    properties = pd.DataFrame(properties)
    
    # Skip files with no detected organoid
    if properties.empty:
        print(f"WARNING: No valid region detected in file {file_path.name}. File excluded from analysis.")
        continue

    # Select the largest region as the main organoid
    max_area=properties['area'].max()
    table=properties[properties['area'] == max_area].copy()

    area_col = table['area']
    perimeter_col = table['perimeter']

    # Morphological metric: circularity
    circularity_col = (4 * np.pi * area_col) / (perimeter_col**2)
    table['circularity'] = circularity_col

    # Intensity measurements in center and periphery
    mean_intensity_CK14_center = np.mean(CK14[inner_mask2])
    mean_intensity_CK8_center = np.mean(CK8[inner_mask2])

    mean_intensity_CK14_periphery = np.mean(CK14[periphery_ring])
    mean_intensity_CK8_periphery = np.mean(CK8[periphery_ring])

    # CK8 / CK14 intensity ratios
    ratio_center=mean_intensity_CK8_center/mean_intensity_CK14_center if mean_intensity_CK14_center > 0 else 0
    ratio_periphery=mean_intensity_CK8_periphery/mean_intensity_CK14_periphery if mean_intensity_CK14_periphery > 0 else 0
    
    table['mean_intensity_CK14_center']=mean_intensity_CK14_center
    table['mean_intensity_CK8_center']=mean_intensity_CK8_center
    table['ratio_center (CK8/CK14)']=ratio_center

    table['mean_intensity_CK14_periphery']=mean_intensity_CK14_periphery
    table['mean_intensity_CK8_periphery']=mean_intensity_CK8_periphery
    table['ratio_periphery (CK8/CK14)']=ratio_periphery

    # Convert measurements from pixels to micrometers
    table.loc[:, 'perimeter_μm'] = table['perimeter'] * μm_per_pixel
    table.loc[:, 'area_μm^2'] = table['area'] * (μm_per_pixel**2)

    # Remove pixel-based measurements
    table.pop('area')
    table.pop('perimeter')

    # Store results
    table['label']=f"Soft-HMW-ΒΜΕ-{cnt}"
    cnt+=1
    results.append(table)
    
    print(table)

# Combine results from all images
final_df = pd.concat(results, ignore_index=True)
print(final_df)

mean_row = final_df.select_dtypes(include=np.number).mean()
mean_row['label'] = 'AVERAGE'
final_df_with_summary = pd.concat([final_df, pd.DataFrame([mean_row])], ignore_index=True)

print(final_df_with_summary)

final_df_with_summary.to_csv(OUTPUT_FILENAME, index=False, sep=';', decimal=',')
