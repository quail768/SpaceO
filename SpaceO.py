#!/usr/bin/env python3

import os
import sys
import time
import json
import uuid
from pathlib import Path
import subprocess
import argparse
import numpy as np
import pandas as pd
import tifffile
import zarr
import cv2

from skimage import io
from skimage.transform import rescale

from valis import registration, micro_rigid_registrar



# =========================
# Preprocessing morphology image
# =========================
def preprocess_morphology(input_path):
    """
    Loads the pyramidal morphology.ome.tif image, creates a maximum projection image and then downsamples 
    the image to match magnfication used to image the CODEX image  in this case, 40X->20X

    """
    tif =  io.imread(input_path)

    if tif.ndim == 3:  # (z, h, w)
        mip = np.max(tif , axis=0)
    elif tif.ndim == 2:  # already 2D
        mip = tif
    else:
         raise ValueError(f"Double check your morpholgy.tif. Current shape is: {tif.shape}")
     # Downsampling with rescale to account for CODEX magnification. 
     # Currently, we image our CODEX data at 20X. Xenium is imaged at 40X. 
     # Change the downsample factor to suit your CODEX microscope configuration.
    downsampled_morphology = rescale(
        mip, 
        0.5 , 
        anti_aliasing=False ,
        preserve_range=True
    ).astype(mip.dtype)

    return downsampled_morphology

# =========================
# Function for performing Registration with Valis
# =========================

#Valis documentation can be found 
def register_affine_single_slide(
    slide_src_dir: str,
    reference_slide: str,
    results_dst_dir: str,
) -> str:

    micro_reg_fraction = 0.2
    
    #processor_dict = {he_img_f: HEDeconvolution, multichannel_img_f:ChannelGetter}


    # Initialize MicroRigidRegistrar resolution image
    registrar = registration.Valis(slide_src_dir, results_dst_dir, micro_rigid_registrar_cls=micro_rigid_registrar.MicroRigidRegistrar,reference_img_f=reference_slide,align_to_reference=True)
    rigid_registrar, non_rigid_registrar, error_df = registrar.register()

    # Calculate max_non_rigid_registration_dim_px to do non-rigid registration on an image that is 25% full resolution (Valis recommendation)
    img_dims = np.array([slide_obj.slide_dimensions_wh[0] for slide_obj in registrar.slide_dict.values()])
    min_max_size = np.min([np.max(d) for d in img_dims])
    img_areas = [np.multiply(*d) for d in img_dims]
    max_img_w, max_img_h = tuple(img_dims[np.argmax(img_areas)])
    micro_reg_size = np.floor(min_max_size*micro_reg_fraction).astype(int)

    # Perform high resolution micro registration
    micro_reg, micro_error = registrar.register_micro(max_non_rigid_registration_dim_px=micro_reg_size,align_to_reference=True)

    # Save non-reference slide as as ome.tiff
    registered_slide_dst_dir = os.path.join(results_dst_dir, "Registerd_slides")
    Path(registered_slide_dst_dir).mkdir(exist_ok=True, parents=True)
    print(registrar)

    registrar.warp_and_save_slides(registered_slide_dst_dir, non_rigid=False)
    return registered_slide_dst_dir





def load_registered_tiff(
    registered_dir: str,
    reference_slide: str,
):
    """
    Load the non-reference registered OME-TIFF from a directory

    """

    registered_dir = Path(registered_dir)
    ref_stem = Path(reference_slide).stem.lower()

    # Collect tif files ensuring compatibility with tiff and tif type images
    tiff_files = [
        p for p in registered_dir.iterdir()
        if  p.is_file()
        and p.suffix.lower()  in {".tif", ".tiff"}
        or  p.name.lower().endswith((".ome.tif", ".ome.tiff"))
    ]

    if not tiff_files:
        raise   FileNotFoundError(
            f"No TIFF files found in {registered_dir}"
        )

    def is_reference(p: Path):
        return p.stem.lower()  == ref_stem  or  p.name.lower().startswith(ref_stem)

    non_ref = [p  for p  in tiff_files if not is_reference(p)]

    if len(non_ref) != 1:
        raise RuntimeError(
            f"Expected exactly 1 non-reference TIFF, found {len(non_ref)}.\n"
            f"Candidates: {[p.name for p in tiff_files]}"
            f"Make sure you only have the runner script and the CODEX .tif file in this directory"
        )

    registered_path = non_ref[0]

    image = tifffile.imread(registered_path)

    if image.ndim != 3:
        raise ValueError(
            f"{registered_path.name} has shape {image.shape}, "
            "expected (C, H, W)"
        )

    return image, str(registered_path)


# =========================
# Use the zarr file to generate masks and converting to geojson
# =========================
 
#Code from https://github.com/gabrielascui/xenium_to_qupath/tree/main was modified to create this function 

def xenium_cells_to_geojson(zarr_dir, pixel_size, downsample_factor):
    # Adjusted scaling factor for coordinate conversion
    scale = pixel_size * downsample_factor

    # Open zarr root
    z = zarr.open(zarr_dir, mode='r')

    # Loading 'cell_id' dataset
    cell_id_data = z['cell_id'][:]

    # Converts hex digits to Xenium-style characters
    def shiftCharacters(c):
        return chr(ord('a') + int(c, 16))

    # Converts NumPy types for JSON output
    def convert_to_native(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (np.integer, np.floating)):
            return obj.item()
        elif isinstance(obj, dict):
            return {k: convert_to_native(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_to_native(i) for i in obj]
        return obj

    # Load cell boundary polygon data from polygon set 0
    vertices = z['polygon_sets/0/vertices'][:]
    num_vertices = z['polygon_sets/0/num_vertices'][:]
    cell_index = z['polygon_sets/0/cell_index'][:]

    # Create GeoJSON feature set
    features = []

    for i in range(len(vertices)):
        n = int(num_vertices[i])
        if n == 0:
            continue

        raw_coords = vertices[i][: n * 2]
        coords = [
            [raw_coords[j] / scale, raw_coords[j + 1] / scale]
            for j in range(0, len(raw_coords), 2)
        ]

        # Ensure polygon is closed
        if coords[0] != coords[-1]:
            coords.append(coords[0])

        # Extract cell ID info
        cid = int(cell_id_data[cell_index[i]][0])
        instance = int(cell_id_data[cell_index[i]][1])
        hexid = format(cid, '08x')
        shifted_id = ''.join(shiftCharacters(c) for c in hexid)
        xenium_name = f"{shifted_id}-{instance}"

        # Construct GeoJSON feature
        feature = {
            "type": "Feature",
            "id": str(uuid.uuid4()),
            "geometry": {
                "type": "Polygon",
                "coordinates": [coords]
            },
            "properties": {
                "cell_id": cid,
                "name": xenium_name,
                "objectType": "cell"
            }
        }

        features.append(feature)

    # Create final geojson object
    geojson_data = {
        "type":  "FeatureCollection",
        "features":  features
    }

    # Convert to native geojson for saving
    return convert_to_native(geojson_data)


# =========================
# Create masks
# =========================
def write_masks_and_lookup(geojson_data, image_hw, output_dir):
    """
    A multistep process that prepares the mask file and a corresponding dictionary for said mask file. 
    This mask file is essential for input to mcquant. 
    The dictionary helps to assign a cell-id from the Xenium ranger experiment to a a quantified cell.
    Each mask in the geojson object as a unique number assigned to it. 
    Step 1 - involves retrieving the unique number and it's corresponding Xenium-Ranger id from the geojson and creating a dictionary.
    Step 2 - involves estimating the dimensions of the CODEX image to be segmented.
    Step 3 - involves  the creation of a mask image tif with dimensions corresponding to the CODEX image(All quantification software requires
    that the image and the mask image tif have identical dimensions). Polygons are drawn using the geojson object and incase the 
    coordinates of the masks are found to be beyond the boundary of the CODEX image, they will be excluded from the newly drawn mask image.
    The cell-ids corresponding to these cells are written to a .txt file. The number of cell ids in this file SHOULD be low if your CODEX image was properly croppped
    and the alignment files from Valis look fine. Each drawn polygon in the mask image will retain a unique number found in the geojson object.
    This means that if a cell was number 300000 in the geojson object and the next cell was 300001, incase cell 300001 was found to be beyon the 
    bounds of the CODEX image and was not drawn, cell 300002 in the geojson object does NOT become 300001 in the mask image tif.
    Step 4 - Involves the writing of the dictionary as a .csv file , the mask image as a .tif , any out of bounds cells as a .txt 
    and a .geojson file which can be imported into qpath to assess how the xenium masks look on top of a CODEX image. 
    Ideally the polygons in this geojson file should fit perfectly/near perfectly over the mask image and should fit with high accuracy over the
    DAPI signal from the CODEX image.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    h, w = image_hw

    features = geojson_data.get("features", [])
    if not features:
        raise ValueError("No features found in GeoJSON")

    # ------------------------------------
    # Step 1 — Build name → label mapping
    # ------------------------------------
    name_to_label = {}
    labels = []
    current_label = 1

    for f in features:
        name = f.get("properties", {}).get("name")
        if name and name not in name_to_label:
            name_to_label[name] = current_label
            labels.append({"label": current_label, "name": name})
            current_label += 1

    # ------------------------------------
    # Step 2 — Create mask canvas
    # ------------------------------------
    mask = np.zeros((h, w), dtype=np.int32)

    out_of_bounds_cells = []

    # ------------------------------------
    # Step 3 — Draw polygons
    # ------------------------------------
    for f in features:
        geom = f.get("geometry", {})
        if geom.get("type") != "Polygon":
            continue

        coords = np.array(geom["coordinates"][0], dtype=np.float32)

        name = f.get("properties", {}).get("name")
        label_value = name_to_label.get(name, 0)

        # Bounds check before clipping
        min_x = np.min(coords[:, 0])
        max_x = np.max(coords[:, 0])
        min_y = np.min(coords[:, 1])
        max_y = np.max(coords[:, 1])

        if min_x < 0 or min_y < 0 or max_x >= w or max_y >= h:
            out_of_bounds_cells.append(name)
            continue  # skip drawing entirely

        poly = np.round(coords).astype(np.int32).reshape((-1, 1, 2))
        cv2.fillPoly(mask, [poly], color=int(label_value))

    # ------------------------------------
    # Step 4 — Write outputs
    # ------------------------------------

    # Save mask
    masks_path = output_dir / "ImageMasks.tif"
    tifffile.imwrite(masks_path, mask)

    # Save label lookup
    pd.DataFrame(labels).to_csv(
        output_dir / "label_lookup.csv", index=False
    )

    # Save GeoJSON copy
    with open(output_dir / "masks.geojson", "w") as f:
        json.dump(geojson_data, f, indent=4, separators=(",", ":"))

    # Save out-of-bounds cells
    if out_of_bounds_cells:
        with open(output_dir / "out_of_bounds_cells.txt", "w") as f:
            for name in out_of_bounds_cells:
                f.write(f"{name}\n")

        print(
            f"WARNING: {len(out_of_bounds_cells)} cells were outside image bounds "
            f"and were excluded from the mask."
        )

    return masks_path



def run_mcquant(
    masks_path, image_path, output_dir, channel_names_csv
):
    """
    Use mcquant from the mcmicro package to quantify the channels. Automatically uses  masks(polygons) from the 
    polygons that were drawn from the geojson. Channels should be be a csv file with no header and each channel on a separate line
    """
    cmd = [
        "mcquant",
        "--masks", str(masks_path),
        "--image", str(image_path),
        "--output", str(output_dir),
        "--channel_names", str(channel_names_csv),
    ]

    subprocess.run(cmd, check=True)



def get_non_lookup_csv(directory):
    """
    Returns the CSV file in the directory that is NOT named 'label_lookup.csv'.
    Assumes there is exactly one such CSV.
    """
    directory = Path(directory)

    # Find all CSV files
    csv_files = list(directory.glob("*.csv"))

    # Filter out label_lookup.csv
    non_lookup = [f for f in csv_files if f.name != "label_lookup.csv"]

    if len(non_lookup) != 1:
        raise ValueError(f"Expected exactly one non-lookup CSV, found {len(non_lookup)}")

    return non_lookup[0] 

def merge_codex_with_labels(mask_csv, lookup_csv):
    """
    Reads two CSV files and performs a left join:
    mask_data.CellID == label_lookup.label
    """
    # Load the tables
    mask_df = pd.read_csv(mask_csv)
    lookup_df = pd.read_csv(lookup_csv)

    # Perform the left join
    merged = mask_df.merge(
        lookup_df,
        how="left",
        left_on="CellID",
        right_on="label"
    )

    return merged

# Create the ArgumentParser object
parser = argparse.ArgumentParser(
    description="Get arguments for SpaceO")
parser.add_argument(
    "Xenium_Image", help="Path to morphology.ome.tif from Xenium ranger output ")
parser.add_argument(
    "Phenocycler_Directory_Home", help="Path to a directory containing the Phenocyler image")
parser.add_argument(
    "Zarr", help="Path to a directory containing the unzipped cells.zarr file from the Xenium ranger output")
parser.add_argument(
    "Channels", help="Path to a .csv file with each row having the name of the channel used for imaging")
parser.add_argument(
    "Spatial_Proteomics_Magnification", help="The magnification used for imaging spatial proteomics data")

# Parse the arguments 
args = parser.parse_args()

# =========================   
# =========================
# MAIN
# =========================
def main():


    input_path = args.Xenium_Image
    SpatialProteomicsImageDirectory = args.Phenocycler_Directory_Home
    zarr_dir = args.Zarr
    channels_csv = args.Channels

    RegistrationOutput = Path(SpatialProteomicsImageDirectory) /"Results"
    reference_slide = Path(SpatialProteomicsImageDirectory) / "Xenium_DAPI.tif"
    Quantification_Output = Path(SpatialProteomicsImageDirectory) / "Quantification"
    lookup_csv = Path(Quantification_Output) / "label_lookup.csv"
    pixel_size = 0.2125
    downsample_factor = 40/args.Spatial_Proteomics_Magnification

    start = time.time()

    print("Preprocessing morphology...")
    ref_img = preprocess_morphology(input_path)
    img16 = ref_img.astype(np.uint16)
    io.imsave(Path(SpatialProteomicsImageDirectory) / "Xenium_DAPI.tif", img16)
                
    print("Registering...")
    registered_path = register_affine_single_slide(
        SpatialProteomicsImageDirectory, reference_slide, RegistrationOutput
    )

    print("Loading registered image...")
    image, registered_path = load_registered_tiff(
        registered_dir= registered_path,
        reference_slide=reference_slide,
    )

    _, h, w = image.shape

    print("Generating GeoJSON...")
    geojson_data = xenium_cells_to_geojson(
        zarr_dir, pixel_size, downsample_factor
    )
    
    print("Rasterizing masks...")
    
    masks_path = write_masks_and_lookup(
        geojson_data, (h, w), Quantification_Output
    )

    print("Running MCQuant...")
    run_mcquant(
        masks_path,
        registered_path,
        Quantification_Output,
        channels_csv,
    )

    print("Adding labels...")
    quantified_unamed_cells=get_non_lookup_csv(Quantification_Output)
    CountMatrixwithLabels = merge_codex_with_labels(
    quantified_unamed_cells,lookup_csv)
    pd.DataFrame(CountMatrixwithLabels).to_csv(
        Quantification_Output / "CountMatrixwithLabels.csv", index=False
    )
    print("Done!")
    print(f"Time: {time.time() - start:.2f}s")



if __name__ == "__main__":
    main()
